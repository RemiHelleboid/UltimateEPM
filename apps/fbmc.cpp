/**
 * @file electron_phonon.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-02-10
 *
 * @copyright Copyright (c) 2024
 *
 */

/**
 * @file electron_phonon.cpp
 * @brief Single-particle FBMC driver for electron-phonon transport.
 */

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <tclap/CmdLine.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "BandStructure.h"
#include "Options.h"
#include "bz_mesh.hpp"
#include "bz_states.hpp"
#include "electron_phonon.hpp"
#include "epm_material.hpp"
#include "single_part_fbmc.hpp"

namespace {

void require_finite(double value, const std::string& option_name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(fmt::format("{} must be finite", option_name));
    }
}

void require_positive(double value, const std::string& option_name) {
    require_finite(value, option_name);
    if (!(value > 0.0)) {
        throw std::invalid_argument(fmt::format("{} must be positive", option_name));
    }
}

void require_positive_int(int value, const std::string& option_name) {
    if (value <= 0) {
        throw std::invalid_argument(fmt::format("{} must be positive", option_name));
    }
}

void require_existing_file(const std::filesystem::path& path, const std::string& what) {
    if (!std::filesystem::exists(path)) {
        throw std::runtime_error(fmt::format("{} does not exist: {}", what, path.string()));
    }
    if (!std::filesystem::is_regular_file(path)) {
        throw std::runtime_error(fmt::format("{} is not a regular file: {}", what, path.string()));
    }
}

bool looks_like_phonon_rates_csv(const std::filesystem::path& path) {
    if (!std::filesystem::is_regular_file(path) || path.extension() != ".csv") {
        return false;
    }

    std::ifstream stream(path);
    std::string   header;
    if (!std::getline(stream, header)) {
        return false;
    }

    return header.find("vertex_index") != std::string::npos && header.find("local_band_index") != std::string::npos &&
           header.find("energy_eV") != std::string::npos && header.find("rate_ac_L_ab") != std::string::npos &&
           header.find("rate_op_T_em") != std::string::npos;
}

std::filesystem::path detect_phonon_rates_file(const std::filesystem::path& directory) {
    std::vector<std::filesystem::path> matches;
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (looks_like_phonon_rates_csv(entry.path())) {
            matches.push_back(entry.path());
        }
    }

    std::sort(matches.begin(), matches.end());

    if (matches.empty()) {
        throw std::runtime_error(
            fmt::format("No phonon-rate CSV was provided and no matching file was found in {}. "
                        "Pass --phononfile explicitly or generate one with elph.epm --export-rates.",
                        directory.string()));
    }

    if (matches.size() > 1) {
        std::ostringstream message;
        message << "Multiple phonon-rate CSV files were found in " << directory.string()
                << ". Pass --phononfile explicitly. Matches:";
        for (const auto& match : matches) {
            message << "\n  " << match.string();
        }
        throw std::runtime_error(message.str());
    }

    return matches.front();
}

std::filesystem::path make_output_directory(const std::string&            requested,
                                            const std::string&            material,
                                            double                        temperature,
                                            const uepm::mesh_bz::vector3& electric_field_v_per_cm) {
    std::filesystem::path outdir;
    if (requested.empty()) {
        if (electric_field_v_per_cm.y() == 0.0 && electric_field_v_per_cm.z() == 0.0) {
            outdir = fmt::format("fbmc_{}_{:.1f}K_{:.3e}Vcm", material, temperature, electric_field_v_per_cm.x());
        } else {
            outdir = fmt::format("fbmc_{}_{:.1f}K_E_{:.3e}_{:.3e}_{:.3e}Vcm",
                                 material,
                                 temperature,
                                 electric_field_v_per_cm.x(),
                                 electric_field_v_per_cm.y(),
                                 electric_field_v_per_cm.z());
        }
    } else {
        outdir = requested;
    }

    if (!std::filesystem::exists(outdir)) {
        std::filesystem::create_directories(outdir);
    }

    if (!std::filesystem::is_directory(outdir)) {
        throw std::runtime_error(fmt::format("Output path is not a directory: {}", outdir.string()));
    }

    return outdir;
}

uepm::fbmc::particle_type parse_particle_type(const std::string& value) {
    if (value == "electron") {
        return uepm::fbmc::particle_type::electron;
    }
    if (value == "hole") {
        return uepm::fbmc::particle_type::hole;
    }
    throw std::invalid_argument("--carrier must be 'electron' or 'hole'");
}

void write_run_metadata(const std::filesystem::path&               outdir,
                        const std::string&                         mesh_file,
                        const std::string&                         rates_file,
                        const std::string&                         material,
                        const std::string&                         phonon_parameter_set,
                        const std::string&                         phonon_parameter_file,
                        const std::string&                         impact_ionization_parameter_set,
                        int                                        nb_particles,
                        int                                        nb_threads,
                        int                                        nb_conduction_bands,
                        int                                        nb_valence_bands,
                        double                                     simulation_time,
                        double                                     temperature,
                        const uepm::mesh_bz::vector3&              electric_field_v_per_cm,
                        double                                     max_energy_eV,
                        double                                     gamma_safety,
                        double                                     warmup_fraction,
                        std::string_view                           bz_domain,
                        std::string_view                           carrier,
                        std::uint64_t                              random_seed,
                        bool                                       enable_impact_ionization,
                        const uepm::fbmc::KeldyshImpactIonization& impact_ionization_model,
                        bool                                       export_history) {
    const auto    meta_file = outdir / "run_info.txt";
    std::ofstream os(meta_file);
    if (!os) {
        throw std::runtime_error(fmt::format("Could not open metadata file for writing: {}", meta_file.string()));
    }

    os << "material = " << material << '\n';
    os << "phonon_parameter_set = " << phonon_parameter_set << '\n';
    os << "phonon_parameter_file = " << (phonon_parameter_file.empty() ? "repository" : phonon_parameter_file) << '\n';
    os << "impact_ionization_parameter_set = "
       << (enable_impact_ionization ? impact_ionization_parameter_set : "disabled") << '\n';
    os << "mesh_file = " << mesh_file << '\n';
    os << "phonon_scattering_rates_file = " << rates_file << '\n';
    os << "n_particles = " << nb_particles << '\n';
    os << "n_threads = " << nb_threads << '\n';
    os << "n_conduction_bands = " << nb_conduction_bands << '\n';
    os << "n_valence_bands = " << nb_valence_bands << '\n';
    os << "simulation_time_s = " << simulation_time << '\n';
    os << "temperature_K = " << temperature << '\n';
    os << "electric_field_x_V_per_cm = " << electric_field_v_per_cm.x() << '\n';
    os << "electric_field_y_V_per_cm = " << electric_field_v_per_cm.y() << '\n';
    os << "electric_field_z_V_per_cm = " << electric_field_v_per_cm.z() << '\n';
    os << "max_energy_eV = " << max_energy_eV << '\n';
    os << "self_scattering_safety_factor = " << gamma_safety << '\n';
    os << "warmup_fraction = " << warmup_fraction << '\n';
    os << "bz_domain = " << bz_domain << '\n';
    os << "carrier = " << carrier << '\n';
    os << "random_seed = " << random_seed << '\n';
    os << "impact_ionization_enabled = " << (enable_impact_ionization ? "true" : "false") << '\n';
    if (enable_impact_ionization) {
        os << "impact_ionization_model = keldysh\n";
        os << "impact_ionization_P0_s_1 = " << impact_ionization_model.m_P0 << '\n';
        os << "impact_ionization_alpha = " << impact_ionization_model.m_alpha << '\n';
        os << "impact_ionization_threshold_eV = " << impact_ionization_model.m_E_threshold << '\n';
    }
    os << "export_history = " << (export_history ? "true" : "false") << '\n';
}

}  // namespace

int main(int argc, const char** argv) try {
    TCLAP::CmdLine cmd("Full-band bulk Monte Carlo simulation.", ' ', "1.2");

    TCLAP::ValueArg<std::string> arg_mesh_file("f",
                                               "meshbandfile",
                                               "File with BZ mesh and band energies.",
                                               true,
                                               "bz.msh",
                                               "string");
    TCLAP::ValueArg<std::string> arg_phonon_file("p",
                                                 "phononfile",
                                                 "File with phonon scattering rates. If omitted, fbmc.epm searches "
                                                 "the current directory for a matching rates CSV.",
                                                 false,
                                                 "",
                                                 "string");
    TCLAP::ValueArg<std::string> arg_material("m",
                                              "material",
                                              "Symbol of the material to use (Si, Ge, GaAs, ...)",
                                              true,
                                              "Si",
                                              "string");
    TCLAP::ValueArg<std::string> arg_phonon_parameter_set("",
                                                          "phonon-params",
                                                          "Electron-phonon parameter set used for phonon dispersion "
                                                          "and scattering final-state selection.",
                                                          false,
                                                          "remi-2026",
                                                          "string");
    TCLAP::ValueArg<std::string> arg_phonon_parameter_file(
        "",
        "phonon-params-file",
        "Load electron-phonon parameters from an explicit YAML file. Mutually exclusive with --phonon-params.",
        false,
        "",
        "path");
    TCLAP::ValueArg<std::string> arg_impact_ionization_parameter_set(
        "",
        "impact-ionization-params",
        "Impact-ionization parameter set. Impact ionization is enabled by default for electrons.",
        false,
        "keldysh",
        "string");
    TCLAP::ValueArg<std::string> arg_outputdir("d", "outdir", "Output directory for results", false, "", "string");
    TCLAP::ValueArg<std::string> arg_carrier("",
                                             "carrier",
                                             "Carrier type: electron or hole.",
                                             false,
                                             "electron",
                                             "string");

    TCLAP::ValueArg<int> arg_nb_part("N", "npart", "Number of particles to simulate", false, 1, "int");
    TCLAP::ValueArg<int> arg_nb_conduction_bands("c",
                                                 "ncbands",
                                                 "Number of conduction bands to consider",
                                                 false,
                                                 -1,
                                                 "int");
    TCLAP::ValueArg<int> arg_nb_valence_bands("v", "nvbands", "Number of valence bands to consider", false, -1, "int");
    TCLAP::ValueArg<int> arg_nb_threads("j", "nthreads", "Number of threads to use", false, 1, "int");

    TCLAP::ValueArg<double> arg_max_energy("e", "maxenergy", "Maximum energy to consider (eV)", false, 10.0, "double");
    TCLAP::ValueArg<double> arg_gamma_safety("",
                                             "gamma-safety",
                                             "Safety factor applied to the maximum total scattering rate.",
                                             false,
                                             1.2,
                                             "double");
    TCLAP::ValueArg<unsigned long long> arg_random_seed("",
                                                        "seed",
                                                        "Base random seed for reproducible per-particle streams. "
                                                        "If omitted, a random seed is generated.",
                                                        false,
                                                        0ULL,
                                                        "integer");
    TCLAP::ValueArg<std::string>        arg_bz_domain("",
                                               "bz-domain",
                                               "Stored BZ domain: full or octant.",
                                               false,
                                               "full",
                                               "string");
    TCLAP::ValueArg<double>             arg_time("t", "time", "Simulation time (s)", false, 1e-12, "double");
    TCLAP::ValueArg<double>             arg_warmup_fraction("",
                                                "warmup",
                                                "Fraction of simulation time ignored for steady-state averages.",
                                                false,
                                                0.2,
                                                "double");
    TCLAP::ValueArg<double> arg_temperature("T", "temperature", "Simulation temperature (K)", false, 300.0, "double");
    TCLAP::ValueArg<double> arg_electric_field_x("",
                                                 "Ex",
                                                 "Electric field in x direction (V/cm)",
                                                 false,
                                                 0.0,
                                                 "double");
    TCLAP::ValueArg<double> arg_electric_field_y("",
                                                 "Ey",
                                                 "Electric field in y direction (V/cm)",
                                                 false,
                                                 0.0,
                                                 "double");
    TCLAP::ValueArg<double> arg_electric_field_z("",
                                                 "Ez",
                                                 "Electric field in z direction (V/cm)",
                                                 false,
                                                 0.0,
                                                 "double");

    TCLAP::SwitchArg arg_plot_with_python(
        "P",
        "plot",
        "Call a python script after the MC runs (currently not wired in this executable).",
        cmd,
        false);
    TCLAP::SwitchArg arg_plot_with_wedge(
        "w",
        "wedge",
        "Consider only the irreducible wedge of the BZ (currently not wired in this executable).",
        cmd,
        false);
    TCLAP::SwitchArg arg_test_elph("", "test-elph", "Run electron-phonon diagnostic before the MC run.", cmd, false);
    TCLAP::SwitchArg arg_export_history(
        "E",
        "export-history",
        "Export per-particle histories. Disabled by default because sweeps can create many CSV files.",
        cmd,
        false);
    TCLAP::SwitchArg arg_disable_impact_ionization("",
                                                   "disable-impact-ionization",
                                                   "Disable impact ionization scattering.",
                                                   cmd,
                                                   false);
    TCLAP::SwitchArg arg_skip_mesh_vtk("",
                                       "skip-mesh-vtk",
                                       "Do not export the static BZ mesh VTK file.",
                                       cmd,
                                       false);

    cmd.add(arg_mesh_file);
    cmd.add(arg_phonon_file);
    cmd.add(arg_material);
    cmd.add(arg_phonon_parameter_set);
    cmd.add(arg_phonon_parameter_file);
    cmd.add(arg_impact_ionization_parameter_set);
    cmd.add(arg_outputdir);
    cmd.add(arg_carrier);
    cmd.add(arg_nb_part);
    cmd.add(arg_nb_conduction_bands);
    cmd.add(arg_nb_valence_bands);
    cmd.add(arg_nb_threads);
    cmd.add(arg_max_energy);
    cmd.add(arg_gamma_safety);
    cmd.add(arg_random_seed);
    cmd.add(arg_bz_domain);
    cmd.add(arg_time);
    cmd.add(arg_warmup_fraction);
    cmd.add(arg_temperature);
    cmd.add(arg_electric_field_x);
    cmd.add(arg_electric_field_y);
    cmd.add(arg_electric_field_z);

    cmd.parse(argc, argv);

    const std::filesystem::path file_mesh                       = arg_mesh_file.getValue();
    std::filesystem::path       file_phonon_scattering          = arg_phonon_file.getValue();
    const std::string           material_symbol                 = arg_material.getValue();
    const std::string           phonon_parameter_set            = arg_phonon_parameter_set.getValue();
    const std::string           phonon_parameter_file           = arg_phonon_parameter_file.getValue();
    const std::string           impact_ionization_parameter_set = arg_impact_ionization_parameter_set.getValue();
    const std::string           init_output_directory           = arg_outputdir.getValue();
    const std::string           carrier_name                    = arg_carrier.getValue();
    const auto                  carrier_type                    = parse_particle_type(carrier_name);
    if (arg_phonon_parameter_file.isSet() && arg_phonon_parameter_set.isSet()) {
        throw std::invalid_argument("--phonon-params-file and --phonon-params are mutually exclusive");
    }

    const int nb_threads          = arg_nb_threads.getValue();
    const int nb_valence_bands    = arg_nb_valence_bands.getValue();
    const int nb_conduction_bands = arg_nb_conduction_bands.getValue();
    const int nb_particles        = arg_nb_part.getValue();

    const double        max_energy_eV             = arg_max_energy.getValue();
    const double        gamma_safety              = arg_gamma_safety.getValue();
    const double        simulation_time_s         = arg_time.getValue();
    const double        warmup_fraction           = arg_warmup_fraction.getValue();
    const double        temperature_K             = arg_temperature.getValue();
    const double        electric_field_x_V_per_cm = arg_electric_field_x.getValue();
    const double        electric_field_y_V_per_cm = arg_electric_field_y.getValue();
    const double        electric_field_z_V_per_cm = arg_electric_field_z.getValue();
    const bool          export_history            = arg_export_history.getValue();
    const bool          enable_impact_ionization  = !arg_disable_impact_ionization.getValue();
    const bool          skip_mesh_vtk             = arg_skip_mesh_vtk.getValue();
    const std::string   bz_domain_name            = arg_bz_domain.getValue();
    const std::uint64_t random_seed               = [&]() {
        if (arg_random_seed.isSet()) {
            return static_cast<std::uint64_t>(arg_random_seed.getValue());
        }
        std::random_device  random_device;
        const std::uint64_t high = static_cast<std::uint64_t>(random_device());
        const std::uint64_t low  = static_cast<std::uint64_t>(random_device());
        return (high << 32U) ^ low;
    }();

    require_positive_int(nb_threads, "--nthreads");
    require_positive_int(nb_particles, "--npart");
    if (nb_conduction_bands < -1) {
        throw std::invalid_argument("--ncbands must be -1 or non-negative");
    }
    if (nb_valence_bands < -1) {
        throw std::invalid_argument("--nvbands must be -1 or non-negative");
    }
    if (carrier_type == uepm::fbmc::particle_type::hole && enable_impact_ionization) {
        throw std::invalid_argument(
            "Impact ionization is not available for holes; pass --disable-impact-ionization");
    }
    require_positive(max_energy_eV, "--maxenergy");
    require_finite(gamma_safety, "--gamma-safety");
    if (gamma_safety < 1.0) {
        throw std::invalid_argument("--gamma-safety must be at least one");
    }
    require_positive(simulation_time_s, "--time");
    require_finite(warmup_fraction, "--warmup");
    if (warmup_fraction < 0.0 || warmup_fraction >= 1.0) {
        throw std::invalid_argument("--warmup must be in [0, 1)");
    }
    require_positive(temperature_K, "--temperature");
    require_finite(electric_field_x_V_per_cm, "--Ex");
    require_finite(electric_field_y_V_per_cm, "--Ey");
    require_finite(electric_field_z_V_per_cm, "--Ez");
    const uepm::mesh_bz::vector3 electric_field_V_per_cm{
        electric_field_x_V_per_cm,
        electric_field_y_V_per_cm,
        electric_field_z_V_per_cm,
    };
    const uepm::mesh_bz::BZDomainMode bz_domain_mode = [&]() {
        if (bz_domain_name == "full") {
            return uepm::mesh_bz::BZDomainMode::full;
        }
        if (bz_domain_name == "octant" || bz_domain_name == "positive-octant") {
            return uepm::mesh_bz::BZDomainMode::positive_octant;
        }
        throw std::invalid_argument("--bz-domain must be 'full' or 'octant'");
    }();
    if (!arg_random_seed.isSet()) {
        fmt::print("Generated random seed: {}\n", random_seed);
    }

    require_existing_file(file_mesh, "Mesh file");
    if (file_phonon_scattering.empty()) {
        file_phonon_scattering = detect_phonon_rates_file(std::filesystem::current_path());
        fmt::print("Auto-detected phonon scattering-rate file: {}\n", file_phonon_scattering.string());
    }
    require_existing_file(file_phonon_scattering, "Phonon scattering-rate file");

    const uepm::physics::material_repository material_repository;
    uepm::fbmc::KeldyshImpactIonization      impact_ionization_model;
    if (enable_impact_ionization) {
        impact_ionization_model = uepm::fbmc::load_keldysh_impact_ionization(material_repository,
                                                                             material_symbol,
                                                                             impact_ionization_parameter_set);
        fmt::print("Loaded Keldysh impact ionization: P0={:.6e} s^-1, alpha={:.6g}, threshold={:.6g} eV\n",
                   impact_ionization_model.m_P0,
                   impact_ionization_model.m_alpha,
                   impact_ionization_model.m_E_threshold);
    }

    const auto output_dir =
        make_output_directory(init_output_directory, material_symbol, temperature_K, electric_field_V_per_cm);

    write_run_metadata(output_dir,
                       file_mesh.string(),
                       file_phonon_scattering.string(),
                       material_symbol,
                       arg_phonon_parameter_file.isSet() ? "external-file" : phonon_parameter_set,
                       phonon_parameter_file,
                       impact_ionization_parameter_set,
                       nb_particles,
                       nb_threads,
                       nb_conduction_bands,
                       nb_valence_bands,
                       simulation_time_s,
                       temperature_K,
                       electric_field_V_per_cm,
                       max_energy_eV,
                       gamma_safety,
                       warmup_fraction,
                       uepm::mesh_bz::bz_domain_mode_name(bz_domain_mode),
                       carrier_name,
                       random_seed,
                       enable_impact_ionization,
                       impact_ionization_model,
                       export_history);

    if (arg_plot_with_python.getValue()) {
        fmt::print(stderr, "[warn] --plot is parsed but not wired in this executable yet.\n");
    }
    if (arg_plot_with_wedge.getValue()) {
        fmt::print(stderr, "[warn] --wedge is parsed but not wired in this executable yet.\n");
    }

    uepm::pseudopotential::Materials materials;
    materials.load_material(material_repository, material_symbol, "chel");
    const uepm::pseudopotential::epm_material current_material = materials.materials.at(material_symbol);

    uepm::mesh_bz::ElectronPhonon mesh(current_material);
    mesh.set_domain_mode(bz_domain_mode);
    if (mesh.stores_positive_octant()) {
        fmt::print(stderr,
                   "[warn] octant mode assumes reflection symmetry under independent x/y/z sign changes "
                   "for band energies and scalar scattering rates.\n");
    }
    mesh.set_number_threads_mesh_ops(nb_threads);
    mesh.set_max_energy_global(max_energy_eV);

    mesh.read_mesh_geometry_from_msh_file(file_mesh.string());
    mesh.build_search_tree();

    const bool shift_conduction_band     = true;
    const bool set_positive_valence_band = carrier_type == uepm::fbmc::particle_type::hole;
    mesh.read_mesh_bands_from_msh_file(file_mesh.string(),
                                       nb_conduction_bands,
                                       nb_valence_bands,
                                       shift_conduction_band,
                                       set_positive_valence_band);

    const auto mesh_carrier_type = carrier_type == uepm::fbmc::particle_type::electron
                                       ? uepm::mesh_bz::MeshParticleType::conduction
                                       : uepm::mesh_bz::MeshParticleType::valence;
    mesh.set_particle_type(mesh_carrier_type);
    const auto nb_elph_bands = mesh.get_number_bands(mesh_carrier_type);
    if (nb_elph_bands == 0) {
        throw std::runtime_error(fmt::format("FBMC requires at least one {} band in the mesh", carrier_name));
    }
    mesh.set_nb_bands_elph(nb_elph_bands);
    mesh.set_temperature(temperature_K);

    const auto vtk_file = output_dir / "mesh_vtk.vtk";
    if (!skip_mesh_vtk && !std::filesystem::exists(vtk_file)) {
        mesh.export_energies_and_gradients_to_vtk(vtk_file.string());
    }

    if (arg_phonon_parameter_file.isSet()) {
        mesh.load_phonon_parameters_from_file(phonon_parameter_file);
    } else {
        mesh.load_phonon_parameters(material_repository, phonon_parameter_set);
    }
    mesh.export_phonon_dispersion((output_dir / "phonon_dispersion.data").string());

    mesh.read_phonon_scattering_rates_from_file(file_phonon_scattering.string());

    if (arg_test_elph.getValue()) {
        mesh.test_elph();
    }

    const std::string           timestamp  = std::to_string(std::time(nullptr));
    const std::filesystem::path fileprefix = output_dir / fmt::format("simulation_results_{}", timestamp);

    constexpr double v_per_cm_to_v_per_m = 1.0e2;

    uepm::fbmc::bulk_fbmc_simulation_config config;
    config.m_number_of_particles           = static_cast<std::size_t>(nb_particles);
    config.m_electric_field                = electric_field_V_per_cm * v_per_cm_to_v_per_m;
    config.m_lattice_temperature           = temperature_K;
    config.m_final_time                    = simulation_time_s;
    config.m_warmup_fraction               = warmup_fraction;
    config.m_max_energy_eV                 = max_energy_eV;
    config.m_self_scattering_safety_factor = gamma_safety;
    config.m_nb_threads                    = static_cast<std::size_t>(nb_threads);
    config.m_enable_impact_ionization      = enable_impact_ionization;
    config.m_impact_ionization_model       = impact_ionization_model;
    config.m_record_history                = export_history;
    config.m_random_seed                   = random_seed;
    config.m_particle_type                 = carrier_type;
    if (export_history) {
        config.m_history_export_prefix = fileprefix.string();
    }

    uepm::fbmc::bulk_fbmc_simulation sim(&mesh, config);

    const auto start = std::chrono::high_resolution_clock::now();
    sim.run_simulation();
    const auto end = std::chrono::high_resolution_clock::now();

    const std::chrono::duration<double> elapsed = end - start;
    fmt::print("Simulation completed in {:.3f} seconds.\n", elapsed.count());

    sim.extract_stats_and_export((fileprefix.string() + "_stats.csv"));
    sim.extract_stats_and_export((output_dir / "observables.csv").string());

    return 0;

} catch (const TCLAP::ArgException& e) {
    std::cerr << "TCLAP error: " << e.error() << " for arg " << e.argId() << '\n';
    return 2;
} catch (const std::out_of_range& e) {
    std::cerr << "Configuration error: " << e.what() << '\n';
    return 3;
} catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << '\n';
    return 1;
}
