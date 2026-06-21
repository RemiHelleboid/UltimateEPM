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

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <tclap/CmdLine.h>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "yaml-cpp/yaml.h"
#include "BandStructure.h"
#include "Options.h"
#include "electron_phonon.hpp"
#include "epm_material.hpp"
#include "fermi_level.hpp"

template <typename Derived>
struct fmt::formatter<Eigen::DenseBase<Derived>> : fmt::ostream_formatter {};

namespace {

void require_positive(int value, const std::string& option_name) {
    if (value <= 0) {
        throw std::invalid_argument(fmt::format("{} must be positive", option_name));
    }
}

void require_positive(double value, const std::string& option_name) {
    if (!(value > 0.0) || !std::isfinite(value)) {
        throw std::invalid_argument(fmt::format("{} must be finite and positive", option_name));
    }
}

void require_band_count(int value, const std::string& option_name) {
    if (value < -1) {
        throw std::invalid_argument(fmt::format("{} must be -1 or non-negative", option_name));
    }
}

uepm::mesh_bz::MeshParticleType parse_carrier_type(const std::string& value) {
    if (value == "electron") {
        return uepm::mesh_bz::MeshParticleType::conduction;
    }
    if (value == "hole") {
        return uepm::mesh_bz::MeshParticleType::valence;
    }
    throw std::invalid_argument("--carrier must be 'electron' or 'hole'");
}

std::filesystem::path make_output_directory(const std::string& requested) {
    std::filesystem::path outdir = requested.empty() ? std::filesystem::path(".") : std::filesystem::path(requested);
    std::filesystem::create_directories(outdir);
    if (!std::filesystem::is_directory(outdir)) {
        throw std::runtime_error(fmt::format("Output path is not a directory: {}", outdir.string()));
    }
    return outdir;
}

void write_kernel_metadata(const std::filesystem::path& kernel_file,
                           const std::filesystem::path& mesh_file,
                           const std::filesystem::path& phonon_parameter_file,
                           const std::string&           material,
                           const std::string&           carrier,
                           int                          nb_conduction_bands,
                           int                          nb_valence_bands,
                           double                       temperature_K,
                           double                       energy_window_eV,
                           const std::string&           bz_domain) {
    const YAML::Node profile = YAML::LoadFile(phonon_parameter_file.string());
    YAML::Node       metadata;
    metadata["schema_version"]          = 1;
    metadata["model"]                   = "electron_phonon_kernel";
    metadata["material"]                = material;
    metadata["carrier"]                 = carrier;
    metadata["mesh_file"]               = std::filesystem::absolute(mesh_file).lexically_normal().string();
    metadata["mesh_size_bytes"]         = std::filesystem::file_size(mesh_file);
    metadata["phonon_parameter_file"]   = std::filesystem::absolute(phonon_parameter_file).lexically_normal().string();
    metadata["parameter_set"]           = profile["parameter_set"];
    metadata["n_conduction_bands"]      = nb_conduction_bands;
    metadata["n_valence_bands"]         = nb_valence_bands;
    metadata["temperature_K"]           = temperature_K;
    metadata["energy_window_eV"]        = energy_window_eV;
    metadata["bz_domain"]               = bz_domain;
    metadata["Radius-WS"]               = profile["Radius-WS"];
    metadata["dispersion"]              = profile["dispersion"];

    const auto metadata_file = std::filesystem::path(kernel_file.string() + ".meta.yaml");
    std::ofstream stream(metadata_file);
    if (!stream) {
        throw std::runtime_error("Could not write kernel metadata file " + metadata_file.string());
    }
    stream << metadata;
    fmt::print("Exported kernel metadata to {}\n", metadata_file.string());
}

}  // namespace

int export_result_mobility(const std::string&     filename,
                           const Eigen::Matrix3d& mu_tensor,
                           double                 mu_iso,
                           double                 Ef,
                           const Options&         my_options,
                           double                 max_energy,
                           double                 temperature,
                           std::size_t            nb_vtx,
                           std::size_t            nb_conduction_bands,
                           std::size_t            nb_valence_bands) {
    std::ofstream    file(filename);
    constexpr double mu_to_cm2Vs = 1e4;  // m^2/(V·s) to cm^2/(V·s)

    if (file) {
        file << "# Mobility tensor computed with EPP\n";
        file << "# epm_material : " << my_options.materialName << "\n";
        file << "# Number of vertices : " << nb_vtx << "\n";
        file << "# Number of conduction bands : " << nb_conduction_bands << "\n";
        file << "# Number of valence bands : " << nb_valence_bands << "\n";
        file << "# Energy range in eV : " << max_energy << "\n";
        file << "# Temperature in Kelvin : " << temperature << "\n";
        file << "# Fermi level in eV : " << Ef << "\n";
        file << "# Mobility tensor in cm^2/(V·s)\n";
        file << mu_tensor * mu_to_cm2Vs << "\n";
        file << "# Isotropic mobility in cm^2/(V·s)\n";
        file << mu_iso * mu_to_cm2Vs << "\n";
        file.close();
        fmt::print("Mobility tensor written to {}\n", filename);
        return 0;
    } else {
        fmt::print(std::cerr, "Error: could not write to file {}\n", filename);
        return 1;
    }
    return 0;
}

int main(int argc, char const* argv[]) {
    TCLAP::CmdLine               cmd("Electron-phonon rate and mobility utility.", ' ', "1.1");
    TCLAP::ValueArg<std::string> arg_mesh_file("f",
                                               "meshbandfile",
                                               "File with BZ mesh and bands energy.",
                                               true,
                                               "bz.msh",
                                               "string");
    TCLAP::ValueArg<std::string> arg_phonon_rates("P",
                                                  "phononrates",
                                                  "File to load BZ phonon rates.",
                                                  false,
                                                  "bz_phonon.csv",
                                                  "string");
    TCLAP::ValueArg<std::string> arg_material("m",
                                              "material",
                                              "Symbol of the material to use (Si, Ge, GaAs, ...)",
                                              true,
                                              "Si",
                                              "string");
    TCLAP::ValueArg<std::string> arg_phonon_parameter_set("",
                                                          "phonon-params",
                                                          "Electron-phonon parameter set, e.g. kamakura, michaillat, "
                                                          "or fischetti.",
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
    TCLAP::ValueArg<std::string> arg_output_dir("d",
                                                "outdir",
                                                "Output directory for generated files.",
                                                false,
                                                "",
                                                "string");
    TCLAP::ValueArg<std::string> arg_carrier("",
                                             "carrier",
                                             "Carrier type for phonon rates: electron or hole.",
                                             false,
                                             "electron",
                                             "string");
    TCLAP::ValueArg<std::string> arg_rates_output("",
                                                  "rates-out",
                                                  "Output CSV for computed electron-phonon rates.",
                                                  false,
                                                  "",
                                                  "string");
    TCLAP::ValueArg<std::string> arg_kernel_file("",
                                                 "kernel-file",
                                                 "Load DP-independent phonon rate kernels from this CSV. The mesh, "
                                                 "bands, temperature, and phonon dispersion must match.",
                                                 false,
                                                 "",
                                                 "string");
    TCLAP::ValueArg<std::string> arg_kernels_output("",
                                                    "kernels-out",
                                                    "Output CSV for DP-independent phonon rate kernels.",
                                                    false,
                                                    "",
                                                    "string");
    TCLAP::ValueArg<int>         arg_nb_energies("e", "nenergy", "Number of energies to compute", false, 250, "int");
    TCLAP::ValueArg<int>         arg_nb_conduction_bands("c",
                                                 "ncbands",
                                                 "Number of conduction bands to consider",
                                                 false,
                                                 -1,
                                                 "int");
    TCLAP::ValueArg<int> arg_nb_valence_bands("v", "nvbands", "Number of valence bands to consider", false, -1, "int");
    TCLAP::ValueArg<int> arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    TCLAP::ValueArg<double> arg_temperature("T", "temperature", "Temperature in Kelvin.", false, 300.0, "double");
    TCLAP::ValueArg<double> arg_band_gap("g", "bandgap", "Band gap energy in eV.", false, 1.12, "double");
    TCLAP::ValueArg<double> arg_energy_range("E",
                                             "energy_window",
                                             "Energy window around the band gap to consider (in eV).",
                                             false,
                                             0.3,
                                             "double");
    TCLAP::SwitchArg        arg_export_rates("X", "export-rates", "Export electron-phonon rates in k.", false);
    TCLAP::SwitchArg arg_export_kernels("", "export-kernels", "Export DP-independent phonon rate kernels in k.", false);
    TCLAP::SwitchArg arg_rates_only("",
                                    "rates-only",
                                    "Stop after exporting rates/kernels; skip diagnostics, mobility, and plotting.",
                                    false);
    TCLAP::SwitchArg arg_skip_mesh_vtk("",
                                       "skip-mesh-vtk",
                                       "Do not export the static BZ mesh VTK file.",
                                       false);
    TCLAP::SwitchArg plot_with_python("p",
                                      "plot",
                                      "Call a python script after the computation to plot the band structure.",
                                      false);
    TCLAP::SwitchArg use_irr_wedge("w", "wedge", "Consider only the irreducible wedge of the BZ.", false);
    TCLAP::ValueArg<std::string> arg_bz_domain("",
                                               "bz-domain",
                                               "Stored BZ domain: full or octant.",
                                               false,
                                               "full",
                                               "string");
    TCLAP::SwitchArg             plot_with_knkpnp("K",
                                      "knkpnp",
                                      "Compute and store the full (n,k) -> (n',k') transition rate matrices.",
                                      false);
    cmd.add(plot_with_python);
    cmd.add(arg_mesh_file);
    cmd.add(arg_material);
    cmd.add(arg_phonon_parameter_set);
    cmd.add(arg_phonon_parameter_file);
    cmd.add(arg_output_dir);
    cmd.add(arg_carrier);
    cmd.add(arg_rates_output);
    cmd.add(arg_kernel_file);
    cmd.add(arg_kernels_output);
    cmd.add(arg_nb_conduction_bands);
    cmd.add(arg_nb_valence_bands);
    cmd.add(arg_nb_energies);
    cmd.add(arg_nb_threads);
    cmd.add(use_irr_wedge);
    cmd.add(arg_bz_domain);
    cmd.add(plot_with_knkpnp);
    cmd.add(arg_temperature);
    cmd.add(arg_energy_range);
    cmd.add(arg_export_rates);
    cmd.add(arg_export_kernels);
    cmd.add(arg_rates_only);
    cmd.add(arg_skip_mesh_vtk);
    cmd.add(arg_phonon_rates);
    cmd.add(arg_band_gap);
    cmd.parse(argc, argv);

    auto start = std::chrono::high_resolution_clock::now();

    uepm::pseudopotential::Materials         materials;
    const uepm::physics::material_repository material_repository;
    materials.load_material(material_repository, arg_material.getValue(), "chel");

    Options my_options;
    my_options.materialName                     = arg_material.getValue();
    my_options.nrLevels                         = arg_nb_conduction_bands.getValue() + arg_nb_valence_bands.getValue();
    my_options.nrThreads                        = arg_nb_threads.getValue();
    const int         number_energies           = arg_nb_energies.getValue();
    const int         nb_conduction_bands       = arg_nb_conduction_bands.getValue();
    const int         nb_valence_bands          = arg_nb_valence_bands.getValue();
    const double      max_energy                = arg_energy_range.getValue();  // eV
    const double      temperature               = arg_temperature.getValue();
    bool              irreducible_wedge_only    = use_irr_wedge.getValue();
    const std::string bz_domain_name            = arg_bz_domain.getValue();
    const std::string mesh_band_input_file      = arg_mesh_file.getValue();
    const std::string phonon_parameter_set      = arg_phonon_parameter_set.getValue();
    const bool        phonon_parameter_file_set = arg_phonon_parameter_file.isSet();
    const std::filesystem::path phonon_parameter_file =
        phonon_parameter_file_set
            ? std::filesystem::path(arg_phonon_parameter_file.getValue())
            : material_repository.parameter_file(arg_material.getValue(), "electron_phonon", phonon_parameter_set);
    const std::string carrier_name              = arg_carrier.getValue();
    const auto        carrier_type              = parse_carrier_type(carrier_name);
    const bool        shift_conduction_band     = true;
    const bool        set_positive_valence_band = carrier_type == uepm::mesh_bz::MeshParticleType::valence;
    const bool        export_rates              = arg_export_rates.getValue();
    const bool        export_kernels            = arg_export_kernels.getValue();
    const bool        rates_only                = arg_rates_only.getValue();
    const bool        skip_mesh_vtk             = arg_skip_mesh_vtk.getValue();
    bool              phonon_rates_provided     = arg_phonon_rates.isSet();
    const bool        kernel_file_provided      = arg_kernel_file.isSet();
    std::string       phonon_rates_file         = "";
    if (phonon_rates_provided) {
        phonon_rates_file = arg_phonon_rates.getValue();
    }
    if (phonon_rates_provided && kernel_file_provided) {
        throw std::invalid_argument("--phononrates and --kernel-file are mutually exclusive");
    }
    if (phonon_parameter_file_set && arg_phonon_parameter_set.isSet()) {
        throw std::invalid_argument("--phonon-params-file and --phonon-params are mutually exclusive");
    }
    if (rates_only && !export_rates && !export_kernels) {
        throw std::invalid_argument("--rates-only requires --export-rates and/or --export-kernels");
    }
    double     band_gap   = arg_band_gap.getValue();
    const auto output_dir = make_output_directory(arg_output_dir.getValue());

    require_positive(my_options.nrThreads, "--nthreads");
    require_positive(number_energies, "--nenergy");
    require_band_count(nb_conduction_bands, "--ncbands");
    require_band_count(nb_valence_bands, "--nvbands");
    require_positive(max_energy, "--energy_window");
    require_positive(temperature, "--temperature");
    require_positive(band_gap, "--bandgap");
    const uepm::mesh_bz::BZDomainMode bz_domain_mode = [&]() {
        if (bz_domain_name == "full") {
            return uepm::mesh_bz::BZDomainMode::full;
        }
        if (bz_domain_name == "octant" || bz_domain_name == "positive-octant") {
            return uepm::mesh_bz::BZDomainMode::positive_octant;
        }
        throw std::invalid_argument("--bz-domain must be 'full' or 'octant'");
    }();
    uepm::pseudopotential::epm_material current_material = materials.materials.at(arg_material.getValue());

    uepm::mesh_bz::ElectronPhonon ElectronPhonon{current_material};
    ElectronPhonon.set_domain_mode(bz_domain_mode);
    if (ElectronPhonon.stores_positive_octant()) {
        fmt::print(stderr,
                   "[warn] octant mode assumes reflection symmetry under independent x/y/z sign changes "
                   "for band energies and scalar scattering rates.\n");
    }
    ElectronPhonon.set_temperature(temperature);
    ElectronPhonon.set_number_threads_mesh_ops(my_options.nrThreads);
    ElectronPhonon.set_max_energy_global(max_energy);

    ElectronPhonon.read_mesh_geometry_from_msh_file(mesh_band_input_file);
    ElectronPhonon.build_search_tree();
    ElectronPhonon.read_mesh_bands_from_msh_file(mesh_band_input_file,
                                                 nb_conduction_bands,
                                                 nb_valence_bands,
                                                 shift_conduction_band,
                                                 set_positive_valence_band);
    const std::string vtk_file = (output_dir / "mesh_vtk.vtk").string();
    if (!skip_mesh_vtk && !std::filesystem::exists(vtk_file)) {
        ElectronPhonon.export_energies_and_gradients_to_vtk(vtk_file);
    }

    if (phonon_parameter_file_set) {
        ElectronPhonon.load_phonon_parameters_from_file(arg_phonon_parameter_file.getValue());
    } else {
        ElectronPhonon.load_phonon_parameters(material_repository, phonon_parameter_set);
    }
    ElectronPhonon.set_particle_type(carrier_type);
    const auto nb_elph_bands = ElectronPhonon.get_number_bands(carrier_type);
    if (nb_elph_bands == 0) {
        throw std::runtime_error(fmt::format("elph.epm requires at least one {} band in the mesh", carrier_name));
    }
    ElectronPhonon.set_nb_bands_elph(nb_elph_bands);

    std::size_t           nb_vtx = ElectronPhonon.get_number_vertices();
    std::filesystem::path name_path(mesh_band_input_file);
    std::string           name_stem = name_path.stem().string();
    auto stamp_params = fmt::format("_T{}K_C{}V{}_N{}", temperature, nb_conduction_bands, nb_valence_bands, nb_vtx);
    const std::filesystem::path prefix_export_path = output_dir / (name_stem + stamp_params);
    std::string                 prefix_export      = prefix_export_path.string();

    const double energy_windows_guard = 10.0 * uepm::constants::k_b_eV * temperature;
    if (max_energy < energy_windows_guard) {
        fmt::print(
            "Warning: energy window {:.3f} eV is small compared to thermal energy scale {:.3f} eV at T = {:.1f} K.\n",
            max_energy,
            energy_windows_guard,
            temperature);
    }
    if (phonon_rates_provided) {
        ElectronPhonon.read_phonon_scattering_rates_from_file(phonon_rates_file);
    } else if (kernel_file_provided) {
        ElectronPhonon.read_phonon_rate_kernels_from_file(arg_kernel_file.getValue());
    } else {
        const bool build_parameterized_rates = !rates_only || export_rates;
        ElectronPhonon.compute_phonon_rates_over_mesh(max_energy, irreducible_wedge_only, build_parameterized_rates);
    }

    if (export_rates && !phonon_rates_provided) {
        const std::string rates_file =
            arg_rates_output.isSet() ? arg_rates_output.getValue() : (output_dir / "phonon_rates.csv").string();
        ElectronPhonon.export_rate_values(rates_file);
    }
    if (export_kernels) {
        if (phonon_rates_provided) {
            throw std::invalid_argument("Cannot export kernels from an ordinary phonon-rate file");
        }
        const std::string kernels_file = arg_kernels_output.isSet() ? arg_kernels_output.getValue()
                                                                    : (output_dir / "phonon_rate_kernels.csv").string();
        ElectronPhonon.export_rate_kernels(kernels_file);
        write_kernel_metadata(kernels_file,
                              mesh_band_input_file,
                              phonon_parameter_file,
                              my_options.materialName,
                              carrier_name,
                              nb_conduction_bands,
                              nb_valence_bands,
                              temperature,
                              max_energy,
                              bz_domain_name);
    }
    if (rates_only) {
        fmt::print("Completed requested phonon export workflow.\n");
        return 0;
    }
    if (carrier_type == uepm::mesh_bz::MeshParticleType::valence) {
        fmt::print("Completed hole-phonon rate workflow.\n");
        return 0;
    }

    ElectronPhonon.test_elph();

    ElectronPhonon.apply_scissor(band_gap);  // eV
    // Solve for Fermi level and export CSV
    uepm::mesh_bz::fermi::Options fermi_options;
    fermi_options.nE                = 250;  // number of energy points for DOS interpolation
    fermi_options.threads           = my_options.nrThreads;
    fermi_options.use_interp        = false;        // use interpolation when computing DOS at given energy
    fermi_options.T_K               = temperature;  // temperature for Fermi-Dirac
    fermi_options.dop.Nd_cm3        = 0.0;          // intrinsic phonon-limited mobility target
    fermi_options.dop.Na_cm3        = 0.0;
    const bool use_iw               = !ElectronPhonon.stores_positive_octant();
    fermi_options.abs_max_energy_eV = 1.0;  // absolute max energy to consider (both conduction and valence)

    auto result = uepm::mesh_bz::fermi::solve_fermi(ElectronPhonon, fermi_options, use_iw);
    if (result.success) {
        fmt::print("Fermi level found: EF = {:.6f} eV\n", result.EF_eV);
        fmt::print("  p = {:.6e} cm^-3\n", result.p_m3 * 1e-6);
        fmt::print("  n = {:.6e} cm^-3\n", result.n_m3 * 1e-6);
    } else {
        fmt::print(std::cerr, "Fermi level not found.\n");
        return 1;
    }
    const double Ef = result.EF_eV;
    const double T  = temperature;

    const auto       mu_tensor   = ElectronPhonon.compute_electron_MRTA_mobility_tensor(Ef, T);
    const double     mu_iso      = ElectronPhonon.compute_electron_MRTA_mobility_isotropic(Ef, T);
    constexpr double mu_to_cm2Vs = 1e4;  // m^2/(V·s) to cm^2/(V·s)
    Eigen::Matrix3d  M           = mu_tensor * mu_to_cm2Vs;
    fmt::print("\n\nAt T = {:.1f} K and EF = {:.6f} eV:\n\n", temperature, Ef);
    fmt::print("μ_iso = {:.12e} cm^2/(V·s)\n\n", mu_iso * mu_to_cm2Vs);
    fmt::print("tensor = \n{} cm^2/(V*s)\n\n\n", fmt::streamed(M));

    double mean_energy = ElectronPhonon.mean_electron_energy_equilibrium(Ef, T, true);
    fmt::print("Mean electron energy above CBM at equilibrium: {:.6f} eV\n", mean_energy);

    std::string output_mobility     = name_stem + "_mobility.txt";
    auto        out                 = name_stem + stamp_params + "_mobility.txt";
    std::string out_rates_vs_energy = prefix_export + "_elph_results.csv";

    export_result_mobility(out_rates_vs_energy,
                           mu_tensor,
                           mu_iso,
                           Ef,
                           my_options,
                           max_energy,
                           temperature,
                           nb_vtx,
                           nb_conduction_bands,
                           nb_valence_bands);

    auto stop     = std::chrono::high_resolution_clock::now();
    auto duration = stop - start;
    fmt::print("\nTotal time : {:.2f} seconds\n\n\n", std::chrono::duration<double>(duration).count());

    if (plot_with_python.getValue()) {
        std::string  rates_vs_Energy_file = prefix_export + "_rates_vs_energy.csv";
        const double energy_step          = 0.001;  // energy step in eV
        ElectronPhonon.compute_plot_electron_phonon_rates_vs_energy_over_mesh(max_energy,
                                                                              energy_step,
                                                                              rates_vs_Energy_file);
        std::string command =
            "python3 " + std::string(PROJECT_SRC_DIR) + "/python/plots/plot_phonon_rate.py -f " + prefix_export;
        fmt::print("Running command: {}\n", command);
        int pyRes = std::system(command.c_str());
        if (pyRes != 0) {
            fmt::print(std::cerr, "Error: Python script returned non-zero exit code {}\n", pyRes);
        }
    }

    return 0;
}
