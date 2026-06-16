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
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
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

    return header.find("vertex_index") != std::string::npos &&
           header.find("local_band_index") != std::string::npos &&
           header.find("energy_eV") != std::string::npos &&
           header.find("rate_ac_L_ab") != std::string::npos &&
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
        throw std::runtime_error(fmt::format(
            "No phonon-rate CSV was provided and no matching file was found in {}. "
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

std::filesystem::path make_output_directory(const std::string& requested,
                                            const std::string& material,
                                            double             temperature,
                                            double             electric_field_v_per_cm) {
    std::filesystem::path outdir;
    if (requested.empty()) {
        outdir = fmt::format("fbmc_{}_{:.1f}K_{:.3e}Vcm", material, temperature, electric_field_v_per_cm);
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

void write_run_metadata(const std::filesystem::path& outdir,
                        const std::string&           mesh_file,
                        const std::string&           rates_file,
                        const std::string&           material,
                        const std::string&           phonon_parameter_set,
                        int                          nb_particles,
                        int                          nb_threads,
                        int                          nb_conduction_bands,
                        int                          nb_valence_bands,
                        double                       simulation_time,
                        double                       temperature,
                        double                       electric_field_v_per_cm,
                        double                       max_energy_eV,
                        double                       warmup_fraction,
                        bool                         export_history) {
    const auto    meta_file = outdir / "run_info.txt";
    std::ofstream os(meta_file);
    if (!os) {
        throw std::runtime_error(fmt::format("Could not open metadata file for writing: {}", meta_file.string()));
    }

    os << "material = " << material << '\n';
    os << "phonon_parameter_set = " << phonon_parameter_set << '\n';
    os << "mesh_file = " << mesh_file << '\n';
    os << "phonon_scattering_rates_file = " << rates_file << '\n';
    os << "n_particles = " << nb_particles << '\n';
    os << "n_threads = " << nb_threads << '\n';
    os << "n_conduction_bands = " << nb_conduction_bands << '\n';
    os << "n_valence_bands = " << nb_valence_bands << '\n';
    os << "simulation_time_s = " << simulation_time << '\n';
    os << "temperature_K = " << temperature << '\n';
    os << "electric_field_V_per_cm = " << electric_field_v_per_cm << '\n';
    os << "max_energy_eV = " << max_energy_eV << '\n';
    os << "warmup_fraction = " << warmup_fraction << '\n';
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
                                                          "kamakura",
                                                          "string");
    TCLAP::ValueArg<std::string> arg_outputdir("d", "outdir", "Output directory for results", false, "", "string");

    TCLAP::ValueArg<int> arg_nb_part("N", "npart", "Number of particles to simulate", false, 1, "int");
    TCLAP::ValueArg<int> arg_nb_conduction_bands("c",
                                                 "ncbands",
                                                 "Number of conduction bands to consider",
                                                 false,
                                                 -1,
                                                 "int");
    TCLAP::ValueArg<int> arg_nb_valence_bands("v", "nvbands", "Number of valence bands to consider", false, -1, "int");
    TCLAP::ValueArg<int> arg_nb_threads("j", "nthreads", "Number of threads to use", false, 1, "int");

    TCLAP::ValueArg<double> arg_max_energy("e",
                                           "maxenergy",
                                           "Maximum energy to consider (eV)",
                                           false,
                                           10.0,
                                           "double");
    TCLAP::ValueArg<double> arg_time("t", "time", "Simulation time (s)", false, 1e-12, "double");
    TCLAP::ValueArg<double> arg_warmup_fraction("",
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

    cmd.add(arg_mesh_file);
    cmd.add(arg_phonon_file);
    cmd.add(arg_material);
    cmd.add(arg_phonon_parameter_set);
    cmd.add(arg_outputdir);
    cmd.add(arg_nb_part);
    cmd.add(arg_nb_conduction_bands);
    cmd.add(arg_nb_valence_bands);
    cmd.add(arg_nb_threads);
    cmd.add(arg_max_energy);
    cmd.add(arg_time);
    cmd.add(arg_warmup_fraction);
    cmd.add(arg_temperature);
    cmd.add(arg_electric_field_x);

    cmd.parse(argc, argv);

    const std::filesystem::path file_mesh              = arg_mesh_file.getValue();
    std::filesystem::path       file_phonon_scattering = arg_phonon_file.getValue();
    const std::string           material_symbol        = arg_material.getValue();
    const std::string           phonon_parameter_set   = arg_phonon_parameter_set.getValue();
    const std::string           init_output_directory  = arg_outputdir.getValue();

    const int nb_threads          = arg_nb_threads.getValue();
    const int nb_valence_bands    = arg_nb_valence_bands.getValue();
    const int nb_conduction_bands = arg_nb_conduction_bands.getValue();
    const int nb_particles        = arg_nb_part.getValue();

    const double max_energy_eV             = arg_max_energy.getValue();
    const double simulation_time_s         = arg_time.getValue();
    const double warmup_fraction           = arg_warmup_fraction.getValue();
    const double temperature_K             = arg_temperature.getValue();
    const double electric_field_x_V_per_cm = arg_electric_field_x.getValue();
    const bool   export_history            = arg_export_history.getValue();

    require_positive_int(nb_threads, "--nthreads");
    require_positive_int(nb_particles, "--npart");
    if (nb_conduction_bands < -1) {
        throw std::invalid_argument("--ncbands must be -1 or non-negative");
    }
    if (nb_valence_bands < -1) {
        throw std::invalid_argument("--nvbands must be -1 or non-negative");
    }
    require_positive(max_energy_eV, "--maxenergy");
    require_positive(simulation_time_s, "--time");
    require_finite(warmup_fraction, "--warmup");
    if (warmup_fraction < 0.0 || warmup_fraction >= 1.0) {
        throw std::invalid_argument("--warmup must be in [0, 1)");
    }
    require_positive(temperature_K, "--temperature");
    require_finite(electric_field_x_V_per_cm, "--Ex");

    require_existing_file(file_mesh, "Mesh file");
    if (file_phonon_scattering.empty()) {
        file_phonon_scattering = detect_phonon_rates_file(std::filesystem::current_path());
        fmt::print("Auto-detected phonon scattering-rate file: {}\n", file_phonon_scattering.string());
    }
    require_existing_file(file_phonon_scattering, "Phonon scattering-rate file");

    const auto output_dir =
        make_output_directory(init_output_directory, material_symbol, temperature_K, electric_field_x_V_per_cm);

    write_run_metadata(output_dir,
                       file_mesh.string(),
                       file_phonon_scattering.string(),
                       material_symbol,
                       phonon_parameter_set,
                       nb_particles,
                       nb_threads,
                       nb_conduction_bands,
                       nb_valence_bands,
                       simulation_time_s,
                       temperature_K,
                       electric_field_x_V_per_cm,
                       max_energy_eV,
                       warmup_fraction,
                       export_history);

    if (arg_plot_with_python.getValue()) {
        fmt::print(stderr, "[warn] --plot is parsed but not wired in this executable yet.\n");
    }
    if (arg_plot_with_wedge.getValue()) {
        fmt::print(stderr, "[warn] --wedge is parsed but not wired in this executable yet.\n");
    }

    uepm::pseudopotential::Materials         materials;
    const uepm::physics::material_repository material_repository;
    materials.load_material(material_repository, material_symbol, "chel");
    const uepm::pseudopotential::epm_material current_material = materials.materials.at(material_symbol);

    uepm::mesh_bz::ElectronPhonon mesh(current_material);
    mesh.set_number_threads_mesh_ops(nb_threads);
    mesh.set_max_energy_global(max_energy_eV);

    mesh.read_mesh_geometry_from_msh_file(file_mesh.string());
    mesh.build_search_tree();

    const bool shift_conduction_band = true;
    mesh.read_mesh_bands_from_msh_file(file_mesh.string(),
                                       nb_conduction_bands,
                                       nb_valence_bands,
                                       shift_conduction_band);

    mesh.set_particle_type(uepm::mesh_bz::MeshParticleType::conduction);
    const auto nb_elph_bands = mesh.get_number_conduction_bands();
    if (nb_elph_bands == 0) {
        throw std::runtime_error("FBMC requires at least one conduction band in the mesh");
    }
    mesh.set_nb_bands_elph(nb_elph_bands);
    mesh.set_temperature(temperature_K);

    const auto vtk_file = output_dir / "mesh_vtk.vtk";
    if (!std::filesystem::exists(vtk_file)) {
        mesh.export_energies_and_gradients_to_vtk(vtk_file.string());
    }

    mesh.load_phonon_parameters(material_repository, phonon_parameter_set);
    mesh.export_phonon_dispersion((output_dir / "phonon_dispersion.data").string());

    mesh.read_phonon_scattering_rates_from_file(file_phonon_scattering.string());

    if (arg_test_elph.getValue()) {
        mesh.test_elph();
    }

    uepm::fbmc::Bulk_environment bulk_env;
    bulk_env.m_temperature = temperature_K;

    constexpr double v_per_cm_to_v_per_m = 1.0e2;
    bulk_env.m_electric_field            = {electric_field_x_V_per_cm * v_per_cm_to_v_per_m, 0.0, 0.0};

    bulk_env.m_doping_concentration = 1.0e10;

    uepm::fbmc::Simulation_parameters sim_params;
    sim_params.m_simulation_time   = simulation_time_s;
    sim_params.m_warmup_fraction   = warmup_fraction;
    sim_params.m_export_frequency  = 10;
    sim_params.m_nb_openmp_threads = nb_threads;

    uepm::fbmc::Single_particle_simulation sim(&mesh, bulk_env, sim_params, nb_particles);

    const auto start = std::chrono::high_resolution_clock::now();
    sim.run_simulation();
    const auto end = std::chrono::high_resolution_clock::now();

    const std::chrono::duration<double> elapsed = end - start;
    fmt::print("Simulation completed in {:.3f} seconds.\n", elapsed.count());

    const std::string           timestamp  = std::to_string(std::time(nullptr));
    const std::filesystem::path fileprefix = output_dir / fmt::format("simulation_results_{}", timestamp);

    if (export_history) {
        sim.export_history(fileprefix.string());
    }
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
