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

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "BandStructure.h"
#include "Material.h"
#include "Options.h"
#include "bz_mesh.hpp"
#include "bz_states.hpp"
#include "electron_phonon.hpp"
#include "single_part_fbmc.hpp"

namespace {

void require_existing_file(const std::filesystem::path& path, const std::string& what) {
    if (!std::filesystem::exists(path)) {
        throw std::runtime_error(fmt::format("{} does not exist: {}", what, path.string()));
    }
    if (!std::filesystem::is_regular_file(path)) {
        throw std::runtime_error(fmt::format("{} is not a regular file: {}", what, path.string()));
    }
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
                        int                          nb_particles,
                        int                          nb_threads,
                        int                          nb_conduction_bands,
                        int                          nb_valence_bands,
                        double                       simulation_time,
                        double                       temperature,
                        double                       electric_field_v_per_cm,
                        double                       max_energy_eV) {
    const auto    meta_file = outdir / "run_info.txt";
    std::ofstream os(meta_file);
    if (!os) {
        throw std::runtime_error(fmt::format("Could not open metadata file for writing: {}", meta_file.string()));
    }

    os << "material = " << material << '\n';
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
}

}  // namespace

int main(int argc, const char** argv) try {
    std::cout << "Hello SinglePartFBMC!" << std::endl;

    TCLAP::CmdLine cmd("FBMC PROGRAM. SINGLE PARTICLE MONTE CARLO SIMULATION.", ' ', "1.1");

    TCLAP::ValueArg<std::string> arg_mesh_file("f", "meshbandfile", "File with BZ mesh and band energies.", true, "bz.msh", "string");
    TCLAP::ValueArg<std::string> arg_phonon_file("p", "phononfile", "File with phonon scattering rates.", true, "rates_all.csv", "string");
    TCLAP::ValueArg<std::string> arg_material("m", "material", "Symbol of the material to use (Si, Ge, GaAs, ...)", true, "Si", "string");
    TCLAP::ValueArg<std::string> arg_outputdir("d", "outdir", "Output directory for results", false, "", "string");

    TCLAP::ValueArg<int> arg_nb_part("N", "npart", "Number of particles to simulate", false, 1, "int");
    TCLAP::ValueArg<int> arg_nb_conduction_bands("c", "ncbands", "Number of conduction bands to consider", false, -1, "int");
    TCLAP::ValueArg<int> arg_nb_valence_bands("v", "nvbands", "Number of valence bands to consider", false, -1, "int");
    TCLAP::ValueArg<int> arg_nb_threads("j", "nthreads", "Number of threads to use", false, 1, "int");

    TCLAP::ValueArg<double> arg_max_energy("e", "maxenergy", "Maximum energy to consider (eV)", false, 1.0e10, "double");
    TCLAP::ValueArg<double> arg_time("t", "time", "Simulation time (s)", false, 1e-12, "double");
    TCLAP::ValueArg<double> arg_temperature("T", "temperature", "Simulation temperature (K)", false, 300.0, "double");
    TCLAP::ValueArg<double> arg_electric_field_x("", "Ex", "Electric field in x direction (V/cm)", false, 0.0, "double");

    TCLAP::SwitchArg arg_plot_with_python("P",
                                          "plot",
                                          "Call a python script after the MC runs (currently not wired in this executable).",
                                          cmd,
                                          false);
    TCLAP::SwitchArg arg_plot_with_wedge("w",
                                         "wedge",
                                         "Consider only the irreducible wedge of the BZ (currently not wired in this executable).",
                                         cmd,
                                         false);
    TCLAP::SwitchArg arg_test_elph("", "test-elph", "Run electron-phonon diagnostic before the MC run.", cmd, false);

    cmd.add(arg_mesh_file);
    cmd.add(arg_phonon_file);
    cmd.add(arg_material);
    cmd.add(arg_outputdir);
    cmd.add(arg_nb_part);
    cmd.add(arg_nb_conduction_bands);
    cmd.add(arg_nb_valence_bands);
    cmd.add(arg_nb_threads);
    cmd.add(arg_max_energy);
    cmd.add(arg_time);
    cmd.add(arg_temperature);
    cmd.add(arg_electric_field_x);

    cmd.parse(argc, argv);

    const std::filesystem::path file_mesh              = arg_mesh_file.getValue();
    const std::filesystem::path file_phonon_scattering = arg_phonon_file.getValue();
    const std::string           material_symbol        = arg_material.getValue();
    const std::string           init_output_directory  = arg_outputdir.getValue();

    const int nb_threads          = std::max(1, arg_nb_threads.getValue());
    const int nb_valence_bands    = arg_nb_valence_bands.getValue();
    const int nb_conduction_bands = arg_nb_conduction_bands.getValue();
    const int nb_particles        = std::max(1, arg_nb_part.getValue());

    const double max_energy_eV             = arg_max_energy.getValue();
    const double simulation_time_s         = arg_time.getValue();
    const double temperature_K             = arg_temperature.getValue();
    const double electric_field_x_V_per_cm = arg_electric_field_x.getValue();

    require_existing_file(file_mesh, "Mesh file");
    require_existing_file(file_phonon_scattering, "Phonon scattering-rate file");

    const auto output_dir = make_output_directory(init_output_directory, material_symbol, temperature_K, electric_field_x_V_per_cm);

    write_run_metadata(output_dir,
                       file_mesh.string(),
                       file_phonon_scattering.string(),
                       material_symbol,
                       nb_particles,
                       nb_threads,
                       nb_conduction_bands,
                       nb_valence_bands,
                       simulation_time_s,
                       temperature_K,
                       electric_field_x_V_per_cm,
                       max_energy_eV);

    if (arg_plot_with_python.getValue()) {
        fmt::print(stderr, "[warn] --plot is parsed but not wired in this executable yet.\n");
    }
    if (arg_plot_with_wedge.getValue()) {
        fmt::print(stderr, "[warn] --wedge is parsed but not wired in this executable yet.\n");
    }

    uepm::pseudopotential::Materials materials;
    const std::filesystem::path      file_material_parameters =
        std::filesystem::path(PROJECT_SRC_DIR) / "parameter_files" / "materials-chel.yaml";
    require_existing_file(file_material_parameters, "Material parameter file");

    materials.load_material_parameters(file_material_parameters.string());
    const uepm::pseudopotential::Material current_material = materials.materials.at(material_symbol);

    uepm::mesh_bz::ElectronPhonon mesh(current_material);
    mesh.set_number_threads_mesh_ops(nb_threads);
    mesh.set_max_energy_global(max_energy_eV);

    mesh.read_mesh_geometry_from_msh_file(file_mesh.string());
    mesh.build_search_tree();

    const bool shift_conduction_band = true;
    mesh.read_mesh_bands_from_msh_file(file_mesh.string(), nb_conduction_bands, nb_valence_bands, shift_conduction_band);

    mesh.set_particle_type(uepm::mesh_bz::MeshParticleType::conduction);
    mesh.set_nb_bands_elph(nb_conduction_bands);
    mesh.set_temperature(temperature_K);

    const auto vtk_file = output_dir / "mesh_vtk.vtk";
    if (!std::filesystem::exists(vtk_file)) {
        mesh.export_energies_and_gradients_to_vtk(vtk_file.string());
    }

    const std::filesystem::path phonon_parameter_file = std::filesystem::path(PROJECT_SRC_DIR) / "parameter_files" / "phonon_kamakura.yaml";
    require_existing_file(phonon_parameter_file, "Phonon parameter file");

    mesh.load_phonon_parameters(phonon_parameter_file.string());
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

    sim.export_history(fileprefix.string());
    sim.extract_stats_and_export((fileprefix.string() + "_stats.csv"));

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