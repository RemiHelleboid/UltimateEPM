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
#include <filesystem>
#include <fstream>
#include <iostream>

#include "BandStructure.h"
#include "Material.h"
#include "Options.h"
#include "bulk_amc_simulation.hpp"
#include "device_amc_simulation.hpp"
#include "mesh.hpp"
#include "msh_file.hpp"
#include "particle_amc.hpp"

int main(int argc, const char** argv) {
    std::cout << "Hello Analytic MC!" << std::endl;

    TCLAP::CmdLine               cmd("Anaytical MC PROGRAM. SINGLE PARTICLE MONTE CARLO SIMULATION.", ' ', "1.0");
    TCLAP::ValueArg<std::string> arg_device_mesh("", "device-mesh", "Path to the device mesh file", true, "", "string");
    TCLAP::ValueArg<std::string> arg_material("m", "material", "Symbol of the material to use (Si, Ge, GaAs, ...)", true, "Si", "string");
    TCLAP::ValueArg<std::string> arg_outputdir("d", "outdir", "Output directory for results", false, "", "string");
    TCLAP::ValueArg<int>         arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    TCLAP::ValueArg<double>      arg_time("t", "time", "Simulation time (s)", false, 1e-12, "double");
    TCLAP::ValueArg<double>      arg_temperature("T", "temperature", "Simulation temperature (K)", false, 300.0, "double");
    TCLAP::ValueArg<double> arg_max_energy("e", "max-energy", "Maximum energy for scattering rate computation (eV)", false, 10.0, "double");
    TCLAP::SwitchArg        plot_with_python("P", "plot", "Call a python script after the MC Runs.", false);
    TCLAP::SwitchArg        arg_export_hist("E", "export", "Export history of all particles to csv files for post-processing.", false);
    cmd.add(plot_with_python);
    cmd.add(arg_material);
    cmd.add(arg_nb_threads);
    cmd.add(arg_temperature);
    cmd.add(arg_export_hist);
    cmd.add(arg_time);
    cmd.add(arg_max_energy);
    cmd.add(arg_outputdir);
    cmd.add(arg_device_mesh);

    cmd.parse(argc, argv);

    const std::string device_mesh_file      = arg_device_mesh.getValue();
    const std::string material_symbol       = arg_material.getValue();
    const std::string init_output_directory = arg_outputdir.getValue();
    int               nb_threads            = arg_nb_threads.getValue();
    int               nb_valence_bands      = 0;
    double            temperature           = arg_temperature.getValue();
    double            max_energy_eV         = arg_max_energy.getValue();

    uepm::file::msh_file fileMSH(device_mesh_file);
    fileMSH.read_mesh();
    fileMSH.read_states();
    uepm::mesh::mesh* p_mesh     = fileMSH.get_p_mesh();
    std::size_t       nbVertices = p_mesh->get_nb_vertices();

    uepm::device::device my_device(p_mesh);

    fmt::print("Device mesh loaded from file: {}\n", device_mesh_file);
    fmt::print("Number of vertices in the mesh: {}\n", nbVertices);

    uepm::pseudopotential::Materials materials;
    const std::string                file_material_parameters = std::string(PROJECT_SRC_DIR) + "/parameter_files/materials-chel.yaml";
    materials.load_material_parameters(file_material_parameters);
    uepm::pseudopotential::Material current_material = materials.materials.at(material_symbol);

    std::string core_filename = std::filesystem::path(device_mesh_file).stem().string();
    std::string output_dir;
    if (init_output_directory.empty()) {
        output_dir = fmt::format("amc_results_{}", core_filename);
    } else {
        output_dir = init_output_directory;
    }
    if (!std::filesystem::exists(output_dir)) {
        std::filesystem::create_directories(output_dir);
    }

    constexpr double m_to_cm = 1e2;
    // uepm::amc::bulk_amc_simulation_config bulk_env;

    // if (particle_type == "electron") {
    //     bulk_env.m_carrier_type = uepm::amc::particle_type::electron;
    // } else if (particle_type == "hole") {
    //     bulk_env.m_carrier_type = uepm::amc::particle_type::hole;
    // } else {
    //     fmt::print(stderr, "Invalid particle type: {}. Should be 'electron' or 'hole'.\n", particle_type);
    //     return 1;
    // }

    // bulk_env.m_record_history       = arg_export_hist.getValue();
    // bulk_env.m_lattice_temperature  = temperature;
    // bulk_env.m_final_time           = arg_time.getValue();
    // bulk_env.m_number_of_particles  = static_cast<std::size_t>(nb_particles);
    // bulk_env.m_electric_field       = {electric_field_x * m_to_cm, 0.0, 0.0};
    // double m_doping_concentration   = 1.0e10;                  // m^-3
    // bulk_env.m_doping_concentration = m_doping_concentration;  // m^-3

    // bulk_env.m_max_energy_eV = max_energy_eV;

    // CREATION OF THE VALLEY MODEL (SI, 6 valleys, non-parabolicity, hardcoded parameters for now, will be read from yaml later)

    fmt::print("Running bulk AMC simulation with the following parameters:\n");
    fmt::print("Material: {}\n", material_symbol);

    uepm::amc::options_device_amc options;
    options.m_t_max                = arg_time.getValue();
    options.m_time_step           = 1e-15;
    // options.m_lattice_temperature = temperature;
    options.m_nb_threads          = nb_threads;
    options.m_max_number_particle = 1000000;
    options.m_avalanche_threshold = 1000;
    options.m_activate_impact_ionization = true;
    options.m_particle_creation_activated = true;

    uepm::amc::device_amc_simulation sim(my_device, options);
    // sim.initialize();
    auto start = std::chrono::high_resolution_clock::now();
    // sim.run();
    // sim.run_self_scattering_emc();

    auto                          end     = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    fmt::print("Simulation completed in {:.3f} seconds.\n", elapsed.count());

    std::string timestamp  = std::to_string(std::time(nullptr));
    std::string fileprefix = fmt::format("{}/simulation_results_{}", output_dir, timestamp);

    // if (arg_export_hist.getValue()) {
    //     sim.export_particles_history_to_csv(fileprefix);
    // }

    // sim.export_observables_to_csv("observables.csv");
    // sim.extract_stats_and_export(fileprefix + "_stats.csv");

    return 0;
}
