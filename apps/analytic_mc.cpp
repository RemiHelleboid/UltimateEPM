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
#include "particle_amc.hpp"

int main(int argc, const char** argv) {
    std::cout << "Hello SinglePartFBMC!" << std::endl;

    TCLAP::CmdLine               cmd("FBMC PROGRAM. SINGLE PARTICLE MONTE CARLO SIMULATION.", ' ', "1.0");
    TCLAP::ValueArg<std::string> arg_material("m", "material", "Symbol of the material to use (Si, Ge, GaAs, ...)", true, "Si", "string");
    TCLAP::ValueArg<std::string> arg_outputdir("d", "outdir", "Output directory for results", false, "", "string");
    TCLAP::ValueArg<int>         arg_nb_part("N", "npart", "Number of particles to simulate", false, 1, "int");
    TCLAP::ValueArg<int>         arg_nb_conduction_bands("c", "ncbands", "Number of conduction bands to consider", false, -1, "int");
    TCLAP::ValueArg<double>      arg_max_energy("e", "maxenergy", "Maximum energy to consider (eV)", false, 1.0e10, "double");
    TCLAP::ValueArg<int>         arg_nb_valence_bands("v", "nvbands", "Number of valence bands to consider", false, -1, "int");
    TCLAP::ValueArg<int>         arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    TCLAP::ValueArg<double>      arg_time("t", "time", "Simulation time (s)", false, 1e-12, "double");
    TCLAP::ValueArg<double>      arg_temperature("T", "temperature", "Simulation temperature (K)", false, 300.0, "double");
    TCLAP::ValueArg<double>      arg_electric_field_x("", "Ex", "Electric field in x direction (V/cm)", false, 0.0, "double");
    TCLAP::SwitchArg             plot_with_python("P", "plot", "Call a python script after the MC Runs.", false);
    TCLAP::SwitchArg             arg_export_hist("E", "export", "Export history of all particles to csv files for post-processing.", false);
    cmd.add(plot_with_python);
    cmd.add(arg_material);
    cmd.add(arg_nb_conduction_bands);
    cmd.add(arg_nb_valence_bands);
    cmd.add(arg_max_energy);
    cmd.add(arg_nb_threads);
    cmd.add(arg_temperature);
    cmd.add(arg_export_hist);
    cmd.add(arg_time);
    cmd.add(arg_nb_part);
    cmd.add(arg_electric_field_x);
    cmd.add(arg_outputdir);

    cmd.parse(argc, argv);

    const std::string material_symbol       = arg_material.getValue();
    const std::string init_output_directory = arg_outputdir.getValue();
    int               nb_threads            = arg_nb_threads.getValue();
    int               nb_valence_bands      = 0;
    int               nb_conduction_bands   = arg_nb_conduction_bands.getValue();
    int               nb_particles          = arg_nb_part.getValue();
    const double      max_energy_eV         = arg_max_energy.getValue();
    double            electric_field_x      = arg_electric_field_x.getValue();
    double            temperature           = arg_temperature.getValue();

    uepm::pseudopotential::Materials materials;
    const std::string                file_material_parameters = std::string(PROJECT_SRC_DIR) + "/parameter_files/materials-chel.yaml";
    materials.load_material_parameters(file_material_parameters);
    uepm::pseudopotential::Material current_material = materials.materials.at(material_symbol);

    std::string output_dir;
    if (init_output_directory.empty()) {
        output_dir = fmt::format("fbmc_{}_{:.1f}_{}", material_symbol, temperature, electric_field_x);
    } else {
        output_dir = init_output_directory;
    }
    if (!std::filesystem::exists(output_dir)) {
        std::filesystem::create_directories(output_dir);
    }

    constexpr double                      m_to_cm = 1e2;
    uepm::amc::bulk_amc_simulation_config bulk_env;
    bulk_env.m_electric_field *= m_to_cm;

    bulk_env.m_lattice_temperature  = temperature;
    bulk_env.m_final_time           = arg_time.getValue();
    bulk_env.m_number_of_particles  = static_cast<std::size_t>(nb_particles);
    bulk_env.m_electric_field       = {electric_field_x * m_to_cm, 0.0, 0.0};
    double m_doping_concentration   = 1.0e10;                  // m^-3
    bulk_env.m_doping_concentration = m_doping_concentration;  // m^-3

    // CREATION OF THE VALLEY MODEL (SI, 6 valleys, non-parabolicity, hardcoded parameters for now, will be read from yaml later)

    uepm::amc::bulk_amc_simulation sim{bulk_env};
    sim.initialize();
    auto start = std::chrono::high_resolution_clock::now();
    sim.run();

    auto                          end     = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    fmt::print("Simulation completed in {:.3f} seconds.\n", elapsed.count());

    std::string timestamp  = std::to_string(std::time(nullptr));
    std::string fileprefix = fmt::format("{}/simulation_results_{}", output_dir, timestamp);

    if (arg_export_hist.getValue()) {
        sim.export_particles_history_to_csv(fileprefix);
    }
    // sim.extract_stats_and_export(fileprefix + "_stats.csv");

    return 0;
}
