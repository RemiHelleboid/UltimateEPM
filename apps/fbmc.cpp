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
#include "bz_mesh.hpp"
#include "bz_states.hpp"
#include "electron_phonon.hpp"
#include "single_part_fbmc.hpp"

int main(int argc, const char** argv) {
    std::cout << "Hello SinglePartFBMC!" << std::endl;

    TCLAP::CmdLine               cmd("FBMC PROGRAM. SINGLE PARTICLE MONTE CARLO SIMULATION.", ' ', "1.0");
    TCLAP::ValueArg<std::string> arg_mesh_file("f", "meshbandfile", "File with BZ mesh and bands energy.", true, "bz.msh", "string");
    TCLAP::ValueArg<std::string> arg_phonon_file("p", "phononfile", "File with phonon scattering rates.", true, "rates_all.csv", "string");
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
    TCLAP::SwitchArg             plot_with_wedge("w", "wedge", "Consider only the irreducible wedge of the BZ.", false);
    cmd.add(plot_with_python);
    cmd.add(arg_mesh_file);
    cmd.add(arg_material);
    cmd.add(arg_nb_conduction_bands);
    cmd.add(arg_nb_valence_bands);
    cmd.add(arg_max_energy);
    cmd.add(arg_nb_threads);
    cmd.add(plot_with_wedge);
    cmd.add(arg_phonon_file);
    cmd.add(arg_time);
    cmd.add(arg_nb_part);
    cmd.add(arg_electric_field_x);
    cmd.add(arg_outputdir);

    cmd.parse(argc, argv);

    const std::string file_mesh              = arg_mesh_file.getValue();
    const std::string material_symbol        = arg_material.getValue();
    const std::string file_phonon_scattering = arg_phonon_file.getValue();
    const std::string init_output_directory  = arg_outputdir.getValue();
    int               nb_threads             = arg_nb_threads.getValue();
    int               nb_valence_bands       = 0;
    int               nb_conduction_bands    = 2;
    int               nb_particles           = arg_nb_part.getValue();
    const double      max_energy_eV          = arg_max_energy.getValue();
    double            electric_field_x       = arg_electric_field_x.getValue();
    double            temperature            = arg_temperature.getValue();

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

    uepm::mesh_bz::ElectronPhonon mesh(current_material);
    mesh.set_number_threads_mesh_ops(nb_threads);
    mesh.read_mesh_geometry_from_msh_file(file_mesh);
    mesh.build_search_tree();
    // mesh.export_octree_to_vtu("octree.vtu");
    mesh.set_nb_bands_elph(nb_conduction_bands);
    mesh.set_max_energy_global(max_energy_eV);

    const bool shift_conduction_band = true;
    mesh.read_mesh_bands_from_msh_file(file_mesh, nb_conduction_bands, nb_valence_bands, shift_conduction_band);

    const std::string vtk_file = "mesh_vtk.vtk";
    if (!std::filesystem::exists(vtk_file)) {
        mesh.export_energies_and_gradients_to_vtk(vtk_file);
    }

    const std::string phonon_file = std::string(PROJECT_SRC_DIR) + "/parameter_files/phonon_kamakura.yaml";
    mesh.load_phonon_parameters(phonon_file);
    mesh.export_phonon_dispersion("phonon_dispersion.data");

    mesh.read_phonon_scattering_rates_from_file(file_phonon_scattering);

    uepm::fbmc::Bulk_environment bulk_env;
    bulk_env.m_temperature    = 300.0;
    bulk_env.m_electric_field = {electric_field_x, 0.0, 0.0};
    constexpr double m_to_cm  = 1e2;
    bulk_env.m_electric_field *= m_to_cm;
    bulk_env.m_doping_concentration = 1e10;
    mesh.set_temperature(bulk_env.m_temperature);
    mesh.set_particle_type(uepm::mesh_bz::MeshParticleType::conduction);
    mesh.set_nb_bands_elph(nb_conduction_bands);

    mesh.test_elph();

    uepm::fbmc::Simulation_parameters sim_params;
    sim_params.m_simulation_time  = arg_time.getValue();
    sim_params.m_export_frequency = 10;

    uepm::fbmc::Single_particle_simulation sim(&mesh, bulk_env, sim_params, nb_particles);

    auto start = std::chrono::high_resolution_clock::now();
    sim.run_simulation();
    auto                          end     = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    fmt::print("Simulation completed in {:.3f} seconds.\n", elapsed.count());

    std::string timestamp  = std::to_string(std::time(nullptr));
    std::string fileprefix = fmt::format("{}/simulation_results_{}", output_dir, timestamp);
    sim.export_history(fileprefix);

    return 0;
}
