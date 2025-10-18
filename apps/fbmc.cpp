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
    TCLAP::ValueArg<int>         arg_nb_energies("e", "nenergy", "Number of energies to compute", false, 250, "int");
    TCLAP::ValueArg<int>         arg_nb_conduction_bands("c", "ncbands", "Number of conduction bands to consider", false, -1, "int");
    TCLAP::ValueArg<int>         arg_nb_valence_bands("v", "nvbands", "Number of valence bands to consider", false, -1, "int");
    TCLAP::ValueArg<int>         arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    TCLAP::SwitchArg plot_with_python("P", "plot", "Call a python script after the computation to plot the band structure.", false);
    TCLAP::SwitchArg plot_with_wedge("w", "wedge", "Consider only the irreducible wedge of the BZ.", false);
    cmd.add(plot_with_python);
    cmd.add(arg_mesh_file);
    cmd.add(arg_material);
    cmd.add(arg_nb_conduction_bands);
    cmd.add(arg_nb_valence_bands);
    cmd.add(arg_nb_energies);
    cmd.add(arg_nb_threads);
    cmd.add(plot_with_wedge);
    cmd.add(arg_phonon_file);

    cmd.parse(argc, argv);

    const std::string file_mesh              = arg_mesh_file.getValue();
    const std::string material_symbol        = arg_material.getValue();
    const std::string file_phonon_scattering = arg_phonon_file.getValue();
    int nb_threads                          = arg_nb_threads.getValue();
    int               nb_valence_bands       = 0;
    int               nb_conduction_bands    = 4;

    uepm::pseudopotential::Materials materials;
    const std::string                file_material_parameters = std::string(PROJECT_SRC_DIR) + "/parameter_files/materials-local.yaml";
    materials.load_material_parameters(file_material_parameters);
    uepm::pseudopotential::Material current_material = materials.materials.at(arg_material.getValue());

    uepm::mesh_bz::ElectronPhonon mesh(current_material);
    mesh.set_number_threads_mesh_ops(nb_threads);
    mesh.read_mesh_geometry_from_msh_file(file_mesh);
    mesh.build_search_tree();
    mesh.export_octree_to_vtu("octree.vtu");
    mesh.read_mesh_bands_from_msh_file(file_mesh, nb_conduction_bands);
    mesh.read_phonon_scattering_rates_from_file(file_phonon_scattering);

    uepm::fbmc::Bulk_environment bulk_env;
    bulk_env.m_temperature          = 300.0;
    bulk_env.m_electric_field       = {1e5, 0.0, 0.0};
    bulk_env.m_doping_concentration = 1e10;

    uepm::fbmc::Simulation_parameters sim_params;
    sim_params.m_simulation_time  = 2000e-12;
    sim_params.m_export_frequency = 1;

    uepm::fbmc::Single_particle_simulation sim(&mesh, bulk_env, sim_params);
    sim.run_simulation();

    std::string timestamp = std::to_string(std::time(nullptr));
    std::string filename  = "particle_history_" + timestamp + ".csv";
    sim.export_history(filename);

    return 0;
}
