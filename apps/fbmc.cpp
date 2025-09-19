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
#include "bz_mesh.hpp"
#include "bz_states.hpp"
#include "electron_phonon.hpp"
#include "single_part_fbmc.hpp"



int main(int argc, const char** argv) {
    std::cout << "Hello SinglePartFBMC!" << std::endl;

    const std::string file_mesh              = "mc_data/mesh_si_cb_1.msh";
    const std::string file_phonon_scattering = "mc_data/rates_all.csv";
    int               nb_conduction_bands    = 4;
    int               nb_valence_bands       = 0;

    const std::string mesh_band_input_file = "bz_si_cb_1.msh";
    const std::string phonon_file          = std::string(CMAKE_SOURCE_DIR) + "/parameter_files/phonon_kamakura.yaml";


    bz_mesh::MeshBZ mesh;
    mesh.read_mesh_geometry_from_msh_file(file_mesh);
    mesh.read_mesh_bands_from_msh_file(file_mesh, nb_valence_bands+nb_conduction_bands);
    mesh.read_phonon_scattering_rates_from_file(file_phonon_scattering);

    fbmc::Bulk_environment bulk_env;
    bulk_env.m_temperature          = 300.0;
    bulk_env.m_electric_field       = {1e5, 0.0, 0.0};
    bulk_env.m_doping_concentration = 1e10;

    fbmc::Simulation_parameters sim_params;
    sim_params.m_simulation_time  = 2000e-12;
    sim_params.m_export_frequency = 1;

    fbmc::Single_particle_simulation sim(&mesh, bulk_env, sim_params);
    sim.run_simulation();

    std::string timestamp = std::to_string(std::time(nullptr));
    std::string filename  = "particle_history_" + timestamp + ".csv";
    sim.export_history(filename);

    return 0;
}
