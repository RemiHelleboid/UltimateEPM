/**
 * @file test_epsilon.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-11-16
 *
 * @copyright Copyright (c) 2022
 *
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

#include "BandStructure.h"
#include "Material.h"
#include "Options.h"
#include "bz_mesh.hpp"
#include "bz_meshfile.hpp"
// #include "dielectric_function.hpp"
#include "doctest/doctest.h"

TEST_CASE("Epsilon_Si") {
    EmpiricalPseudopotential::Materials materials;
    // const std::string                   file_material_parameters = std::string(CMAKE_SOURCE_DIR) + "/parameter_files/materials.yaml";
    // materials.load_material_parameters(file_material_parameters);
    // EmpiricalPseudopotential::Material current_material     = materials.materials.at("Si");
    // const std::string                  mesh_band_input_file = std::string(CMAKE_SOURCE_DIR) + "/examples/bz_mesh_1_50e-3.msh_1.msh";
    // bz_mesh::MeshBZ                    MeshGeometry(current_material);
    // MeshGeometry.read_mesh_geometry_from_msh_file(mesh_band_input_file, true    );

    // bz_mesh::DielectricFunction MyDielectricFunction(std::make_unique<bz_mesh::MeshBZ>(MeshGeometry), current_material);

    // const bz_mesh::vector3 q_small = {0.01, 0.01, 0.01};
    // int nb_energies = 100;
    // double energy_min = 0.0;
    // double energy_max = 10.0;
    // std::vector<double> list_energies;
    // for (int i = 0; i < nb_energies; i++) {
    //     list_energies.push_back(energy_min + i * (energy_max - energy_min) / (nb_energies - 1));
    // }
    // std::vector<double> list_epsilon_real;
    // for (auto energy : list_energies) {
    //     list_epsilon_real.push_back(MyDielectricFunction.compute_epsilon_real(q_small, energy));
    //     std::cout << "energy = " << energy << " epsilon_real = " << list_epsilon_real.back() << std::endl; 
    // }
    // for (int i = 0; i < nb_energies; i++) {
    //     std::cout << list_energies[i] << "," << list_epsilon_real[i] << std::endl;
    // }
    
}
