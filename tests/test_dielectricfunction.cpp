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
#include "DielectricFunction.hpp"
#include "Material.h"
#include "Options.h"
#include "bz_mesh.hpp"
#include "bz_meshfile.hpp"
#include "doctest/doctest.h"

TEST_CASE("Epsilon_Si") {
    // Eigen::MatrixXcd Rmat = Eigen::MatrixXcd::Random(139,139);
    // Eigen::MatrixXcd Adj = Rmat.adjoint();
    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es;

    // es.compute(Rmat);

    // for (int i=0; i < es.eigenvalues().size(); ++i) {
    //     auto eigenvect = es.eigenvectors().col(i);
    //     std::cout << "Norm eigenvect ::" << i << ": " << eigenvect.norm() << std::endl;
    // }
    // exit(0);

    EmpiricalPseudopotential::Materials materials;
    const std::string                   file_material_parameters = std::string(CMAKE_SOURCE_DIR) + "/parameter_files/materials-local.yaml";
    materials.load_material_parameters(file_material_parameters);
    EmpiricalPseudopotential::Material current_material = materials.materials.at("Si");

    EmpiricalPseudopotential::BandStructure band_structure{};
    const std::size_t                       nb_bands           = 18;
    const std::size_t                       nearest_neightbors = 10;
    const bool                              non_local_corr     = false;
    band_structure.Initialize(current_material, nb_bands, {}, nearest_neightbors, non_local_corr);

    EmpiricalPseudopotential::DielectricFunction MyDielectricFunc(current_material, band_structure.get_basis_vectors(), nb_bands);
    // const std::size_t                            nb_kpoints = 2000;
    // MyDielectricFunc.generate_k_points_random(nb_kpoints);

    std::size_t Nxyz              = 100;
    bool        irreducible_wedge = true;
    MyDielectricFunc.generate_k_points_grid(Nxyz, Nxyz, Nxyz, 0.0, irreducible_wedge);
    std::cout << "Number of kpoints in the irreducible wedge: " << MyDielectricFunc.get_kpoints().size() << std::endl;
    MyDielectricFunc.export_kpoints("TestKpoints.csv");

    int    nb_threads   = 32;
    double eta_smearing = 2.0e-2;

    double              q_xyz      = 1.0e-12;
    double              dqx        = 0.005;
    const double        min_energy = 0.0;
    const double        max_energy = 20.0;
    const double        d_energy   = 0.01;
    std::vector<double> list_energy;
    for (double energy = min_energy; energy <= max_energy; energy += d_energy) {
        list_energy.push_back(energy);
    }

    for (double qx = q_xyz; qx <= 1.0+dqx; qx += dqx) {
        std::cout << "q_x : " << qx << std::endl;

        Vector3D<double> q_vect{qx, 0, 0};
        std::vector<double> list_epsilon = MyDielectricFunc.compute_dielectric_function(q_vect, list_energy, eta_smearing, nb_threads);
        std::string filename = std::string("Experiment_qx/") + "epsilon_Smearing" + std::to_string(eta_smearing) + "_Qx" +
                               std::to_string(qx) + "Nxyz" + std::to_string(Nxyz) + ".csv";
        std::ofstream file_dielectric_function(filename);
        file_dielectric_function << "energy,epsilon" << std::endl;
        for (std::size_t i = 0; i < list_energy.size(); ++i) {
            file_dielectric_function << list_energy[i] << "," << list_epsilon[i] << std::endl;
        }
        file_dielectric_function.close();
        const std::string python_plot_band_structure_script = std::string(CMAKE_SOURCE_DIR) + "/python/plots/plot_eps_vs_energy.py";
        std::string       python_call                       = "python3 " + python_plot_band_structure_script + " --filename " + filename;
        bool              call_python_plot                  = false;
        // bool call_python_plot = true;
        if (call_python_plot) {
            std::cout << "Executing: " << python_call << std::endl;
            int succes_plot = system(python_call.c_str());
            std::cout << "Succes plot: " << succes_plot << std::endl;
        }
    }
}
