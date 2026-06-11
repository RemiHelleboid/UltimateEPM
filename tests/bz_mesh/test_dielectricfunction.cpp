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
#include "dielectric_mesh.hpp"
#include "doctest/doctest.h"

namespace {
class SyntheticDielectricMesh : public uepm::mesh_bz::DielectricMesh {
 public:
    void initialize() {
        m_list_vertices = {
            uepm::mesh_bz::Vertex(0, 0.0, 0.0, 0.0),
            uepm::mesh_bz::Vertex(1, 1.0, 0.0, 0.0),
            uepm::mesh_bz::Vertex(2, 0.0, 1.0, 0.0),
            uepm::mesh_bz::Vertex(3, 0.0, 0.0, 1.0),
        };
        std::array<uepm::mesh_bz::Vertex*, 4> vertices = {
            &m_list_vertices[0], &m_list_vertices[1], &m_list_vertices[2], &m_list_vertices[3]};
        m_list_tetrahedra.emplace_back(0, vertices);
        m_energies = {1.0, 3.0};
        m_dielectric_function = {
            {{0.0, 0.0}, {10.0, 20.0}},
            {{1.0, 2.0}, {11.0, 22.0}},
            {{2.0, 4.0}, {12.0, 24.0}},
            {{3.0, 6.0}, {13.0, 26.0}},
        };
        build_search_tree();
    }
};
}  // namespace

TEST_CASE("dielectric interpolation uses node-major storage and clamps energy") {
    SyntheticDielectricMesh mesh;
    mesh.initialize();
    const uepm::mesh_bz::vector3 centroid(0.25, 0.25, 0.25);

    const auto below = mesh.interpolate_dielectric_function(centroid, 0.0);
    CHECK(below.real() == doctest::Approx(1.5));
    CHECK(below.imag() == doctest::Approx(3.0));

    const auto middle = mesh.interpolate_dielectric_function(centroid, 2.0);
    CHECK(middle.real() == doctest::Approx(6.5));
    CHECK(middle.imag() == doctest::Approx(13.0));

    const auto above = mesh.interpolate_dielectric_function(centroid, 4.0);
    CHECK(above.real() == doctest::Approx(11.5));
    CHECK(above.imag() == doctest::Approx(23.0));
}

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

    // uepm::pseudopotential::Materials materials;
    // const std::string                   file_material_parameters = std::string(PROJECT_SRC_DIR) +
    // "/parameter_files/materials-chel.yaml"; materials.load_material_parameters(file_material_parameters);
    // uepm::pseudopotential::Material current_material = materials.materials.at("Si");

    // uepm::pseudopotential::BandStructure band_structure{};
    // const std::size_t                       nb_bands           = 18;
    // const std::size_t                       nearest_neightbors = 10;
    // const bool                              non_local_corr     = false;
    // band_structure.Initialize(current_material, nb_bands, {}, nearest_neightbors, non_local_corr);

    // uepm::pseudopotential::DielectricFunction MyDielectricFunc(current_material, band_structure.get_basis_vectors(), nb_bands);
    // // const std::size_t                            nb_kpoints = 2000;
    // // MyDielectricFunc.generate_k_points_random(nb_kpoints);

    // std::size_t Nxyz              = 100;
    // bool        irreducible_wedge = true;
    // MyDielectricFunc.generate_k_points_grid(Nxyz, Nxyz, Nxyz, 0.0, irreducible_wedge);
    // std::cout << "Number of kpoints in the irreducible wedge: " << MyDielectricFunc.get_kpoints().size() << std::endl;
    // MyDielectricFunc.export_kpoints("TestKpoints.csv");

    // int    nb_threads   = 32;
    // double eta_smearing = 2.0e-2;

    // double              q_xyz      = 1.0e-12;
    // double              dqx        = 0.005;
    // const double        min_energy = 0.0;
    // const double        max_energy = 20.0;
    // const double        d_energy   = 0.01;
    // std::vector<double> list_energy;
    // for (double energy = min_energy; energy <= max_energy; energy += d_energy) {
    //     list_energy.push_back(energy);
    // }

    // for (double qx = q_xyz; qx <= 1.0+dqx; qx += dqx) {
    //     std::cout << "q_x : " << qx << std::endl;

    //     Vector3D<double> q_vect{qx, 0, 0};
    //     std::vector<double> list_epsilon = MyDielectricFunc.compute_dielectric_function(q_vect, list_energy, eta_smearing);
    //     std::string filename = std::string("Experiment_qx/") + "epsilon_Smearing" + std::to_string(eta_smearing) + "_Qx" +
    //                            std::to_string(qx) + "Nxyz" + std::to_string(Nxyz) + ".csv";
    //     std::ofstream file_dielectric_function(filename);
    //     file_dielectric_function << "energy,epsilon" << std::endl;
    //     for (std::size_t i = 0; i < list_energy.size(); ++i) {
    //         file_dielectric_function << list_energy[i] << "," << list_epsilon[i] << std::endl;
    //     }
    //     file_dielectric_function.close();
    //     const std::string python_plot_band_structure_script = std::string(PROJECT_SRC_DIR) + "/python/plots/plot_eps_vs_energy.py";
    //     std::string       python_call                       = "python3 " + python_plot_band_structure_script + " --filename " + filename;
    //     bool              call_python_plot                  = false;
    //     // bool call_python_plot = true;
    //     if (call_python_plot) {
    //         std::cout << "Executing: " << python_call << std::endl;
    //         int succes_plot = system(python_call.c_str());
    //         std::cout << "Succes plot: " << succes_plot << std::endl;
    //     }
    // }
}
