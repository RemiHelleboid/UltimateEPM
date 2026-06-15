/**
 * @file test_fem_2d.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-04
 *
 * @copyright Copyright (c) 2021
 *
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <cmath>

#include "doctest/doctest.h"
#include "finite_element2d.hpp"
#include "materials.hpp"
#include "msh_file.hpp"
#include "vtkWriter.hpp"

TEST_CASE("Test the stiffness elementary matrix on ref element") {
    // Construct the element.
    uepm::mesh::vertex                     V1(0, 0.0, 0.0, 0.0);
    uepm::mesh::vertex                     V2(1, 1.0, 0.0, 0.0);
    uepm::mesh::vertex                     V3(2, 0.0, 1.0, 0.0);
    std::shared_ptr<uepm::mesh::element2d> sp_reference_element =
        std::make_shared<uepm::mesh::element2d>(&V1, &V2, &V3);
    const Eigen::Matrix3d MatrixElementRef =
        uepm::fem::FiniteElementP1System2d::compute_elementary_stiffness_matrix(sp_reference_element);
    const Eigen::Matrix3d THEORETICAL_MATRIX{{2.0, -1.0, -1.0}, {-1.0, 1.0, 0.0}, {-1.0, 0.0, 1.0}};
    const Eigen::Matrix3d DIFFERENCE_MATRIX = 0.5 * THEORETICAL_MATRIX - MatrixElementRef;
    CHECK(DIFFERENCE_MATRIX.squaredNorm() < 1e-9);
}

TEST_CASE("Testing Poisson 2d on a unit disk.") {
    static const std::string file_input_test_msh = PROJECT_SRC_DIR + std::string("/tests/test_data/disk.msh");
    uepm::file::msh_file     fileMSH(file_input_test_msh);
    fileMSH.read_mesh();
    uepm::mesh::mesh*                  p_mesh = fileMSH.get_p_mesh();
    uepm::fem::FiniteElementP1System2d MyPoissonTest(p_mesh, p_mesh->get_nb_vertices());
    MyPoissonTest.compute_stiffness_matrix();
    MyPoissonTest.compute_second_member(4.0);
    MyPoissonTest.apply_dirichlet_condition("edge", 0.0);
    MyPoissonTest.decompose_matrix();
    MyPoissonTest.solve_system();
    MyPoissonTest.export_solution_csv("Poisson2d_disk.csv");
    MyPoissonTest.add_solution_to_mesh_functions("Poisson_Solution");

    auto   solution     = MyPoissonTest.get_solution();
    double min_solution = *std::min_element(solution.begin(), solution.end());
    double max_solution = *std::max_element(solution.begin(), solution.end());

    const double min_test_disk_poisson = 0.0;
    const double max_test_disk_poisson = 1.0;

    CHECK(min_solution == doctest::Approx(min_test_disk_poisson));
    CHECK_EQ(max_solution, doctest::Approx(max_test_disk_poisson).epsilon(1e-2));

    const std::string FileName = "result_poisson_2d_disk_dirichlet.vtk";
    uepm::file::export_as_vtk(*(p_mesh), FileName, {}, {}, true);
}

TEST_CASE("Testing Poisson 2d on a unit square.") {
    static const std::string file_input_test_msh = PROJECT_SRC_DIR + std::string("/tests/test_data/square.msh");
    uepm::file::msh_file     fileMSH(file_input_test_msh);
    fileMSH.read_mesh();
    uepm::mesh::mesh*                  p_mesh = fileMSH.get_p_mesh();
    uepm::fem::FiniteElementP1System2d MyPoissonTest(p_mesh, p_mesh->get_nb_vertices());
    MyPoissonTest.compute_stiffness_matrix();
    MyPoissonTest.compute_second_member(4.0);
    MyPoissonTest.apply_dirichlet_condition("edge_0", 0.0);
    MyPoissonTest.apply_dirichlet_condition("edge_1", 0.0);
    MyPoissonTest.apply_dirichlet_condition("edge_2", 0.0);
    MyPoissonTest.apply_dirichlet_condition("edge_3", 0.0);
    MyPoissonTest.decompose_matrix();
    MyPoissonTest.solve_system();
    MyPoissonTest.export_solution_csv("Poisson2d_square.csv");
    MyPoissonTest.add_solution_to_mesh_functions("Poisson_Solution");

    auto   solution     = MyPoissonTest.get_solution();
    double min_solution = *std::min_element(solution.begin(), solution.end());
    double max_solution = *std::max_element(solution.begin(), solution.end());

    const double min_test_square_poisson = 0.0;
    const double max_test_square_poisson = 0.294581;

    CHECK(min_solution == doctest::Approx(min_test_square_poisson).epsilon(1e-2));
    CHECK_EQ(max_solution, doctest::Approx(max_test_square_poisson).epsilon(1e-2));

    const std::string FileName = "result_poisson_2d_square_dirichlet.vtk";
    uepm::file::export_as_vtk(*(p_mesh), FileName, {}, {}, true);
}

TEST_CASE("Testing Poisson 2d on a unit square with Neuman BC.") {
    static const std::string file_input_test_msh = PROJECT_SRC_DIR + std::string("/tests/test_data/square.msh");
    uepm::file::msh_file     fileMSH(file_input_test_msh);
    fileMSH.read_mesh();
    uepm::mesh::mesh* p_mesh = fileMSH.get_p_mesh();

    uepm::fem::FiniteElementP1System2d MyPoissonTest(p_mesh, p_mesh->get_nb_vertices());
    MyPoissonTest.compute_stiffness_matrix();
    MyPoissonTest.compute_second_member(1.0);
    MyPoissonTest.apply_dirichlet_condition("edge_0", 0.0);
    MyPoissonTest.apply_dirichlet_condition("edge_2", 0.0);
    MyPoissonTest.apply_neuman_condition("edge_1", -5.0);
    MyPoissonTest.apply_neuman_condition("edge_3", 5.0);
    MyPoissonTest.decompose_matrix();
    MyPoissonTest.solve_system();
    MyPoissonTest.export_solution_csv("Poisson2d_square_Neuman.csv");
    MyPoissonTest.add_solution_to_mesh_functions("Poisson_Solution");

    auto   solution     = MyPoissonTest.get_solution();
    double min_solution = *std::min_element(solution.begin(), solution.end());
    double max_solution = *std::max_element(solution.begin(), solution.end());

    const double min_test_square_poisson = -1.56298;
    const double max_test_square_poisson = 1.81292;

    CHECK(min_solution == doctest::Approx(min_test_square_poisson).epsilon(1e-2));
    CHECK_EQ(max_solution, doctest::Approx(max_test_square_poisson).epsilon(1e-2));

    const std::string FileName = "result_poisson_2d_square_neuman.vtk";
    uepm::file::export_as_vtk(*(p_mesh), FileName, {}, {}, true);
}
TEST_CASE("Testing Poisson 2d linear patch test") {
    static const std::string file_input_test_msh = PROJECT_SRC_DIR + std::string("/tests/test_data/square.msh");

    uepm::file::msh_file fileMSH(file_input_test_msh);
    fileMSH.read_mesh();
    uepm::mesh::mesh* p_mesh = fileMSH.get_p_mesh();
    uepm::fem::FiniteElementP1System2d poisson(p_mesh, p_mesh->get_nb_vertices());
    poisson.compute_stiffness_matrix();
    auto zero_source = [](double, double) { return 0.0; };
    poisson.compute_second_member(zero_source);
    poisson.apply_dirichlet_condition("edge_3", 0.0);
    poisson.apply_dirichlet_condition("edge_1", 20.0);
    poisson.decompose_matrix();
    poisson.solve_system();
    const auto& solution = poisson.get_solution();

    double max_error = 0.0;

    for (std::size_t i = 0; i < p_mesh->get_nb_vertices(); ++i) {
        const auto *vertex = p_mesh->get_p_vertex(i);

        const double x     = vertex->x();
        const double exact = 20.0 * x;

        max_error = std::max(max_error, std::abs(solution[i] - exact));
    }

    CHECK(max_error < 1.0e-10);
}