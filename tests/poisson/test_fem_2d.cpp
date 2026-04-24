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
    mesh::vertex V1(0, 0.0, 0.0, 0.0);
    mesh::vertex V2(1, 1.0, 0.0, 0.0);
    mesh::vertex V3(2, 0.0, 1.0, 0.0);
    std::shared_ptr<mesh::element2d> sp_reference_element = std::make_shared<mesh::element2d>(&V1, &V2, &V3);
    const Eigen::Matrix3d MatrixElementRef = fem::FiniteElementP1System2d::compute_elementary_stiffness_matrix(sp_reference_element);
    const Eigen::Matrix3d THEORETICAL_MATRIX{{2.0, -1.0, -1.0}, {-1.0, 1.0, 0.0}, {-1.0, 0.0, 1.0}};
    const Eigen::Matrix3d DIFFERENCE_MATRIX = 0.5 * THEORETICAL_MATRIX - MatrixElementRef;
    CHECK(DIFFERENCE_MATRIX.squaredNorm() < 1e-9);
}

TEST_CASE("Testing Poisson 2d on a unit circle.") {
    static const std::string file_input_test_msh = CMAKE_SOURCE_DIR + std::string("/example/data/circle_r1.msh");
    file::msh_file           fileMSH(file_input_test_msh);
    fileMSH.read_mesh();
    mesh::mesh*                  p_mesh = fileMSH.get_p_mesh();
    fem::FiniteElementP1System2d MyPoissonTest(p_mesh, p_mesh->get_nb_vertices());
    MyPoissonTest.compute_stiffness_matrix();
    MyPoissonTest.compute_second_member(4.0);
    MyPoissonTest.apply_dirichlet_condition("Contact_1", 0.0);
    MyPoissonTest.decompose_matrix();
    MyPoissonTest.solve_system();
    MyPoissonTest.export_solution_csv("Poisson2d_circle.csv");
    MyPoissonTest.add_solution_to_mesh_functions("Armin_Solution");

    auto   solution     = MyPoissonTest.get_solution();
    double min_solution = *std::min_element(solution.begin(), solution.end());
    double max_solution = *std::max_element(solution.begin(), solution.end());

    const double min_test_circle_poisson = 0.0;
    const double max_test_circle_poisson = 1.0;

    CHECK(min_solution == doctest::Approx(min_test_circle_poisson));
    CHECK_EQ(max_solution, doctest::Approx(max_test_circle_poisson));

    const std::string FileName = "TEST_POISSON_2D_CIRCLE.vtk";
    file::export_as_vtk(*(p_mesh), FileName, {}, {}, true);
}

TEST_CASE("Testing Poisson 2d on a unit square.") {
    static const std::string file_input_test_msh = CMAKE_SOURCE_DIR + std::string("/example/data/square_test.msh");
    file::msh_file           fileMSH(file_input_test_msh);
    fileMSH.read_mesh();
    mesh::mesh*                  p_mesh = fileMSH.get_p_mesh();
    fem::FiniteElementP1System2d MyPoissonTest(p_mesh, p_mesh->get_nb_vertices());
    MyPoissonTest.compute_stiffness_matrix();
    MyPoissonTest.compute_second_member(1.0);
    MyPoissonTest.apply_dirichlet_condition("Contact_1", 0.0);
    MyPoissonTest.apply_dirichlet_condition("Contact_2", 0.0);
    MyPoissonTest.apply_dirichlet_condition("Contact_3", 0.0);
    MyPoissonTest.apply_dirichlet_condition("Contact_4", 0.0);
    MyPoissonTest.decompose_matrix();
    MyPoissonTest.solve_system();
    MyPoissonTest.export_solution_csv("Poisson2d_square.csv");
    MyPoissonTest.add_solution_to_mesh_functions("Armin_Solution");

    auto   solution     = MyPoissonTest.get_solution();
    double min_solution = *std::min_element(solution.begin(), solution.end());
    double max_solution = *std::max_element(solution.begin(), solution.end());

    const double min_test_square_poisson = 0.0;
    const double max_test_square_poisson = 0.0736688;

    CHECK(min_solution == doctest::Approx(min_test_square_poisson));
    CHECK_EQ(max_solution, doctest::Approx(max_test_square_poisson));

    const std::string FileName = "TEST_POISSON_2D_SQUARE.vtk";
    file::export_as_vtk(*(p_mesh), FileName, {}, {}, true);
}

TEST_CASE("Testing Poisson 2d on a unit square with Neuman BC.") {
    static const std::string file_input_test_msh = CMAKE_SOURCE_DIR + std::string("/example/data/square_1234.msh");
    file::msh_file           fileMSH(file_input_test_msh);
    fileMSH.read_mesh();
    mesh::mesh* p_mesh = fileMSH.get_p_mesh();

    fem::FiniteElementP1System2d MyPoissonTest(p_mesh, p_mesh->get_nb_vertices());
    MyPoissonTest.compute_stiffness_matrix();
    MyPoissonTest.compute_second_member(1.0);
    MyPoissonTest.apply_dirichlet_condition("Contact_1", 0.0);
    MyPoissonTest.apply_dirichlet_condition("Contact_3", 0.0);
    MyPoissonTest.apply_neuman_condition("Contact_2", -5.0);
    MyPoissonTest.apply_neuman_condition("Contact_4", 5.0);
    MyPoissonTest.decompose_matrix();
    MyPoissonTest.solve_system();
    MyPoissonTest.export_solution_csv("Poisson2d_square_Neuman.csv");
    MyPoissonTest.add_solution_to_mesh_functions("Armin_Solution");

    auto   solution     = MyPoissonTest.get_solution();
    double min_solution = *std::min_element(solution.begin(), solution.end());
    double max_solution = *std::max_element(solution.begin(), solution.end());

    const double min_test_square_poisson = -1.56298;
    const double max_test_square_poisson = 1.81292;

    CHECK(min_solution == doctest::Approx(min_test_square_poisson));
    CHECK_EQ(max_solution, doctest::Approx(max_test_square_poisson));

    const std::string FileName = "TEST_POISSON_2D_SQUARE_NEUMAN.vtk";
    file::export_as_vtk(*(p_mesh), FileName, {}, {}, true);
}

TEST_CASE("Testing Poisson 2d with a strong arctan profile.") {
    static const std::string file_input_test_msh = CMAKE_SOURCE_DIR + std::string("/example/data/square_1234.msh");
    file::msh_file           fileMSH(file_input_test_msh);
    fileMSH.read_mesh();
    mesh::mesh* p_mesh = fileMSH.get_p_mesh();

    fem::FiniteElementP1System2d MyPoissonTest(p_mesh, p_mesh->get_nb_vertices());
    MyPoissonTest.compute_stiffness_matrix();

    auto arctan_profile = [&](double x, double y) {return (- 1000.0 * atan(x-0.5));};
    MyPoissonTest.compute_second_member(arctan_profile);
    MyPoissonTest.apply_dirichlet_condition("Contact_2", 20);
    MyPoissonTest.apply_dirichlet_condition("Contact_4", 0);
    MyPoissonTest.decompose_matrix();
    MyPoissonTest.solve_system();
    MyPoissonTest.export_solution_csv("Poisson2d_Arctan.csv");
    MyPoissonTest.add_solution_to_mesh_functions("Armin_Solution");

    auto   solution     = MyPoissonTest.get_solution();
    double min_solution = *std::min_element(solution.begin(), solution.end());
    double max_solution = *std::max_element(solution.begin(), solution.end());

    const double min_test_square_poisson = -1.56298;
    const double max_test_square_poisson = 1.81292;



    // CHECK(min_solution == doctest::Approx(min_test_square_poisson));
    // CHECK_EQ(max_solution, doctest::Approx(max_test_square_poisson));

    const std::string FileName = "TEST_POISSON_2D_ARCTAN.vtk";
    file::export_as_vtk(*(p_mesh), FileName, {}, {}, true);
}