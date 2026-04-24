/**
 * @file test_fem_3d.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-09
 *
 * @copyright Copyright (c) 2021
 *
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <cmath>

#include "doctest/doctest.h"
#include "finite_element3d.hpp"
#include "materials.hpp"
#include "msh_file.hpp"
#include "vtkWriter.hpp"


TEST_CASE("Test the stiffness elementary matrix on ref element") {
    // Construct the element.
    mesh::vertex                     V1(0, 0.0, 0.0, 0.0);
    mesh::vertex                     V2(1, 1.0, 0.0, 0.0);
    mesh::vertex                     V3(2, 0.0, 1.0, 0.0);
    mesh::vertex                     V4(3, 0.0, 0.0, 1.0);
    std::shared_ptr<mesh::element3d> sp_reference_element = std::make_shared<mesh::element3d>(&V1, &V2, &V3, &V4);
    const Eigen::Matrix4d MatrixElementRef = fem::FiniteElementP1System3d::compute_elementary_stiffness_matrix(sp_reference_element);
    Eigen::Matrix4d THEORETICAL_MATRIX{{3.0, -1.0, -1.0, -1.0}, {-1.0, 1.0, 0.0, 0.0}, {-1.0, 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0, 1.0}};
    THEORETICAL_MATRIX = THEORETICAL_MATRIX * -1.0; 
    const Eigen::Matrix4d DIFFERENCE_MATRIX = (1.0 / 6.0) * THEORETICAL_MATRIX - MatrixElementRef;
    std::cout << THEORETICAL_MATRIX << std::endl << std::endl << std::endl;
    std::cout << MatrixElementRef << std::endl;
    std::cout << "VOLUME : " << sp_reference_element->get_measure() << std::endl;
    CHECK(DIFFERENCE_MATRIX.squaredNorm() < 1e-9);
}

TEST_CASE("Testing Poisson 3d on a unit sphere.") {
    static const std::string file_input_test_msh = CMAKE_SOURCE_DIR + std::string("/example/data/sphere_r1.msh");
    file::msh_file           fileMSH(file_input_test_msh);
    fileMSH.read_mesh();
    mesh::mesh*                  p_mesh = fileMSH.get_p_mesh();
    fem::FiniteElementP1System3d MyPoissonTest(p_mesh, p_mesh->get_nb_vertices());
    MyPoissonTest.compute_stiffness_matrix();
    MyPoissonTest.compute_second_member(-6.0);
    MyPoissonTest.apply_dirichlet_condition("Contact_1", 0.0);
    MyPoissonTest.decompose_matrix();
    MyPoissonTest.solve_system();
    MyPoissonTest.export_solution_csv("Poisson3d_sphere.csv");
    MyPoissonTest.add_solution_to_mesh_functions("Armin_Solution");

    auto   solution     = MyPoissonTest.get_solution();
    double min_solution = *std::min_element(solution.begin(), solution.end());
    double max_solution = *std::max_element(solution.begin(), solution.end());

    const double min_test_circle_sphere = 0.0;
    const double max_test_circle_sphere = 1.0;

    // 1% error accepted
    CHECK_EQ(min_solution, doctest::Approx(min_test_circle_sphere).epsilon(0.01));
    CHECK_EQ(max_solution, doctest::Approx(max_test_circle_sphere).epsilon(0.01));

    const std::string FileName = "TEST_POISSON_3D_SPHERE.vtk";
    file::export_as_vtk(*(p_mesh), FileName, {}, {}, true);
}

// TEST_CASE("Testing Poisson 3d on a unit cube.") {
//     static const std::string file_input_test_msh = CMAKE_SOURCE_DIR + std::string("/example/data/simple_cube.msh");
//     file::msh_file           fileMSH(file_input_test_msh);
//     fileMSH.read_mesh();
//     mesh::mesh*                  p_mesh = fileMSH.get_p_mesh();
//     fem::FiniteElementP1System3d MyPoissonTest(p_mesh, p_mesh->get_nb_vertices());
//     MyPoissonTest.compute_stiffness_matrix();

//     auto func3d = [](double x, double y, double z) {return (2.0*M_PI) * (2.0*M_PI) * std::sin((2.0*M_PI)*x); };
//     MyPoissonTest.compute_second_member(func3d);
//     MyPoissonTest.apply_dirichlet_condition("Contact_1", 0.0);
//     MyPoissonTest.apply_dirichlet_condition("Contact_2", 0.0);
//     MyPoissonTest.apply_dirichlet_condition("Contact_3", 0.0);
//     MyPoissonTest.apply_dirichlet_condition("Contact_4", 0.0);
//     MyPoissonTest.apply_dirichlet_condition("Contact_5", 0.0);
//     MyPoissonTest.apply_dirichlet_condition("Contact_6", 0.0);
//     MyPoissonTest.decompose_matrix();
//     MyPoissonTest.solve_system();
//     MyPoissonTest.export_solution_csv("Poisson3d_cube.csv");
//     auto   solution     = MyPoissonTest.get_solution();
//     double min_solution = *std::min_element(solution.begin(), solution.end());
//     double max_solution = *std::max_element(solution.begin(), solution.end());

//     const double min_test_circle_sphere = 0.0;
//     const double max_test_circle_sphere = 0.0558769;

//     CHECK_EQ(min_solution, doctest::Approx(min_test_circle_sphere).epsilon(0.01));
//     CHECK_EQ(max_solution, doctest::Approx(max_test_circle_sphere).epsilon(0.01));

//     MyPoissonTest.add_solution_to_mesh_functions("ArminxPoissonxSolution");

//     const std::string FileName = "test_3d_FEM.vtk";
//     file::export_as_vtk(*(p_mesh), FileName, {}, {}, true);
// }

// TEST_CASE("Testing Poisson 3d on a unit cube with Neuman BC.") {
//     static const std::string file_input_test_msh = CMAKE_SOURCE_DIR + std::string("/example/data/simple_cube.msh");
//     file::msh_file           fileMSH(file_input_test_msh);
//     fileMSH.read_mesh();
//     mesh::mesh*                  p_mesh = fileMSH.get_p_mesh();
//     fem::FiniteElementP1System3d MyPoissonTest(p_mesh, p_mesh->get_nb_vertices());
//     MyPoissonTest.compute_stiffness_matrix();
//     MyPoissonTest.compute_second_member(1.0);
//     MyPoissonTest.apply_neuman_condition("Contact_1", -5.0);
//     MyPoissonTest.apply_neuman_condition("Contact_2", 5.0);
//     MyPoissonTest.apply_dirichlet_condition("Contact_3", 0.0);
//     MyPoissonTest.apply_dirichlet_condition("Contact_4", 0.0);
//     MyPoissonTest.apply_dirichlet_condition("Contact_5", 0.0);
//     MyPoissonTest.apply_dirichlet_condition("Contact_6", 0.0);
//     MyPoissonTest.decompose_matrix();
//     MyPoissonTest.solve_system();
//     MyPoissonTest.export_solution_csv("Poisson3d_neuman_cube.csv");
//     auto   solution     = MyPoissonTest.get_solution();
//     double min_solution = *std::min_element(solution.begin(), solution.end());
//     double max_solution = *std::max_element(solution.begin(), solution.end());

//     const double min_test_neuman = -1.32943;
//     const double max_test_neuman = 1.47619;

//     CHECK_EQ(min_solution, doctest::Approx(min_test_neuman).epsilon(0.01));
//     CHECK_EQ(max_solution, doctest::Approx(max_test_neuman).epsilon(0.01));
// }
