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
    uepm::mesh::vertex                     V1(0, 0.0, 0.0, 0.0);
    uepm::mesh::vertex                     V2(1, 1.0, 0.0, 0.0);
    uepm::mesh::vertex                     V3(2, 0.0, 1.0, 0.0);
    uepm::mesh::vertex                     V4(3, 0.0, 0.0, 1.0);
    std::shared_ptr<uepm::mesh::element3d> sp_reference_element =
        std::make_shared<uepm::mesh::element3d>(&V1, &V2, &V3, &V4);
    const Eigen::Matrix4d MatrixElementRef =
        uepm::fem::FiniteElementP1System3d::compute_elementary_stiffness_matrix(sp_reference_element);
    Eigen::Matrix4d       THEORETICAL_MATRIX{{3.0, -1.0, -1.0, -1.0},
                                             {-1.0, 1.0, 0.0, 0.0},
                                             {-1.0, 0.0, 1.0, 0.0},
                                             {-1.0, 0.0, 0.0, 1.0}};
    const Eigen::Matrix4d DIFFERENCE_MATRIX = (1.0 / 6.0) * THEORETICAL_MATRIX - MatrixElementRef;
    std::cout << THEORETICAL_MATRIX << std::endl << std::endl << std::endl;
    std::cout << MatrixElementRef << std::endl;
    std::cout << "VOLUME : " << sp_reference_element->get_measure() << std::endl;
    CHECK(DIFFERENCE_MATRIX.squaredNorm() < 1e-9);
}

TEST_CASE("Testing Poisson 3d on a unit sphere.") {
    static const std::string file_input_test_msh = PROJECT_SRC_DIR + std::string("/tests/test_data/sphere.msh");
    uepm::file::msh_file     fileMSH(file_input_test_msh);
    fileMSH.read_mesh();
    uepm::mesh::mesh*                  p_mesh = fileMSH.get_p_mesh();
    uepm::fem::FiniteElementP1System3d MyPoissonTest(p_mesh, p_mesh->get_nb_vertices());
    MyPoissonTest.compute_stiffness_matrix();
    MyPoissonTest.compute_second_member(6.0);
    MyPoissonTest.apply_dirichlet_condition("boundary", 0.0);
    MyPoissonTest.decompose_matrix();
    MyPoissonTest.solve_system();
    MyPoissonTest.export_solution_csv("Poisson3d_sphere.csv");
    MyPoissonTest.add_solution_to_mesh_functions("Poisson_Solution");

    auto   solution     = MyPoissonTest.get_solution();
    double min_solution = *std::min_element(solution.begin(), solution.end());
    double max_solution = *std::max_element(solution.begin(), solution.end());

    const double min_test_circle_sphere = 0.0;
    const double max_test_circle_sphere = 1.0;

    CHECK_EQ(min_solution, doctest::Approx(min_test_circle_sphere).epsilon(1e-2));
    CHECK_EQ(max_solution, doctest::Approx(max_test_circle_sphere).epsilon(1e-2));

    const std::string FileName = "result_poisson_3d_sphere_dirichlet.vtk";
    uepm::file::export_as_vtk(*(p_mesh), FileName, {}, {}, true);
}

TEST_CASE("Testing Poisson 3d on a unit cube.") {
    static const std::string file_input_test_msh = PROJECT_SRC_DIR + std::string("/tests/test_data/cube.msh");
    uepm::file::msh_file           fileMSH(file_input_test_msh);
    fileMSH.read_mesh();
    uepm::mesh::mesh*                  p_mesh = fileMSH.get_p_mesh();
    uepm::fem::FiniteElementP1System3d MyPoissonTest(p_mesh, p_mesh->get_nb_vertices());
    MyPoissonTest.compute_stiffness_matrix();

    auto func3d = [](double x, double y, double z) {return 6.0;};
    MyPoissonTest.compute_second_member(func3d);
    MyPoissonTest.apply_dirichlet_condition("x_min", 0.0);
    MyPoissonTest.apply_dirichlet_condition("x_max", 0.0);
    MyPoissonTest.apply_dirichlet_condition("y_min", 0.0);
    MyPoissonTest.apply_dirichlet_condition("y_max", 0.0);
    MyPoissonTest.apply_dirichlet_condition("z_min", 0.0);
    MyPoissonTest.apply_dirichlet_condition("z_max", 0.0);
    MyPoissonTest.decompose_matrix();
    MyPoissonTest.solve_system();
    MyPoissonTest.export_solution_csv("Poisson3d_cube.csv");
    auto   solution     = MyPoissonTest.get_solution();
    double min_solution = *std::min_element(solution.begin(), solution.end());
    double max_solution = *std::max_element(solution.begin(), solution.end());

    const double min_test_circle_sphere = 0.0;
    const double max_test_circle_sphere = 0.336626;

    CHECK_EQ(min_solution, doctest::Approx(min_test_circle_sphere).epsilon(1e-2));
    CHECK_EQ(max_solution, doctest::Approx(max_test_circle_sphere).epsilon(1e-2));

    MyPoissonTest.add_solution_to_mesh_functions("Poisson_Solution");

    const std::string FileName = "result_poisson_3d_cube_dirichlet.vtk";
    uepm::file::export_as_vtk(*(p_mesh), FileName, {}, {}, true);
}

TEST_CASE("Testing Poisson 3d on a unit cube with Neuman BC.") {
    static const std::string file_input_test_msh = PROJECT_SRC_DIR + std::string("/tests/test_data/cube.msh");
    uepm::file::msh_file           fileMSH(file_input_test_msh);
    fileMSH.read_mesh();
    uepm::mesh::mesh*                  p_mesh = fileMSH.get_p_mesh();
    uepm::fem::FiniteElementP1System3d MyPoissonTest(p_mesh, p_mesh->get_nb_vertices());
    MyPoissonTest.compute_stiffness_matrix();
    MyPoissonTest.compute_second_member(0.0);
    MyPoissonTest.apply_neuman_condition("x_min", 5.0);
    MyPoissonTest.apply_neuman_condition("x_max", 5.0);
    MyPoissonTest.apply_dirichlet_condition("y_min", 0.0);
    MyPoissonTest.apply_dirichlet_condition("y_max", 0.0);
    MyPoissonTest.apply_dirichlet_condition("z_min", 0.0);
    MyPoissonTest.apply_dirichlet_condition("z_max", 0.0);
    MyPoissonTest.decompose_matrix();
    MyPoissonTest.solve_system();
    MyPoissonTest.export_solution_csv("Poisson3d_neuman_cube.csv");
    auto   solution     = MyPoissonTest.get_solution();
    double min_solution = *std::min_element(solution.begin(), solution.end());
    double max_solution = *std::max_element(solution.begin(), solution.end());

    const double min_test_neuman = 0;
    const double max_test_neuman = 1.47619;

    CHECK_EQ(min_solution, doctest::Approx(min_test_neuman).epsilon(1e-2));
    CHECK_EQ(max_solution, doctest::Approx(max_test_neuman).epsilon(1e-2));
    MyPoissonTest.add_solution_to_mesh_functions("Poisson_Solution");

    const std::string FileName = "result_poisson_3d_cube_dirichlet.vtk";
    uepm::file::export_as_vtk(*(p_mesh), FileName, {}, {}, true);
}
