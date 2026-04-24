/**
 * @file test_poisson_2d.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-06
 *
 * @copyright Copyright (c) 2021
 *
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest/doctest.h"
#include "msh_file.hpp"
#include "poisson_solver_2d.hpp"
#include "vtkWriter.hpp"



TEST_CASE("Testing Poisson 2d on a PN Junction.") {

    std::string material_file = CMAKE_SOURCE_DIR + std::string("/data/materials.yaml"); 
    std::cout << "Material file : " << material_file << std::endl;
    physic::material::list_materials list_of_materials;
    list_of_materials.load_materials_from_file(material_file);

    static const std::string file_input_test_msh = CMAKE_SOURCE_DIR + std::string("/tests/test_data/2D_pn_diode_5V.msh");
    file::msh_file           fileMSH(file_input_test_msh);
    fileMSH.read_mesh();
    fileMSH.read_states();
    mesh::mesh* p_mesh     = fileMSH.get_p_mesh();
    std::size_t nbVertices = p_mesh->get_nb_vertices();

    std::cout << "Start building Poisson system ..." << std::endl;
    fem::poisson_solver_2d MyPoissonSolver(p_mesh, p_mesh->get_nb_vertices(), list_of_materials);
    MyPoissonSolver.compute_stiffness_matrix();
    MyPoissonSolver.update_second_member();
    
    MyPoissonSolver.apply_dirichlet_condition("kathode", 5.36);
    MyPoissonSolver.apply_dirichlet_condition("anode", -0.432793);
    // MyPoissonSolver.apply_dirichlet_condition("kathode", 0.355292);
    std::cout << "Start computing poisson solution ..." << std::endl;

    MyPoissonSolver.decompose_matrix();
    MyPoissonSolver.solve_system();
    MyPoissonSolver.export_solution_csv("5_PN_JUNCTION_SOLUTION.csv");
    MyPoissonSolver.add_solution_to_mesh_functions("Armin_Solution");
    fileMSH.export_as_msh("TEST_POISSON_DIODE_PN_5V.msh", {}, 1);
    const std::string FileName = "TEST_POISSON_DIODE_PN_5V.vtk";
    file::export_as_vtk(*(p_mesh), FileName, {}, {}, true);

    // p_mesh->export_all_vertices_data_to_csv("POISSON_PN_CSV.csv");

    auto   solution     = MyPoissonSolver.get_solution();
    double min_solution = *std::min_element(solution.begin(), solution.end());
    double max_solution = *std::max_element(solution.begin(), solution.end());

    std::cout << "MIN SOL = " << min_solution << std::endl;
    std::cout << "MAX SOL = " << max_solution << std::endl;

    const double min_test_si_ge = -0.765151;
    const double max_test_si_ge = 5.28987;

    // CHECK(min_solution == doctest::Approx(min_test_si_ge));
    // CHECK_EQ(max_solution, doctest::Approx(max_test_si_ge));
}


