/**
 * @file test_poisson_3d.cpp
 * @author Rémi Helleboid (remi.helleboid@st.com)
 * @brief 
 * @version 0.1
 * @date 2021-10-19
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest/doctest.h"
#include "msh_file.hpp"
#include "poisson_solver_3d.hpp"
#include "materials.hpp"
#include "msh_file.hpp"
#include "vtkWriter.hpp"



TEST_CASE("Testing Poisson 3d on a PN Junction.") {
    std::string material_file = PROJECT_SRC_DIR + std::string("/data/materials.yaml");
    uepm::physic::material::list_materials list_of_materials;
    list_of_materials.load_materials_from_file(material_file);


    static const std::string file_input_test_STF = PROJECT_SRC_DIR + std::string("/tests/test_data/diode_test_admc.msh");
    uepm::file::msh_file           fileMSH(file_input_test_STF);
    fileMSH.read_mesh();
    fileMSH.read_states();
    uepm::mesh::mesh* p_mesh     = fileMSH.get_p_mesh();
    std::size_t nbVertices = p_mesh->get_nb_vertices();

    std::cout << "Start building Poisson system ..." << std::endl;
    uepm::fem::poisson_solver_3d MyPoissonSolver(p_mesh, p_mesh->get_nb_vertices(), list_of_materials);
    MyPoissonSolver.compute_stiffness_matrix();
    std::cout << "Start computing second member ..." << std::endl;
    MyPoissonSolver.update_second_member();
    
    // MyPoissonSolver.apply_dirichlet_condition("cathode", -0.36);
    MyPoissonSolver.apply_dirichlet_condition("anode", 0.49);
    std::cout << "Start computing poisson solution ..." << std::endl;

    MyPoissonSolver.decompose_matrix();
    MyPoissonSolver.solve_system();
    MyPoissonSolver.export_solution_csv("3D_PIN_JUNCTION_SOLUTION.csv");
    MyPoissonSolver.add_solution_to_mesh_functions("ArminxPoissonxSolution");
    fileMSH.export_as_msh("3D_TEST_POISSON_DIODE_PIN.msh", {}, 1);
    const std::string FileName = "TEST_POISSON_DIODE_PN_5V.vtk";

    uepm::file::export_as_vtk(*(p_mesh), FileName, {}, {}, true);

    p_mesh->export_all_vertices_data_to_csv("POISSON_PIN_CSV.csv");

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

// TEST_CASE("Testing Poisson 3d on a FANCY PN Junction.") {
//     static const std::string file_input_test_msh = PROJECT_SRC_DIR + std::string("/example/data/FANCY_3D_PN_DIODE.STF");
//     uepm::file::STF_file           fileMSH(file_input_test_msh);
//     fileMSH.read_mesh();
//     fileMSH.read_states();
//     uepm::mesh::mesh* p_mesh     = fileMSH.get_p_mesh();
//     std::size_t nbVertices = p_mesh->get_nb_vertices();

//     std::cout << "Start building Poisson system ..." << std::endl;
//     uepm::fem::poisson_solver_3d MyPoissonSolver(p_mesh, p_mesh->get_nb_vertices());
//     MyPoissonSolver.compute_stiffness_matrix();
//     std::cout << "Start computing second member ..." << std::endl;
//     MyPoissonSolver.update_second_member();
    
//     MyPoissonSolver.apply_dirichlet_condition("kathode", 1.0);
//     // MyPoissonSolver.apply_dirichlet_condition("kathode", 0.0);
//     // MyPoissonSolver.apply_dirichlet_condition("anode", -0.432793);
//     std::cout << "Start computing poisson solution ..." << std::endl;

//     MyPoissonSolver.decompose_matrix();
//     MyPoissonSolver.solve_system();
//     // MyPoissonSolver.export_solution_csv("3D_PN_JUNCTION_SOLUTION.csv");
//     MyPoissonSolver.add_solution_to_mesh_functions("ArminxSolution");
//     fileMSH.export_as_STF_from_template("FANCY_3D_PN_DIODE_0V.STF",file_input_test_msh);

//     auto   solution     = MyPoissonSolver.get_solution();
//     double min_solution = *std::min_element(solution.begin(), solution.end());
//     double max_solution = *std::max_element(solution.begin(), solution.end());

//     LOG_ERROR << "MIN SOL = " << min_solution;
//     LOG_ERROR << "MAX SOL = " << max_solution;

//     const double min_test_si_ge = -0.765151;
//     const double max_test_si_ge = 5.28987;

    // CHECK(min_solution == doctest::Approx(min_test_si_ge));
    // CHECK_EQ(max_solution, doctest::Approx(max_test_si_ge));
// }

