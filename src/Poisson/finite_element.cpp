/**
 * @file poisson2d.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-09-30
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "finite_element.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <fstream>

namespace uepm::fem {

std::string get_status_string(fem_status status) {
    if (status == fem_status::None) {
        return "None.";
    }
    if (status == fem_status::Assembled) {
        return "Assembled.";
    }
    if (status == fem_status::SolvedFailed) {
        return "Solver failed.";
    }
    if (status == fem_status::SolvedSuccess) {
        return "Solver success.";
    }
    return "Unknown";
}

void FiniteElementSystem::reset_system() {
    m_matrix_lhs    = EigenSparseMatrix(0, 0);
    m_solution      = EigenVector(0);
    m_second_member = EigenVector(0);
    m_status        = fem_status::None;
}

const double FiniteElementSystem::very_large_value = 1e120;

void FiniteElementSystem::decompose_matrix() {
    Eigen::SparseLU<EigenSparseMatrix, Eigen::COLAMDOrdering<int> > solver;
    m_solver.analyzePattern(m_matrix_lhs);
    m_solver.factorize(m_matrix_lhs);
    // m_iter_solver.setTolerance(1e-6);
    // m_iter_solver.compute(m_matrix_lhs);
    if (m_solver.info() != Eigen::Success) {
        LOG_ERROR << "The matrix decomposition failed.";
        throw std::runtime_error("The matrix decomposition failed.");
    }
}

void FiniteElementSystem::solve_system() {
    // auto start = std::chrono::high_resolution_clock::now();
    m_solution = m_solver.solve(m_second_member);
    if (m_solver.info() != Eigen::Success) {
        LOG_ERROR << "The resolution of the linear system failed. ";
        throw std::runtime_error("The resolution of the linear system failed. ");
    }
    // auto stop     = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    // print_time_different_units(duration.count());
}

double FiniteElementSystem::solver_L2_error() const {
    double error_l2 = (m_matrix_lhs * m_solution - m_second_member).norm();  //.norm() is an Eigen method that compute l2 norm of an array
    return (error_l2);
}

void FiniteElementSystem::export_solution_csv(const std::string &filename) {
    std::ofstream file_csv(filename);
    file_csv << "X,Y,Z,Solution" << std::endl;
    for (const auto &vtx : m_p_mesh->get_list_vertices()) {
        file_csv << vtx.x() << ',';
        file_csv << vtx.y() << ',';
        file_csv << vtx.z() << ',';
        file_csv << m_solution(vtx.get_index()) << "\n";
    }
    file_csv.close();
}

void FiniteElementSystem::add_solution_to_mesh_functions(const std::string &function_name, bool add_gradient) {
    std::vector<double> solution_values;
    solution_values.resize(m_solution.size());
    EigenVector::Map(&solution_values[0], m_solution.size()) = m_solution;
    m_p_mesh->create_scalar_function_from_values_on_vertex(function_name, solution_values);
    if (add_gradient) {
        m_p_mesh->create_gradient_function(function_name, function_name + "_gradient");
        const double gradient_scale = 1.0;
        m_p_mesh->add_electric_field_to_vertices(function_name + "_gradient", gradient_scale);
    }
}

void   FiniteElementSystem::add_second_member_to_mesh_functions(const std::string& function_name) {
    std::vector<double> second_member_values;
    second_member_values.resize(m_second_member.size());
    EigenVector::Map(&second_member_values[0], m_second_member.size()) = m_second_member;
    m_p_mesh->create_scalar_function_from_values_on_vertex(function_name, second_member_values);    
}


void FiniteElementSystem::print_infos() {
    std::cout << "FINITE ELEMENT SYSTEM INFOS" << std::endl;
    std::cout << "System size                       : " << m_solution.size() << std::endl;
    std::cout << "System status                     : " << get_status_string(m_status) << std::endl;

    std::cout << "Solution minimum                  : " << m_solution.minCoeff() << std::endl;
    std::cout << "Solution maximum                  : " << m_solution.maxCoeff() << std::endl;
    std::cout << "L2 Error of Solver                : " << solver_L2_error() << std::endl;
}

}  // namespace fem