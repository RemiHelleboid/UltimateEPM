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

#include <fmt/core.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <cmath>
#include <fstream>
#include <stdexcept>

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
    m_solution      = EigenVector::Zero(0);
    m_second_member = EigenVector::Zero(0);
    m_status        = fem_status::None;
}

const double FiniteElementSystem::very_large_value = 1e120;

void FiniteElementSystem::decompose_matrix() {
    Eigen::SparseLU<EigenSparseMatrix, Eigen::COLAMDOrdering<int> > solver;
    m_solver.analyzePattern(m_matrix_lhs);
    m_solver.factorize(m_matrix_lhs);
    if (m_solver.info() != Eigen::Success) {
        fmt::print("The matrix decomposition failed.\n");
        throw std::runtime_error("The matrix decomposition failed.");
    }
}

void FiniteElementSystem::solve_system() {
    m_solution = m_solver.solve(m_second_member);
    if (m_solver.info() != Eigen::Success) {
        fmt::print("The resolution of the linear system failed.\n");
        throw std::runtime_error("The resolution of the linear system failed. ");
    }
    if (!m_solution.allFinite()) {
        throw std::runtime_error("The Poisson linear system produced a non-finite solution.");
    }
}

void FiniteElementSystem::set_solution(const EigenVector &solution) {
    if (solution.size() != m_solution.size()) {
        throw std::invalid_argument("FiniteElementSystem::set_solution size mismatch.");
    }
    if (!solution.allFinite()) {
        throw std::invalid_argument("FiniteElementSystem::set_solution received a non-finite solution.");
    }
    m_solution = solution;
}

void FiniteElementSystem::mix_solution_with(const EigenVector &old_solution, double old_solution_fraction) {
    if (old_solution.size() != m_solution.size()) {
        throw std::invalid_argument("FiniteElementSystem::mix_solution_with size mismatch.");
    }
    if (!std::isfinite(old_solution_fraction) || old_solution_fraction < 0.0 || old_solution_fraction > 1.0) {
        throw std::invalid_argument("Poisson mixing old_solution_fraction must be in [0, 1].");
    }
    if (!old_solution.allFinite()) {
        throw std::invalid_argument("Poisson mixing old solution contains non-finite values.");
    }
    if (!m_solution.allFinite()) {
        throw std::invalid_argument("Poisson mixing new solution contains non-finite values.");
    }
    m_solution = old_solution_fraction * old_solution + (1.0 - old_solution_fraction) * m_solution;
    if (!m_solution.allFinite()) {
        throw std::runtime_error("Poisson mixing produced a non-finite solution.");
    }
}

double FiniteElementSystem::solver_L2_error() const {
    double error_l2 = (m_matrix_lhs * m_solution - m_second_member)
                          .norm();  //.norm() is an Eigen method that compute l2 norm of an array
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
        constexpr double gradient_scale = 1.0;
        m_p_mesh->add_electric_field_to_vertices(function_name + "_gradient", gradient_scale);
    }
}

void FiniteElementSystem::update_mesh_electric_field_from_solution() {
    const mesh::vector3 null_vector{0.0, 0.0, 0.0};
    std::vector<mesh::vector3> vector_values(m_p_mesh->get_nb_vertices(), null_vector);
    std::vector<double>        vector_average_renormalization(m_p_mesh->get_nb_vertices(), 0.0);

    m_p_mesh->for_each_bulk_element([&](const mesh::element &element) {
        const auto &vertices = element.get_vertices();
        if (vertices.size() != 3) {
            return;
        }

        const double x_0     = vertices[0]->x();
        const double y_0     = vertices[0]->y();
        const double x_1     = vertices[1]->x();
        const double y_1     = vertices[1]->y();
        const double x_2     = vertices[2]->x();
        const double y_2     = vertices[2]->y();
        const double value_0 = m_solution(vertices[0]->get_index());
        const double value_1 = m_solution(vertices[1]->get_index());
        const double value_2 = m_solution(vertices[2]->get_index());

        const double surface = x_0 * y_1 - x_0 * y_2 - x_1 * y_0 + x_1 * y_2 + x_2 * y_0 - x_2 * y_1;
        if (surface == 0.0) {
            return;
        }

        const double grad_x =
            (value_0 * y_1 - value_0 * y_2 - value_1 * y_0 + value_1 * y_2 + value_2 * y_0 - value_2 * y_1) /
            surface;
        const double grad_y =
            (-value_0 * x_1 + value_0 * x_2 + value_1 * x_0 - value_1 * x_2 - value_2 * x_0 + value_2 * x_1) /
            surface;
        mesh::vector3 gradient_at_element{grad_x, grad_y, 0.0};
        if (std::isnan(gradient_at_element.norm())) {
            gradient_at_element = null_vector;
        }

        for (const auto *vertex : vertices) {
            const auto index = vertex->get_index();
            vector_values[index] += gradient_at_element;
            vector_average_renormalization[index] += 1.0;
        }
    });

    constexpr double electric_field_scale = -1.0e4;
    for (std::size_t index_vtx = 0; index_vtx < m_p_mesh->get_nb_vertices(); ++index_vtx) {
        auto *vertex = m_p_mesh->get_p_vertex(index_vtx);
        if (vertex == nullptr || vector_average_renormalization[index_vtx] == 0.0) {
            continue;
        }
        vector_values[index_vtx] *= electric_field_scale / vector_average_renormalization[index_vtx];
        vertex->set_electric_field(vector_values[index_vtx]);
    }
}

void FiniteElementSystem::add_second_member_to_mesh_functions(const std::string &function_name) {
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

}  // namespace uepm::fem
