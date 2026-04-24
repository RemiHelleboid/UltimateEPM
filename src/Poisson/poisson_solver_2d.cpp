/**
 * @file poisson_solver_2d.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-06
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "poisson_solver_2d.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

// #include "config.h"
#include "finite_element.hpp"
#include "finite_element2d.hpp"

#ifdef USE_OPENMP_ACCELERATION
#include <omp.h>
constexpr int USE_OPEN_MP_POISSON = 1;
#endif

namespace uepm::fem {

void poisson_solver_2d::compute_stiffness_matrix() {
    LOG_INFO << "Computing the stiffness matrix of the Poisson 2D system.";
    std::size_t number_vertices = m_p_mesh->get_nb_vertices();
    m_matrix_lhs.resize(number_vertices, number_vertices);
    constexpr int number_element_per_line_matrix = 12;
    m_matrix_lhs.reserve(Eigen::VectorXd::Constant(number_vertices, number_element_per_line_matrix));

    for (const auto &bulk_region : m_p_mesh->get_all_p_bulk_region()) {
        std::string                material_name   = bulk_region->get_material();
        physic::material::material region_material = m_list_materials.get_material(material_name);
        if (region_material.m_name == "Unknown") {
            LOG_ERROR << "UNKNOWN MATERIAL IN POISSON SOLVER : " << material_name;
            throw std::invalid_argument("The material given in Poisson solver is unknown : " + material_name);
        }
        double relative_permittivity = region_material.m_parameters["dielectric-constant"];
        auto   list_sp_elements      = bulk_region->get_list_elements();
        for (auto &&sp_element : list_sp_elements) {
            m_list_bulk_elements.push_back(sp_element);
            std::vector<mesh::vertex *> p_vertices_list       = sp_element->get_vertices();
            double                      absolute_permittivity = relative_permittivity * physic::constant::vacuum_permittivity;
            Eigen::Matrix3d             elementary_stiffness_matrix =
                (absolute_permittivity / physic::constant::elementary_charge) * compute_elementary_stiffness_matrix(sp_element);
            for (int index_row = 0; index_row < 3; ++index_row) {
                for (int index_col = 0; index_col < 3; ++index_col) {
                    m_matrix_lhs.coeffRef(p_vertices_list[index_row]->get_index(), p_vertices_list[index_col]->get_index()) +=
                        elementary_stiffness_matrix(index_row, index_col);
                }
            }
        }
    }
    m_matrix_lhs.makeCompressed();
}

Eigen::Vector3d poisson_solver_2d::compute_charge_density_elementary_second_member(const std::shared_ptr<mesh::element> &sp_triangle) {
    const auto       p_vtx     = sp_triangle->get_vertices();
    constexpr double cm3_to_m3 = 1e-6;

    const double    value_p1 = p_vtx[0]->get_space_charge() * cm3_to_m3;
    const double    value_p2 = p_vtx[1]->get_space_charge() * cm3_to_m3;
    const double    value_p3 = p_vtx[2]->get_space_charge() * cm3_to_m3;
    Eigen::Vector3d ElementarySecondMember{value_p1, value_p2, value_p3};

    constexpr double one_third = 1.0 / 3.0;
    ElementarySecondMember *= sp_triangle->get_measure() * one_third;

    return ElementarySecondMember;
}

void poisson_solver_2d::update_second_member() {
    m_second_member = EigenVector::Constant(get_system_size(), 0.0);
    for (std::size_t index_triangle = 0; index_triangle < m_list_bulk_elements.size(); ++index_triangle) {
        const std::vector<mesh::vertex *> p_vertices_list = m_list_bulk_elements[index_triangle]->get_vertices();
        Eigen::Vector3d ElementarySecondMember = -compute_charge_density_elementary_second_member(m_list_bulk_elements[index_triangle]);
        for (int index_row = 0; index_row < 3; ++index_row) {
            m_second_member(p_vertices_list[index_row]->get_index()) += ElementarySecondMember(index_row);
        }
    }
}

}  // namespace fem
