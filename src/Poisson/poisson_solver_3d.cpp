            /**
 * @file poisson_solver_3d.cpp
 * @author Rémi Helleboid (remi.helleboid@st.com)
 * @brief
 * @version 0.1
 * @date 2021-10-19
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "poisson_solver_3d.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

// #include "config.h"
#include "finite_element.hpp"
#include "finite_element3d.hpp"

#ifdef USE_OPENMP_ACCELERATION
#include <omp.h>
constexpr int USE_OPEN_MP_POISSON = 1;
#endif

namespace uepm::fem {

void poisson_solver_3d::compute_stiffness_matrix() {
    LOG_INFO << "Computing the stiffness matrix of the POISSON system 3D.";
    constexpr int number_element_per_line_matrix = 60;
    std::size_t   number_vertices                = m_p_mesh->get_nb_vertices();
    m_matrix_lhs.resize(number_vertices, number_vertices);
    m_matrix_lhs.reserve(Eigen::VectorXd::Constant(number_vertices, number_element_per_line_matrix));

    int index_element = 0;
    for (const auto &bulk_region : m_p_mesh->get_all_p_bulk_region()) {
        std::string                material_name   = bulk_region->get_material();
        physic::material::material region_material = m_list_materials.get_material(material_name);
        if (region_material.m_name == "Unknown") {
            LOG_ERROR << "UNKNOWN MATERIAL IN POISSON SOLVER : " << material_name;
            throw std::invalid_argument("The material given in Poisson solver is unknown : " + material_name);
        }
        if (region_material.m_name == "Gas") {
            continue;
        }
        double relative_permittivity = region_material.m_parameters["dielectric-constant"];
        auto   list_sp_elements      = bulk_region->get_list_elements();
        double absolute_permittivity = relative_permittivity * physic::constant::vacuum_permittivity;
        for (auto &&sp_element : list_sp_elements) {
            ++index_element;
            m_list_bulk_elements.push_back(sp_element);
            std::vector<mesh::vertex *> p_vertices_list = sp_element->get_vertices();
            Eigen::Matrix4d             elementary_stiffness_matrix =
                (absolute_permittivity / physic::constant::elementary_charge) * compute_elementary_stiffness_matrix(sp_element);
            for (int index_row = 0; index_row < 4; ++index_row) {
                for (int index_col = 0; index_col < 4; ++index_col) {
                    m_matrix_lhs.coeffRef(p_vertices_list[index_row]->get_index(), p_vertices_list[index_col]->get_index()) +=
                        elementary_stiffness_matrix(index_row, index_col);
                }
            }
        }
    }
    m_matrix_lhs.makeCompressed();
}

Eigen::Vector4d poisson_solver_3d::compute_charge_density_elementary_second_member(const std::shared_ptr<mesh::element> &sp_tetra) {
    const auto       p_vtx       = sp_tetra->get_vertices();
    constexpr double cmm3_to_mm3 = 1e-6;
    const double     value_p1    = p_vtx[0]->get_space_charge() * cmm3_to_mm3;
    const double     value_p2    = p_vtx[1]->get_space_charge() * cmm3_to_mm3;
    const double     value_p3    = p_vtx[2]->get_space_charge() * cmm3_to_mm3;
    const double     value_p4    = p_vtx[3]->get_space_charge() * cmm3_to_mm3;
    Eigen::Vector4d  ElementarySecondMember{value_p1, value_p2, value_p3, value_p4};

    constexpr double one_fourth = 1.0 / 4.0;
    ElementarySecondMember *= sp_tetra->get_measure() * one_fourth;

    return ElementarySecondMember;
}

void poisson_solver_3d::update_second_member() {
    m_second_member = EigenVector::Constant(get_system_size(), 0.0);
    for (std::size_t index_tetra = 0; index_tetra < m_list_bulk_elements.size(); ++index_tetra) {
        const std::vector<mesh::vertex *> p_vertices_list = m_list_bulk_elements[index_tetra]->get_vertices();
        Eigen::Vector4d ElementarySecondMember = -compute_charge_density_elementary_second_member(m_list_bulk_elements[index_tetra]);
        for (int index_row = 0; index_row < 4; ++index_row) {
            m_second_member(p_vertices_list[index_row]->get_index()) += ElementarySecondMember(index_row);
        }
    }
}

}  // namespace uepm::fem
