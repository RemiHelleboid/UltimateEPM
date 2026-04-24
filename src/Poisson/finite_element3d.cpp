/**
 * @file finite_element3d.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-06
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "finite_element3d.hpp"

namespace uepm::fem {

Eigen::Matrix4d FiniteElementP1System3d::compute_elementary_stiffness_matrix(std::shared_ptr<mesh::element> tetra) {
    mesh::vertex* p_vtx1 = tetra->get_vertex(0);
    mesh::vertex* p_vtx2 = tetra->get_vertex(1);
    mesh::vertex* p_vtx3 = tetra->get_vertex(2);
    mesh::vertex* p_vtx4 = tetra->get_vertex(3);

    double X14 = (p_vtx4->x() - p_vtx1->x());
    double X13 = (p_vtx3->x() - p_vtx1->x());
    double X12 = (p_vtx2->x() - p_vtx1->x());

    double Y14 = (p_vtx4->y() - p_vtx1->y());
    double Y13 = (p_vtx3->y() - p_vtx1->y());
    double Y12 = (p_vtx2->y() - p_vtx1->y());

    double Z14 = (p_vtx4->z() - p_vtx1->z());
    double Z13 = (p_vtx3->z() - p_vtx1->z());
    double Z12 = (p_vtx2->z() - p_vtx1->z());

    //  The elementary matrix is symmetric so we just compute 10 coefficients.

    double M11 = pow(X12, 2) * pow(Y13, 2) - 2 * pow(X12, 2) * Y13 * Y14 + pow(X12, 2) * pow(Y14, 2) + pow(X12, 2) * pow(Z13, 2) -
                 2 * pow(X12, 2) * Z13 * Z14 + pow(X12, 2) * pow(Z14, 2) - 2 * X12 * X13 * Y12 * Y13 + 2 * X12 * X13 * Y12 * Y14 +
                 2 * X12 * X13 * Y13 * Y14 - 2 * X12 * X13 * pow(Y14, 2) - 2 * X12 * X13 * Z12 * Z13 + 2 * X12 * X13 * Z12 * Z14 +
                 2 * X12 * X13 * Z13 * Z14 - 2 * X12 * X13 * pow(Z14, 2) + 2 * X12 * X14 * Y12 * Y13 - 2 * X12 * X14 * Y12 * Y14 -
                 2 * X12 * X14 * pow(Y13, 2) + 2 * X12 * X14 * Y13 * Y14 + 2 * X12 * X14 * Z12 * Z13 - 2 * X12 * X14 * Z12 * Z14 -
                 2 * X12 * X14 * pow(Z13, 2) + 2 * X12 * X14 * Z13 * Z14 + pow(X13, 2) * pow(Y12, 2) - 2 * pow(X13, 2) * Y12 * Y14 +
                 pow(X13, 2) * pow(Y14, 2) + pow(X13, 2) * pow(Z12, 2) - 2 * pow(X13, 2) * Z12 * Z14 + pow(X13, 2) * pow(Z14, 2) -
                 2 * X13 * X14 * pow(Y12, 2) + 2 * X13 * X14 * Y12 * Y13 + 2 * X13 * X14 * Y12 * Y14 - 2 * X13 * X14 * Y13 * Y14 -
                 2 * X13 * X14 * pow(Z12, 2) + 2 * X13 * X14 * Z12 * Z13 + 2 * X13 * X14 * Z12 * Z14 - 2 * X13 * X14 * Z13 * Z14 +
                 pow(X14, 2) * pow(Y12, 2) - 2 * pow(X14, 2) * Y12 * Y13 + pow(X14, 2) * pow(Y13, 2) + pow(X14, 2) * pow(Z12, 2) -
                 2 * pow(X14, 2) * Z12 * Z13 + pow(X14, 2) * pow(Z13, 2) + pow(Y12, 2) * pow(Z13, 2) - 2 * pow(Y12, 2) * Z13 * Z14 +
                 pow(Y12, 2) * pow(Z14, 2) - 2 * Y12 * Y13 * Z12 * Z13 + 2 * Y12 * Y13 * Z12 * Z14 + 2 * Y12 * Y13 * Z13 * Z14 -
                 2 * Y12 * Y13 * pow(Z14, 2) + 2 * Y12 * Y14 * Z12 * Z13 - 2 * Y12 * Y14 * Z12 * Z14 - 2 * Y12 * Y14 * pow(Z13, 2) +
                 2 * Y12 * Y14 * Z13 * Z14 + pow(Y13, 2) * pow(Z12, 2) - 2 * pow(Y13, 2) * Z12 * Z14 + pow(Y13, 2) * pow(Z14, 2) -
                 2 * Y13 * Y14 * pow(Z12, 2) + 2 * Y13 * Y14 * Z12 * Z13 + 2 * Y13 * Y14 * Z12 * Z14 - 2 * Y13 * Y14 * Z13 * Z14 +
                 pow(Y14, 2) * pow(Z12, 2) - 2 * pow(Y14, 2) * Z12 * Z13 + pow(Y14, 2) * pow(Z13, 2);

    double M22 = pow(X13, 2) * pow(Y14, 2) + pow(X13, 2) * pow(Z14, 2) - 2 * X13 * X14 * Y13 * Y14 - 2 * X13 * X14 * Z13 * Z14 +
                 pow(X14, 2) * pow(Y13, 2) + pow(X14, 2) * pow(Z13, 2) + pow(Y13, 2) * pow(Z14, 2) - 2 * Y13 * Y14 * Z13 * Z14 +
                 pow(Y14, 2) * pow(Z13, 2);

    double M33 = pow(X12, 2) * pow(Y14, 2) + pow(X12, 2) * pow(Z14, 2) - 2 * X12 * X14 * Y12 * Y14 - 2 * X12 * X14 * Z12 * Z14 +
                 pow(X14, 2) * pow(Y12, 2) + pow(X14, 2) * pow(Z12, 2) + pow(Y12, 2) * pow(Z14, 2) - 2 * Y12 * Y14 * Z12 * Z14 +
                 pow(Y14, 2) * pow(Z12, 2);

    double M44 = pow(X12, 2) * pow(Y13, 2) + pow(X12, 2) * pow(Z13, 2) - 2 * X12 * X13 * Y12 * Y13 - 2 * X12 * X13 * Z12 * Z13 +
                 pow(X13, 2) * pow(Y12, 2) + pow(X13, 2) * pow(Z12, 2) + pow(Y12, 2) * pow(Z13, 2) - 2 * Y12 * Y13 * Z12 * Z13 +
                 pow(Y13, 2) * pow(Z12, 2);

    double M12 = -X12 * X13 * Y13 * Y14 + X12 * X13 * pow(Y14, 2) - X12 * X13 * Z13 * Z14 + X12 * X13 * pow(Z14, 2) +
                 X12 * X14 * pow(Y13, 2) - X12 * X14 * Y13 * Y14 + X12 * X14 * pow(Z13, 2) - X12 * X14 * Z13 * Z14 +
                 pow(X13, 2) * Y12 * Y14 - pow(X13, 2) * pow(Y14, 2) + pow(X13, 2) * Z12 * Z14 - pow(X13, 2) * pow(Z14, 2) -
                 X13 * X14 * Y12 * Y13 - X13 * X14 * Y12 * Y14 + 2 * X13 * X14 * Y13 * Y14 - X13 * X14 * Z12 * Z13 - X13 * X14 * Z12 * Z14 +
                 2 * X13 * X14 * Z13 * Z14 + pow(X14, 2) * Y12 * Y13 - pow(X14, 2) * pow(Y13, 2) + pow(X14, 2) * Z12 * Z13 -
                 pow(X14, 2) * pow(Z13, 2) - Y12 * Y13 * Z13 * Z14 + Y12 * Y13 * pow(Z14, 2) + Y12 * Y14 * pow(Z13, 2) -
                 Y12 * Y14 * Z13 * Z14 + pow(Y13, 2) * Z12 * Z14 - pow(Y13, 2) * pow(Z14, 2) - Y13 * Y14 * Z12 * Z13 -
                 Y13 * Y14 * Z12 * Z14 + 2 * Y13 * Y14 * Z13 * Z14 + pow(Y14, 2) * Z12 * Z13 - pow(Y14, 2) * pow(Z13, 2);

    double M13 = pow(X12, 2) * Y13 * Y14 - pow(X12, 2) * pow(Y14, 2) + pow(X12, 2) * Z13 * Z14 - pow(X12, 2) * pow(Z14, 2) -
                 X12 * X13 * Y12 * Y14 + X12 * X13 * pow(Y14, 2) - X12 * X13 * Z12 * Z14 + X12 * X13 * pow(Z14, 2) - X12 * X14 * Y12 * Y13 +
                 2 * X12 * X14 * Y12 * Y14 - X12 * X14 * Y13 * Y14 - X12 * X14 * Z12 * Z13 + 2 * X12 * X14 * Z12 * Z14 -
                 X12 * X14 * Z13 * Z14 + X13 * X14 * pow(Y12, 2) - X13 * X14 * Y12 * Y14 + X13 * X14 * pow(Z12, 2) - X13 * X14 * Z12 * Z14 -
                 pow(X14, 2) * pow(Y12, 2) + pow(X14, 2) * Y12 * Y13 - pow(X14, 2) * pow(Z12, 2) + pow(X14, 2) * Z12 * Z13 +
                 pow(Y12, 2) * Z13 * Z14 - pow(Y12, 2) * pow(Z14, 2) - Y12 * Y13 * Z12 * Z14 + Y12 * Y13 * pow(Z14, 2) -
                 Y12 * Y14 * Z12 * Z13 + 2 * Y12 * Y14 * Z12 * Z14 - Y12 * Y14 * Z13 * Z14 + Y13 * Y14 * pow(Z12, 2) -
                 Y13 * Y14 * Z12 * Z14 - pow(Y14, 2) * pow(Z12, 2) + pow(Y14, 2) * Z12 * Z13;

    double M14 = -pow(X12, 2) * pow(Y13, 2) + pow(X12, 2) * Y13 * Y14 - pow(X12, 2) * pow(Z13, 2) + pow(X12, 2) * Z13 * Z14 +
                 2 * X12 * X13 * Y12 * Y13 - X12 * X13 * Y12 * Y14 - X12 * X13 * Y13 * Y14 + 2 * X12 * X13 * Z12 * Z13 -
                 X12 * X13 * Z12 * Z14 - X12 * X13 * Z13 * Z14 - X12 * X14 * Y12 * Y13 + X12 * X14 * pow(Y13, 2) - X12 * X14 * Z12 * Z13 +
                 X12 * X14 * pow(Z13, 2) - pow(X13, 2) * pow(Y12, 2) + pow(X13, 2) * Y12 * Y14 - pow(X13, 2) * pow(Z12, 2) +
                 pow(X13, 2) * Z12 * Z14 + X13 * X14 * pow(Y12, 2) - X13 * X14 * Y12 * Y13 + X13 * X14 * pow(Z12, 2) -
                 X13 * X14 * Z12 * Z13 - pow(Y12, 2) * pow(Z13, 2) + pow(Y12, 2) * Z13 * Z14 + 2 * Y12 * Y13 * Z12 * Z13 -
                 Y12 * Y13 * Z12 * Z14 - Y12 * Y13 * Z13 * Z14 - Y12 * Y14 * Z12 * Z13 + Y12 * Y14 * pow(Z13, 2) -
                 pow(Y13, 2) * pow(Z12, 2) + pow(Y13, 2) * Z12 * Z14 + Y13 * Y14 * pow(Z12, 2) - Y13 * Y14 * Z12 * Z13;

    double M23 = -X12 * X13 * pow(Y14, 2) - X12 * X13 * pow(Z14, 2) + X12 * X14 * Y13 * Y14 + X12 * X14 * Z13 * Z14 +
                 X13 * X14 * Y12 * Y14 + X13 * X14 * Z12 * Z14 - pow(X14, 2) * Y12 * Y13 - pow(X14, 2) * Z12 * Z13 -
                 Y12 * Y13 * pow(Z14, 2) + Y12 * Y14 * Z13 * Z14 + Y13 * Y14 * Z12 * Z14 - pow(Y14, 2) * Z12 * Z13;

    double M24 = X12 * X13 * Y13 * Y14 + X12 * X13 * Z13 * Z14 - X12 * X14 * pow(Y13, 2) - X12 * X14 * pow(Z13, 2) -
                 pow(X13, 2) * Y12 * Y14 - pow(X13, 2) * Z12 * Z14 + X13 * X14 * Y12 * Y13 + X13 * X14 * Z12 * Z13 + Y12 * Y13 * Z13 * Z14 -
                 Y12 * Y14 * pow(Z13, 2) - pow(Y13, 2) * Z12 * Z14 + Y13 * Y14 * Z12 * Z13;

    double M34 = -pow(X12, 2) * Y13 * Y14 - pow(X12, 2) * Z13 * Z14 + X12 * X13 * Y12 * Y14 + X12 * X13 * Z12 * Z14 +
                 X12 * X14 * Y12 * Y13 + X12 * X14 * Z12 * Z13 - X13 * X14 * pow(Y12, 2) - X13 * X14 * pow(Z12, 2) -
                 pow(Y12, 2) * Z13 * Z14 + Y12 * Y13 * Z12 * Z14 + Y12 * Y14 * Z12 * Z13 - Y13 * Y14 * pow(Z12, 2);

    const double    multiplication_factor = -1.0 / (36 * fabs(tetra->get_measure()));
    Eigen::Matrix4d ElementaryStiffnessMatrix{{M11, M12, M13, M14}, {M12, M22, M23, M24}, {M13, M23, M33, M34}, {M14, M24, M34, M44}};
    ElementaryStiffnessMatrix *= multiplication_factor;
    return ElementaryStiffnessMatrix;
}

Eigen::Vector4d FiniteElementP1System3d::compute_elementary_second_member(std::shared_ptr<mesh::element> tetra, double constant_value) {
    constexpr double one_fourth = 1.0 / 4.0;
    double           value      = fabs(tetra->get_measure()) * constant_value * one_fourth;
    Eigen::Vector4d  ElementarySecondMember{value, value, value, value};
    return ElementarySecondMember;
}

Eigen::Vector4d FiniteElementP1System3d::compute_elementary_second_member(std::shared_ptr<mesh::element>                tetra,
                                                                          std::function<double(double, double, double)> function) {
    constexpr double      one_fourth = 1.0 / 4.0;
    const mesh::vector3   barycenter = tetra->get_barycenter();
    const double          value      = fabs(tetra->get_measure()) * function(barycenter.x(), barycenter.y(), barycenter.z()) * one_fourth;
    const Eigen::Vector4d ElementarySecondMember{value, value, value, value};
    return ElementarySecondMember;
}

void FiniteElementP1System3d::compute_mass_matrix() {}

void FiniteElementP1System3d::compute_stiffness_matrix() {
    LOG_INFO << "Computing the stiffness matrix of the FEM system 3D.";
    std::size_t                                 number_vertices = m_p_mesh->get_nb_vertices();
    std::vector<std::shared_ptr<mesh::element>> list_p_tetra    = m_p_mesh->get_list_bulk_element();
    m_matrix_lhs.resize(number_vertices, number_vertices);
    constexpr int number_element_per_line_matrix = 16;
    m_matrix_lhs.reserve(Eigen::VectorXd::Constant(number_vertices, number_element_per_line_matrix));

    for (auto&& p_tetra : list_p_tetra) {
        std::vector<mesh::vertex*> p_vertices_list             = p_tetra->get_vertices();
        Eigen::Matrix4d            elementary_stiffness_matrix = compute_elementary_stiffness_matrix(p_tetra);
        for (int index_row = 0; index_row < 4; ++index_row) {
            for (int index_col = 0; index_col < 4; ++index_col) {
                m_matrix_lhs.coeffRef(p_vertices_list[index_row]->get_index(), p_vertices_list[index_col]->get_index()) +=
                    elementary_stiffness_matrix(index_row, index_col);
            }
        }
    }
    m_matrix_lhs.makeCompressed();
}

void FiniteElementP1System3d::compute_second_member(double constant_value) {
    LOG_INFO << "Computing second member of FEM system.";
    std::size_t                                 number_vertices = m_p_mesh->get_nb_vertices();
    std::vector<std::shared_ptr<mesh::element>> list_p_tetra    = m_p_mesh->get_list_bulk_element();
    m_second_member                                             = EigenVector::Constant(number_vertices, 0.0);
    for (auto&& p_tetra : list_p_tetra) {
        std::vector<mesh::vertex*> p_vertices_list        = p_tetra->get_vertices();
        Eigen::Vector4d            ElementarySecondMember = compute_elementary_second_member(p_tetra, constant_value);
        for (int index_row = 0; index_row < 4; ++index_row) {
            m_second_member(p_vertices_list[index_row]->get_index()) += ElementarySecondMember(index_row);
        }
    }
}

void FiniteElementP1System3d::compute_second_member(std::function<double(double, double, double)> function) {
    LOG_INFO << "Computing second member of FEM system.";
    std::size_t                                 number_vertices = m_p_mesh->get_nb_vertices();
    std::vector<std::shared_ptr<mesh::element>> list_p_tetra    = m_p_mesh->get_list_bulk_element();
    m_second_member                                             = EigenVector::Constant(number_vertices, 0.0);
    for (auto&& p_tetra : list_p_tetra) {
        std::vector<mesh::vertex*> p_vertices_list        = p_tetra->get_vertices();
        Eigen::Vector4d            ElementarySecondMember = compute_elementary_second_member(p_tetra, function);
        for (int index_row = 0; index_row < 4; ++index_row) {
            m_second_member(p_vertices_list[index_row]->get_index()) += ElementarySecondMember(index_row);
        }
    }
}

void FiniteElementP1System3d::apply_dirichlet_condition(const std::string& region_name, const double boundary_value) {
    LOG_INFO << "Applying Dirichlet BC on region : " << region_name;
    std::cout << "Applying Dirichlet BC on region : " << region_name << " with value : " << boundary_value << std::endl;
    const mesh::region*      boundary_region              = m_p_mesh->get_p_region(region_name);
    std::vector<std::size_t> region_unique_vertices_index = boundary_region->get_unique_vertices_as_vector();
    for (auto&& index_vertex : region_unique_vertices_index) {
        m_matrix_lhs.coeffRef(index_vertex, index_vertex) += FiniteElementSystem::very_large_value;
        m_second_member(index_vertex) += FiniteElementSystem::very_large_value * boundary_value;
    }
}

void FiniteElementP1System3d::apply_dirichlet_condition_second_member(const std::string& region_name, const double boundary_value) {
    const mesh::region*      boundary_region              = m_p_mesh->get_p_region(region_name);
    std::vector<std::size_t> region_unique_vertices_index = boundary_region->get_unique_vertices_as_vector();
    for (auto&& index_vertex : region_unique_vertices_index) {
        m_second_member(index_vertex) = FiniteElementSystem::very_large_value * boundary_value;
    }
}

void FiniteElementP1System3d::apply_neuman_condition(const std::string& region_name, const double boundary_value) {
    LOG_INFO << "Applying Neuman BC on region : " << region_name;
    constexpr double    one_third       = 1.0 / 3.0;
    const mesh::region* boundary_region = m_p_mesh->get_p_region(region_name);
    auto                list_elements   = boundary_region->get_list_elements();
    for (auto&& sp_element : list_elements) {
        double triangle_surface = sp_element->get_measure();
        auto   list_vtx         = sp_element->get_vertices();
        for (auto&& p_vtx : list_vtx) {
            m_second_member(p_vtx->get_index()) += one_third * triangle_surface * boundary_value;
        }
    }
}

}  // namespace uepm::fem
