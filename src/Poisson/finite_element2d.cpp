/**
 * @file finite_element2d.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-02
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "finite_element2d.hpp"

namespace uepm::fem {

const Eigen::Matrix3d FiniteElementP1System2d::ElementaryMassMatrixRefElement{{2.0, 1.0, 1.0}, {1.0, 2.0, 1.0}, {1.0, 1.0, 2.0}};

Eigen::Matrix3d FiniteElementP1System2d::compute_elementary_stiffness_matrix(std::shared_ptr<mesh::element> triangle) {
    mesh::vertex* p_vtx1 = triangle->get_vertex(0);
    mesh::vertex* p_vtx2 = triangle->get_vertex(1);
    mesh::vertex* p_vtx3 = triangle->get_vertex(2);
    double        X13    = p_vtx3->x() - p_vtx1->x();
    double        X23    = p_vtx3->x() - p_vtx2->x();
    double        X12    = p_vtx2->x() - p_vtx1->x();
    double        Y13    = p_vtx3->y() - p_vtx1->y();
    double        Y23    = p_vtx3->y() - p_vtx2->y();
    double        Y12    = p_vtx2->y() - p_vtx1->y();

    //  The elementary matrix is symmetric so we just compute 6 coefficients.

    double M11 = Y23 * Y23 + X23 * X23;
    double M22 = Y13 * Y13 + X13 * X13;
    double M33 = Y12 * Y12 + X12 * X12;
    double M12 = (Y23 * (-Y13)) + (X23 * (-X13));
    double M13 = ((Y12)*Y23) + ((X12)*X23);
    double M23 = ((-Y13) * (Y12)) + ((-X13) * (X12));

    const double multiplication_factor = -1.0 / (4 * fabs(triangle->get_measure()));

    Eigen::Matrix3d ElementaryStiffnessMatrix{{M11 * multiplication_factor, M12 * multiplication_factor, M13 * multiplication_factor},
                                              {M12 * multiplication_factor, M22 * multiplication_factor, M23 * multiplication_factor},
                                              {M13 * multiplication_factor, M23 * multiplication_factor, M33 * multiplication_factor}};
    // ElementaryStiffnessMatrix *= multiplication_factor;
    return ElementaryStiffnessMatrix;
}

Eigen::Matrix3d FiniteElementP1System2d::compute_elementary_mass_matrix(std::shared_ptr<mesh::element> triangle) {
    const double    multiplication_factor = (triangle->get_measure() / 12.0);
    Eigen::Matrix3d ElementaryMassMatrix  = multiplication_factor * ElementaryMassMatrixRefElement;
    return ElementaryMassMatrix;
}

Eigen::Vector3d FiniteElementP1System2d::compute_elementary_second_member(std::shared_ptr<mesh::element> triangle, double constant_value) {
    constexpr double one_third = 1.0 / 3.0;
    double           value     = fabs(triangle->get_measure()) * constant_value * one_third;
    Eigen::Vector3d  ElementarySecondMember{value, value, value};
    return ElementarySecondMember;
}

Eigen::Vector3d FiniteElementP1System2d::compute_elementary_second_member(std::shared_ptr<mesh::element>        triangle,
                                                                          std::function<double(double, double)> function) {
    mesh::vector3    barycenter = triangle->get_barycenter();
    constexpr double one_third  = 1.0 / 3.0;
    double           value      = fabs(triangle->get_measure()) * function(barycenter.x(), barycenter.y()) * one_third;
    Eigen::Vector3d  ElementarySecondMember{value, value, value};
    return ElementarySecondMember;
}

void FiniteElementP1System2d::compute_stiffness_matrix() {
    LOG_INFO << "Computing the stiffness matrix of the FEM system.";
    std::size_t                                 number_vertices  = m_p_mesh->get_nb_vertices();
    std::vector<std::shared_ptr<mesh::element>> list_p_triangles = m_p_mesh->get_list_bulk_element();
    m_matrix_lhs.resize(number_vertices, number_vertices);
    constexpr int number_element_per_line_matrix = 12;
    m_matrix_lhs.reserve(Eigen::VectorXd::Constant(number_vertices, number_element_per_line_matrix));

    for (auto&& p_triangle : list_p_triangles) {
        std::vector<mesh::vertex*> p_vertices_list             = p_triangle->get_vertices();
        Eigen::Matrix3d            elementary_stiffness_matrix = compute_elementary_stiffness_matrix(p_triangle);

        for (int index_row = 0; index_row < 3; ++index_row) {
            for (int index_col = 0; index_col < 3; ++index_col) {
                m_matrix_lhs.coeffRef(p_vertices_list[index_row]->get_index(), p_vertices_list[index_col]->get_index()) +=
                    elementary_stiffness_matrix(index_row, index_col);
            }
        }
    }
    m_matrix_lhs.makeCompressed();
}

void FiniteElementP1System2d::compute_mass_matrix() {}

void FiniteElementP1System2d::compute_second_member(double constant_value) {
    LOG_INFO << "Computing second member of FEM system.";
    std::size_t                                 number_vertices  = m_p_mesh->get_nb_vertices();
    std::vector<std::shared_ptr<mesh::element>> list_p_triangles = m_p_mesh->get_list_bulk_element();
    m_second_member                                              = EigenVector::Constant(number_vertices, 0.0);
    for (auto&& p_triangle : list_p_triangles) {
        std::vector<mesh::vertex*> p_vertices_list        = p_triangle->get_vertices();
        Eigen::Vector3d            ElementarySecondMember = compute_elementary_second_member(p_triangle, constant_value);
        for (int index_row = 0; index_row < 3; ++index_row) {
            m_second_member(p_vertices_list[index_row]->get_index()) += ElementarySecondMember(index_row);
        }
    }
}

void FiniteElementP1System2d::compute_second_member(std::function<double(double, double)> function) {
    LOG_INFO << "Computing second member of FEM system.";
    std::size_t                                 number_vertices  = m_p_mesh->get_nb_vertices();
    std::vector<std::shared_ptr<mesh::element>> list_p_triangles = m_p_mesh->get_list_bulk_element();
    m_second_member                                              = EigenVector::Constant(number_vertices, 0.0);
    for (auto&& p_triangle : list_p_triangles) {
        std::vector<mesh::vertex*> p_vertices_list        = p_triangle->get_vertices();
        Eigen::Vector3d            ElementarySecondMember = compute_elementary_second_member(p_triangle, function);
        for (int index_row = 0; index_row < 3; ++index_row) {
            m_second_member(p_vertices_list[index_row]->get_index()) += ElementarySecondMember(index_row);
        }
    }
}

void FiniteElementP1System2d::apply_dirichlet_condition(const std::string& region_name, const double boundary_value) {
    LOG_INFO << "Applying Dirichlet BC on region : " << region_name;
    const mesh::region*      boundary_region              = m_p_mesh->get_p_region(region_name);
    std::vector<std::size_t> region_unique_vertices_index = boundary_region->get_unique_vertices_as_vector();
    for (auto&& index_vertex : region_unique_vertices_index) {
        m_matrix_lhs.coeffRef(index_vertex, index_vertex) += FiniteElementSystem::very_large_value;
        m_second_member(index_vertex) += FiniteElementSystem::very_large_value * boundary_value;
    }
}

void FiniteElementP1System2d::apply_dirichlet_condition_second_member(const std::string& region_name, const double boundary_value) {
    const mesh::region*      boundary_region              = m_p_mesh->get_p_region(region_name);
    std::vector<std::size_t> region_unique_vertices_index = boundary_region->get_unique_vertices_as_vector();
    for (auto&& index_vertex : region_unique_vertices_index) {
        m_second_member(index_vertex) = FiniteElementSystem::very_large_value * boundary_value;
    }
}

void FiniteElementP1System2d::apply_neuman_condition(const std::string& region_name, const double boundary_value) {
    LOG_INFO << "Applying Neuman BC on region : " << region_name;
    constexpr double    one_half        = 1.0 / 2.0;
    const mesh::region* boundary_region = m_p_mesh->get_p_region(region_name);
    auto                list_elements   = boundary_region->get_list_elements();
    for (auto&& sp_element : list_elements) {
        double edge_length   = sp_element->get_measure();
        auto   list_vtx_pair = sp_element->get_edges_as_index_pair();
        for (auto&& vtx_pair : list_vtx_pair) {
            m_second_member(vtx_pair[0]) += one_half * edge_length * boundary_value;
            m_second_member(vtx_pair[1]) += one_half * edge_length * boundary_value;
        }
    }
}

}  // namespace uepm::fem