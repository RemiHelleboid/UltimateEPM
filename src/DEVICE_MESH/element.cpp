/**
 * @file element.cpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief
 * @version 0.1
 * @date 2021-07-10
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "element.hpp"

#include <algorithm>
#include <cassert>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "bbox.hpp"

namespace uepm {

namespace mesh {

vertex *element::get_vertex(const std::size_t index) const {
    if (index >= m_vertices.size()) {
        throw std::invalid_argument("Try to access to a vertex with an index superior to the number of vertices in the element.");
    }
    return m_vertices[index];
}

std::vector<std::size_t> element::get_vertices_index() const {
    std::vector<std::size_t> index_vertices(m_vertices.size());
    std::transform(m_vertices.begin(), m_vertices.end(), index_vertices.begin(), [](const vertex * v) { return v->get_index(); });
    return (index_vertices);
}

bool element::contains_vertex_with_index(std::size_t index_vtx) const {
    return std::find_if(m_vertices.begin(), m_vertices.end(), [index_vtx](const vertex *vtx) { return vtx->get_index() == index_vtx; }) !=
           m_vertices.end();
}

/**
 * @brief Two tetrahedron are adjacent if they have 3 vertices in common.
 *
 * @param second_element
 * @return true
 * @return false
 */
std::optional<std::vector<std::size_t>> element::is_adjacent_to_element(const element &second_element) const {
    const std::vector<std::size_t> list_index_element_1 = get_vertices_index();
    const std::vector<std::size_t> list_index_element_2 = second_element.get_vertices_index();
    const std::set<std::size_t>    set_vertices_1(list_index_element_1.begin(), list_index_element_1.end());
    const std::set<std::size_t>    set_vertices_2(list_index_element_2.begin(), list_index_element_2.end());
    std::vector<std::size_t>       vector_common_vtx;
    std::set_intersection(set_vertices_1.begin(),
                          set_vertices_1.end(),
                          set_vertices_2.begin(),
                          set_vertices_2.end(),
                          std::back_inserter(vector_common_vtx));
    if (vector_common_vtx.empty()) {
        return std::nullopt;
    }
    return vector_common_vtx;
}

unsigned int element::get_element_type() const {
    constexpr unsigned int element_type_edge        = 1;
    constexpr unsigned int element_type_triangle    = 2;
    constexpr unsigned int element_type_tetrahedron = 4;
    if (m_vertices.size() == 2) {
        return element_type_edge;
    }
    if (m_vertices.size() == 3) {
        return element_type_triangle;
    }
    if (m_vertices.size() == 4) {
        return element_type_tetrahedron;
    }
    assert("Error with element type : number of vertices not in {2, 3, 4}.");
    return -1;
}

vector3 element::get_barycenter() const {
    vector3 barycenter{0.0, 0.0, 0.0};
    for (auto &&p_vtx : m_vertices) {
        barycenter += *p_vtx;
    }
    barycenter *= 1.0 / m_vertices.size();
    return barycenter;
}

bbox element::compute_bounding_box() const {
    std::vector<double> X_coords(m_vertices.size());
    std::vector<double> Y_coords(m_vertices.size());
    std::vector<double> Z_coords(m_vertices.size());
    std::transform(m_vertices.begin(), m_vertices.end(), X_coords.begin(), [&](const auto &p_vtx) { return p_vtx->x(); });
    std::transform(m_vertices.begin(), m_vertices.end(), Y_coords.begin(), [&](const auto &p_vtx) { return p_vtx->y(); });
    std::transform(m_vertices.begin(), m_vertices.end(), Z_coords.begin(), [&](const auto &p_vtx) { return p_vtx->z(); });
    const double x_min = *std::min_element(X_coords.begin(), X_coords.end());
    const double x_max = *std::max_element(X_coords.begin(), X_coords.end());
    const double y_min = *std::min_element(Y_coords.begin(), Y_coords.end());
    const double y_max = *std::max_element(Y_coords.begin(), Y_coords.end());
    const double z_min = *std::min_element(Z_coords.begin(), Z_coords.end());
    const double z_max = *std::max_element(Z_coords.begin(), Z_coords.end());
    return bbox(x_min, x_max, y_min, y_max, z_min, z_max);
}

double element::interpolate_scalar_at_location(const std::string &name, const vector3 &location) const {
    const auto barycentric_coords = compute_barycentric_coordinate(location);
    double     interpolated_value = 0;
    for (unsigned int index_vtx = 0; index_vtx < m_vertices.size(); ++index_vtx) {
        interpolated_value += barycentric_coords[index_vtx] * m_vertices[index_vtx]->get_scalar_data(name);
    }
    return interpolated_value;
}

vector3 element::interpolate_electric_field_at_location(const vector3 &location) const {
    const auto barycentric_coords = compute_barycentric_coordinate(location);
    vector3    interpolated_value(0.0, 0.0, 0.0);
    for (unsigned int index_vtx = 0; index_vtx < m_vertices.size(); ++index_vtx) {
        interpolated_value += barycentric_coords[index_vtx] * m_vertices[index_vtx]->get_electric_field();
    }
    return interpolated_value;
}

double element::interpolate_doping_at_location(const vector3 &location) const {
    const auto barycentric_coords = compute_barycentric_coordinate(location);
    double     interpolated_value = 0;
    for (unsigned int index_vtx = 0; index_vtx < m_vertices.size(); ++index_vtx) {
        interpolated_value += barycentric_coords[index_vtx] * m_vertices[index_vtx]->get_doping_concentration();
    }
    return interpolated_value;
}

vector3 element::interpolate_e_grad_diffusion_at_location(const vector3 &location) const {
    const auto barycentric_coords = compute_barycentric_coordinate(location);
    vector3    interpolated_value(0.0, 0.0, 0.0);
    for (unsigned int index_vtx = 0; index_vtx < m_vertices.size(); ++index_vtx) {
        interpolated_value += barycentric_coords[index_vtx] * m_vertices[index_vtx]->get_e_grad_diffusion();
    }
    return interpolated_value;
}

vector3 element::interpolate_h_grad_diffusion_at_location(const vector3 &location) const {
    const auto barycentric_coords = compute_barycentric_coordinate(location);
    vector3    interpolated_value(0.0, 0.0, 0.0);
    for (unsigned int index_vtx = 0; index_vtx < m_vertices.size(); ++index_vtx) {
        interpolated_value += barycentric_coords[index_vtx] * m_vertices[index_vtx]->get_h_grad_diffusion();
    }
    return interpolated_value;
}

vector3 element::interpolate_vector_at_location(const std::string &name, const vector3 &location) const {
    const auto barycentric_coords = compute_barycentric_coordinate(location);
    vector3    interpolated_value(0.0, 0.0, 0.0);
    for (unsigned int index_vtx = 0; index_vtx < m_vertices.size(); ++index_vtx) {
        interpolated_value += barycentric_coords[index_vtx] * m_vertices[index_vtx]->get_vector_data(name);
    }
    return interpolated_value;
}

void element::distribute_charge_quantity_over_vertices(double charge_quantity) {
    constexpr double micro_m3_to_cm3 = 1.0e-12;
    const double     volume_element  = fabs(this->get_measure());
    for (auto &&p_vtx : m_vertices) {
        p_vtx->add_charge_to_vertex(charge_quantity / (volume_element * micro_m3_to_cm3));
    }
}

}  // namespace mesh

}  // namespace uepm