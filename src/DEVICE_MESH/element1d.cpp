/**
 * @file element1d.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-13
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "element1d.hpp"

#include <random>

namespace uepm {

namespace mesh {

std::vector<std::array<std::size_t, 2>> element1d::get_edges_as_index_pair() const {
    const std::vector<std::size_t>   vect_vtx_idx = get_vertices_index();
    const std::array<std::size_t, 2> edge_1       = {vect_vtx_idx[0], vect_vtx_idx[1]};
    return {edge_1};
}

std::vector<double> element1d::compute_barycentric_coordinate(const vector3 &location) const {
    const double lambda_0 = distance(*m_vertices[0], location) / get_length();
    return {lambda_0, 1 - lambda_0};
}

vector3 element1d::compute_surface_normal() const {
    const vector3 vtx_0 = *m_vertices[0];
    const vector3 vtx_1 = *m_vertices[1];
    const vector3 edge_1{vtx_1 - vtx_0};
    vector3       normal_vector{-edge_1.y(), edge_1.x(), 0.0};
    normal_vector *= 1.0 / normal_vector.norm();
    return normal_vector;
}

bool element1d::is_location_inside_element(const vector3 &location) const {
    const std::vector<double> barycentric_coords = compute_barycentric_coordinate(location);
    return (barycentric_coords[0] >= 0 && barycentric_coords[1] >= 0);
}

vector3 element1d::compute_gradient(const std::string &scalar_field_name) const {
    const double micron_to_cm = 1e-4;
    const double value_0      = m_vertices[0]->get_scalar_data(scalar_field_name);
    const double value_1      = m_vertices[1]->get_scalar_data(scalar_field_name);
    const double x_0          = m_vertices[0]->x() * micron_to_cm;
    const double x_1          = m_vertices[1]->x() * micron_to_cm;
    const double distance     = x_1 - x_0;
    const double grad_x       = (value_1 - value_0) / (distance * micron_to_cm);
    return {grad_x, 0.0, 0.0};
}

double element1d::integrate_scalar(const std::string &scalar_field_name) const {
    const double micron_to_cm = 1e-4;
    const double value_0      = m_vertices[0]->get_scalar_data(scalar_field_name);
    const double value_1      = m_vertices[1]->get_scalar_data(scalar_field_name);
    const double x_0          = m_vertices[0]->x() * micron_to_cm;
    const double x_1          = m_vertices[1]->x() * micron_to_cm;
    const double distance     = x_1 - x_0;
    return (value_0 + value_1) * distance * micron_to_cm / 2.0;
}

vector3 element1d::integrate_vector(const std::string &vector_field_name) const {
    const double micron_to_cm = 1e-4;
    const double value_0_x    = m_vertices[0]->get_vector_data(vector_field_name).x();
    const double value_1_x    = m_vertices[1]->get_vector_data(vector_field_name).x();
    const double x_0          = m_vertices[0]->x() * micron_to_cm;
    const double x_1          = m_vertices[1]->x() * micron_to_cm;
    const double distance     = x_1 - x_0;
    return {(value_0_x + value_1_x) * distance * micron_to_cm / 2.0, 0.0, 0.0};
}

std::map<std::shared_ptr<element>, vector3> element1d::compute_element_line_intersection(const vector3 &point_A,
                                                                                         const vector3 &point_B) const {
    std::map<std::shared_ptr<element>, vector3> map_position_face_intersections;
    std::optional<vector3>                      opt_intersection_point =
        compute_line_line_intersection(point_A, point_B, *m_vertices[0], *m_vertices[1]);
    if (opt_intersection_point.has_value()) {
        map_position_face_intersections[std::make_shared<element1d>(*this)] = opt_intersection_point.value();
    }
    return map_position_face_intersections;
}

/**
 * @brief Return a uniformly drawn point on the 1d element.
 *
 * @return vector3
 */
vector3 element1d::draw_uniform_random_point_inside_element() const {
    std::mt19937                           random_generator{std::random_device()()};
    std::uniform_real_distribution<double> uniform_unit_distribution{std::uniform_real_distribution<double>{0.0, 1.0}};
    const double                           random_distance_to_vtxA = uniform_unit_distribution(random_generator);
    return random_distance_to_vtxA * *m_vertices[0] + (1 - random_distance_to_vtxA) * *m_vertices[1];
}

vector3 element1d::draw_uniform_random_point_inside_element(std::minstd_rand &random_generator) const {
    std::uniform_real_distribution<double> uniform_unit_distribution{std::uniform_real_distribution<double>{0.0, 1.0}};
    const double                           random_distance_to_vtxA = uniform_unit_distribution(random_generator);
    return random_distance_to_vtxA * *m_vertices[0] + (1 - random_distance_to_vtxA) * *m_vertices[1];
}


}  // namespace mesh

}  // namespace uepm