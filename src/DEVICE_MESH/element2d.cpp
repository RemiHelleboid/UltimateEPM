/**
 * @file element2d.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2021-09-08
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "element2d.hpp"

#include <array>
#include <cassert>
#include <fstream>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace uepm {

namespace mesh {

inline std::vector<std::array<std::size_t, 2>> element2d::get_edges_as_index_pair() const {
    const auto &indices = get_vertices_index();
    return {{indices[0], indices[1]}, {indices[1], indices[2]}, {indices[2], indices[0]}};
}

std::vector<std::pair<vertex *, vertex *>> element2d::get_edges_as_vertex_pair() const {
    return {{m_vertices[0], m_vertices[1]}, {m_vertices[1], m_vertices[2]}, {m_vertices[2], m_vertices[0]}};
}

std::vector<element1d> element2d::get_list_edges() const {
    return {element1d(m_vertices[0], m_vertices[1]), element1d(m_vertices[1], m_vertices[2]), element1d(m_vertices[2], m_vertices[0])};
}

double element2d::get_perimeters() const {
    const double l1 = distance_between(m_vertices[0], m_vertices[1]);
    const double l2 = distance_between(m_vertices[1], m_vertices[2]);
    const double l3 = distance_between(m_vertices[2], m_vertices[0]);
    return l1 + l2 + l3;
}

double element2d::compute_area() const {
    const vector3 v0 = *m_vertices[0];
    const vector3 v1 = *m_vertices[1];
    const vector3 v2 = *m_vertices[2];

    const double l1 = (v1 - v0).norm();
    const double l2 = (v2 - v1).norm();
    const double l3 = (v0 - v2).norm();

    const double s = 0.5 * (l1 + l2 + l3);

    return std::sqrt(s * (s - l1) * (s - l2) * (s - l3));
}

vector3 element2d::compute_surface_normal() const {
    const vector3 vtx_0 = *m_vertices[0];
    const vector3 vtx_1 = *m_vertices[1];
    const vector3 vtx_2 = *m_vertices[2];
    const vector3 edge_1{vtx_1 - vtx_0};
    const vector3 edge_2{vtx_2 - vtx_0};
    vector3       normal_vector = cross_product(edge_1, edge_2);
    normal_vector *= 1.0 / normal_vector.norm();
    return normal_vector;
}

void element2d::compute_precomputed_field_for_barycentric_coordinates() {
    m_bary_coord_precomputed_v00 = m_precomputed_v0.dot(m_precomputed_v0);
    m_bary_coord_precomputed_v01 = m_precomputed_v0.dot(m_precomputed_v1);
    m_bary_coord_precomputed_v11 = m_precomputed_v1.dot(m_precomputed_v1);
    m_inverse_bary_coord_precomputed_divisor = 1.0 / (
        m_bary_coord_precomputed_v00 * m_bary_coord_precomputed_v11 - m_bary_coord_precomputed_v01 * m_bary_coord_precomputed_v01);
}

std::vector<double> element2d::compute_barycentric_coordinate(const vector3 &location) const {
    const vector3 v2       = location - *m_vertices[0];
    const double  value_20 = v2.dot(m_precomputed_v0);
    const double  value_21 = v2.dot(m_precomputed_v1);

    const double lambda_2 =
        (m_bary_coord_precomputed_v11 * value_20 - m_bary_coord_precomputed_v01 * value_21) * m_inverse_bary_coord_precomputed_divisor;
    const double lambda_3 =
        (m_bary_coord_precomputed_v00 * value_21 - m_bary_coord_precomputed_v01 * value_20) * m_inverse_bary_coord_precomputed_divisor;
    const double lambda_1 = 1 - lambda_2 - lambda_3;

    return {lambda_1, lambda_2, lambda_3};
}

bool element2d::is_location_inside_element(const vector3 &location) const {
    const std::vector<double> barycentric_coords = compute_barycentric_coordinate(location);
    return std::none_of(barycentric_coords.begin(), barycentric_coords.end(), [](const double &value) { return value < 0.0; });
}

vector3 element2d::compute_gradient(const std::string &scalar_field_name) const {
    const double micron_to_cm = 1.0;
    const double x_0          = m_vertices[0]->x();
    const double y_0          = m_vertices[0]->y();
    const double x_1          = m_vertices[1]->x();
    const double y_1          = m_vertices[1]->y();
    const double x_2          = m_vertices[2]->x();
    const double y_2          = m_vertices[2]->y();
    const double value_0      = m_vertices[0]->get_scalar_data(scalar_field_name);
    const double value_1      = m_vertices[1]->get_scalar_data(scalar_field_name);
    const double value_2      = m_vertices[2]->get_scalar_data(scalar_field_name);

    const double surface = x_0 * y_1 - x_0 * y_2 - x_1 * y_0 + x_1 * y_2 + x_2 * y_0 - x_2 * y_1;
    if (surface == 0.0) {
        std::cerr << "Error : surface is null" << std::endl;
        return {0.0, 0.0, 0.0};
    }

    const double grad_x =
        (value_0 * y_1 - value_0 * y_2 - value_1 * y_0 + value_1 * y_2 + value_2 * y_0 - value_2 * y_1) / (surface * micron_to_cm);
    const double grad_y =
        (-value_0 * x_1 + value_0 * x_2 + value_1 * x_0 - value_1 * x_2 - value_2 * x_0 + value_2 * x_1) / (surface * micron_to_cm);

    return {grad_x, grad_y, 0.0};
}

double element2d::integrate_scalar(const std::string &scalar_field_name) const {
    const double micron2_to_cm2 = 1e-8;
    const double x_0            = m_vertices[0]->x();
    const double y_0            = m_vertices[0]->y();
    const double x_1            = m_vertices[1]->x();
    const double y_1            = m_vertices[1]->y();
    const double x_2            = m_vertices[2]->x();
    const double y_2            = m_vertices[2]->y();
    const double value_0        = m_vertices[0]->get_scalar_data(scalar_field_name);
    const double value_1        = m_vertices[1]->get_scalar_data(scalar_field_name);
    const double value_2        = m_vertices[2]->get_scalar_data(scalar_field_name);

    const double surface = x_0 * y_1 - x_0 * y_2 - x_1 * y_0 + x_1 * y_2 + x_2 * y_0 - x_2 * y_1;
    if (surface == 0.0) {
        std::cerr << "Error : surface is null" << std::endl;
        return 0.0;
    }

    double integral = (1.0 / 3.0) * (value_0 + value_1 + value_2) * get_measure() * micron2_to_cm2;

    return integral;
}

vector3 element2d::integrate_vector(const std::string &vector_field_name) const {
    const double  micron2_to_cm2 = 1e-8;
    const double  x_0            = m_vertices[0]->x();
    const double  y_0            = m_vertices[0]->y();
    const double  x_1            = m_vertices[1]->x();
    const double  y_1            = m_vertices[1]->y();
    const double  x_2            = m_vertices[2]->x();
    const double  y_2            = m_vertices[2]->y();
    const vector3 value_0        = m_vertices[0]->get_vector_data(vector_field_name);
    const vector3 value_1        = m_vertices[1]->get_vector_data(vector_field_name);
    const vector3 value_2        = m_vertices[2]->get_vector_data(vector_field_name);

    const double surface = x_0 * y_1 - x_0 * y_2 - x_1 * y_0 + x_1 * y_2 + x_2 * y_0 - x_2 * y_1;
    if (surface == 0.0) {
        std::cerr << "Error : surface is null" << std::endl;
        return {0.0, 0.0, 0.0};
    }

    const double int_x = (value_0.x() + value_1.x() + value_2.x()) * surface * micron2_to_cm2 / 3.0;
    const double int_y = (value_0.y() + value_1.y() + value_2.y()) * surface * micron2_to_cm2 / 3.0;
    const double int_z = (value_0.z() + value_1.z() + value_2.z()) * surface * micron2_to_cm2 / 3.0;

    return {int_x, int_y, int_z};
}

/**
 * @brief
 *
 * @param point_A
 * @param point_B
 * @return std::map<std::shared_ptr<element>, vector3>
 */
std::map<std::shared_ptr<element>, vector3> element2d::compute_element_line_intersection(const vector3 &point_A,
                                                                                         const vector3 &point_B) const {
    std::map<std::shared_ptr<element>, vector3> map_position_face_intersections;
    const auto                                  list_edges_as_vtx_pair = get_edges_as_vertex_pair();
    for (const auto &pair_vtx : list_edges_as_vtx_pair) {
        auto intersection_result = compute_line_line_intersection(point_A, point_B, *(pair_vtx.first), *(pair_vtx.second));
        if (intersection_result.has_value()) {
            mesh::element1d face{pair_vtx.first, pair_vtx.second};
            map_position_face_intersections[std::make_shared<element1d>(face)] = intersection_result.value();
        }
    }
    return map_position_face_intersections;
}

/**
 * @brief This function use the Möller–Trumbore intersection algorithm to compute the potential intersection
 * beteween a segment and a triangle in 3D.
 * If there is an intersection, the function return the intersection point.
 * If there is not intersection, it returns an empty std::optional object.
 *
 * @param point_A
 * @param point_B
 * @return std::optional<vector3>
 */
std::optional<vector3> element2d::compute_line_triangle_intersection_3d(const vector3 &point_A, const vector3 &point_B) const {
    const double  epsilon_intersection = 1e-12;
    const vector3 vtx_0                = *m_vertices[0];
    const vector3 vtx_1                = *m_vertices[1];
    const vector3 vtx_2                = *m_vertices[2];
    const vector3 edge_1{vtx_1 - vtx_0};
    const vector3 edge_2{vtx_2 - vtx_0};
    const vector3 ray_segment{point_B - point_A};
    const vector3 cross_vector_1 = cross_product(ray_segment, edge_2);
    const double  product_1      = edge_1.dot(cross_vector_1);
    // Case 1 ; product_1 == 0 -> ray and triangle are parallel.
    if (product_1 < epsilon_intersection && product_1 > -epsilon_intersection) {
        return {};
    }
    const double  inverse_product_1 = 1.0 / product_1;
    const vector3 triangle_ray_origin_edge{point_A - vtx_0};
    const double  product_2 = inverse_product_1 * triangle_ray_origin_edge.dot(cross_vector_1);
    if (product_2 < 0.0 || product_2 > 1.0) {
        return {};
    }
    const vector3 cross_vector_2 = cross_product(triangle_ray_origin_edge, edge_1);
    const double  product_3      = inverse_product_1 * ray_segment.dot(cross_vector_2);
    if (product_3 < 0.0 || product_2 + product_3 > 1.0) {
        return {};
    }
    const double intersection_barycentric_coordinate = inverse_product_1 * edge_2.dot(cross_vector_2);
    if (intersection_barycentric_coordinate >= epsilon_intersection) {
        const vector3 point_intersection =
            point_A * (1 - intersection_barycentric_coordinate) + intersection_barycentric_coordinate * point_B;
        if (is_point_between_two_others(point_A, point_B, point_intersection)) {
            return point_intersection;
        }
        return std::nullopt;
    }
    return std::nullopt;
}

/**
 * @brief Draw a random point inside the triangle.
 *
 * It is done by using the barycentric coordinates
 *
 * @return vector3
 */
vector3 element2d::draw_uniform_random_point_inside_element() const {
    std::mt19937                           random_generator{std::random_device()()};
    std::uniform_real_distribution<double> uniform_unit_distribution{std::uniform_real_distribution<double>{0.0, 1.0}};

    // Generate two random numbers for barycentric coordinates
    double r1 = uniform_unit_distribution(random_generator);
    double r2 = uniform_unit_distribution(random_generator);

    // Transform to ensure uniform distribution inside the triangle
    double sqrt_r1            = std::sqrt(r1);
    double lambda_random_vtxA = 1 - sqrt_r1;
    double lambda_random_vtxB = sqrt_r1 * (1 - r2);
    double lambda_vtxC        = sqrt_r1 * r2;

    vector3 random_point = lambda_random_vtxA * *m_vertices[0] + lambda_random_vtxB * *m_vertices[1] + lambda_vtxC * *m_vertices[2];

    // Check if the point is inside the triangle.
    if (!is_location_inside_element(random_point)) {
        std::cerr << "Error : Random point is not inside the triangle." << std::endl;
        throw std::invalid_argument("Random point is not inside the triangle.");
    }
    return random_point;
}

vector3 element2d::draw_uniform_random_point_inside_element(std::minstd_rand &random_generator) const {
    std::uniform_real_distribution<double> uniform_unit_distribution{std::uniform_real_distribution<double>{0.0, 1.0}};

    // Generate two random numbers for barycentric coordinates
    double r1 = uniform_unit_distribution(random_generator);
    double r2 = uniform_unit_distribution(random_generator);

    // Transform to ensure uniform distribution inside the triangle
    double sqrt_r1            = std::sqrt(r1);
    double lambda_random_vtxA = 1 - sqrt_r1;
    double lambda_random_vtxB = sqrt_r1 * (1 - r2);
    double lambda_vtxC        = sqrt_r1 * r2;

    vector3 random_point = lambda_random_vtxA * *m_vertices[0] + lambda_random_vtxB * *m_vertices[1] + lambda_vtxC * *m_vertices[2];

    // Check if the point is inside the triangle.
    if (!is_location_inside_element(random_point)) {
        std::cerr << "Error : Random point is not inside the triangle." << std::endl;
        throw std::invalid_argument("Random point is not inside the triangle.");
    }
    return random_point;
}


}  // namespace mesh

}  // namespace uepm