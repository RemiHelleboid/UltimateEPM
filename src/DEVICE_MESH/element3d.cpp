/**
 * @file element3d.cpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief
 * @version 0.1
 * @date 2021-09-02
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "element3d.hpp"

#include <algorithm>
#include <memory>
#include <random>
#include <set>

#include "bbox.hpp"
#include "element1d.hpp"
#include "element2d.hpp"

namespace uepm {

namespace mesh {

std::vector<std::array<std::size_t, 2>> element3d::get_edges_as_index_pair() const {
    const std::vector<std::size_t>   vector_vtx_idx = get_vertices_index();
    const std::array<std::size_t, 2> edge_1         = {vector_vtx_idx[0], vector_vtx_idx[1]};
    const std::array<std::size_t, 2> edge_2         = {vector_vtx_idx[1], vector_vtx_idx[2]};
    const std::array<std::size_t, 2> edge_3         = {vector_vtx_idx[2], vector_vtx_idx[0]};
    const std::array<std::size_t, 2> edge_4         = {vector_vtx_idx[3], vector_vtx_idx[0]};
    const std::array<std::size_t, 2> edge_5         = {vector_vtx_idx[3], vector_vtx_idx[1]};
    const std::array<std::size_t, 2> edge_6         = {vector_vtx_idx[3], vector_vtx_idx[2]};
    return {edge_1, edge_2, edge_3, edge_4, edge_5, edge_6};
}

std::vector<element1d> element3d::get_list_edges() const {
    return {element1d(m_vertices[0], m_vertices[1]), element1d(m_vertices[1], m_vertices[2]), element1d(m_vertices[2], m_vertices[0]),
            element1d(m_vertices[3], m_vertices[0]), element1d(m_vertices[3], m_vertices[1]), element1d(m_vertices[3], m_vertices[2])};
}

std::vector<element2d> element3d::get_faces_as_element2d() const {
    return {element2d(m_vertices[0], m_vertices[1], m_vertices[3]), element2d(m_vertices[0], m_vertices[3], m_vertices[2]),
            element2d(m_vertices[1], m_vertices[2], m_vertices[3]), element2d(m_vertices[2], m_vertices[1], m_vertices[0])};
}

vector3 element3d::compute_surface_normal() const {
    // Return an exception.
    throw std::invalid_argument("Normal surface is not available for tetrahedrons.");
}

std::vector<double> element3d::compute_barycentric_coordinate(const vector3 &location) const {
    const vector3 v_loc1        = location - *m_vertices[0];
    const double  signed_volume = 6.0 * m_signed_volume;
    const double  lambda_2      = scalar_triple_product(v_loc1, m_v1, m_v2) / signed_volume;
    const double  lambda_3      = scalar_triple_product(v_loc1, m_v2, m_v0) / signed_volume;
    const double  lambda_4      = scalar_triple_product(v_loc1, m_v0, m_v1) / signed_volume;
    const double  lambda_1      = 1.0 - lambda_2 - lambda_3 - lambda_4;
    return {lambda_1, lambda_2, lambda_3, lambda_4};
}

bool element3d::is_location_inside_element(const vector3 &location) const {
    const vector3 v_loc1        = location - *m_vertices[0];
    const vector3 v_loc2        = location - *m_vertices[1];
    const double  signed_volume = 6.0 * m_signed_volume;
    const double  lambda_2      = scalar_triple_product(v_loc1, m_v1, m_v2) / signed_volume;
    const double  lambda_3      = scalar_triple_product(v_loc1, m_v2, m_v0) / signed_volume;
    const double  lambda_4      = scalar_triple_product(v_loc1, m_v0, m_v1) / signed_volume;
    const double  lambda_1      = scalar_triple_product(v_loc2, m_v4, m_v3) / signed_volume;
    return (lambda_1 >= 0 && lambda_2 >= 0 && lambda_3 >= 0 && lambda_4 >= 0);
}

vector3 element3d::compute_gradient(const std::string &scalar_field_name) const {
    // Positions
    const double x0 = m_vertices[0]->x(), y0 = m_vertices[0]->y(), z0 = m_vertices[0]->z();
    const double x1 = m_vertices[1]->x(), y1 = m_vertices[1]->y(), z1 = m_vertices[1]->z();
    const double x2 = m_vertices[2]->x(), y2 = m_vertices[2]->y(), z2 = m_vertices[2]->z();
    const double x3 = m_vertices[3]->x(), y3 = m_vertices[3]->y(), z3 = m_vertices[3]->z();

    // Edges from vertex 0
    const vector3 a{x1 - x0, y1 - y0, z1 - z0};  // v1 - v0
    const vector3 b{x2 - x0, y2 - y0, z2 - z0};  // v2 - v0
    const vector3 c{x3 - x0, y3 - y0, z3 - z0};  // v3 - v0

    // Scalar values
    const double u0 = m_vertices[0]->get_scalar_data(scalar_field_name);
    const double u1 = m_vertices[1]->get_scalar_data(scalar_field_name);
    const double u2 = m_vertices[2]->get_scalar_data(scalar_field_name);
    const double u3 = m_vertices[3]->get_scalar_data(scalar_field_name);

    // Differences relative to vertex 0
    const double du1 = u1 - u0;
    const double du2 = u2 - u0;
    const double du3 = u3 - u0;

    // det(R) = a · (b × c)  (6 * signed volume)
    const double     det = a.dot(cross_product(b, c));
    constexpr double eps = 1e-14;
    if (std::abs(det) < eps) {
        // Degenerate or nearly degenerate tetrahedron – choose a policy:
        // return zero, assert, or throw. Here we return zero.
        return vector3{0.0, 0.0, 0.0};
    }

    // // R^{-T} Δu => gradient
    const vector3 grad = (cross_product(b, c) * du1 + cross_product(c, a) * du2 + cross_product(a, b) * du3) / det;
    return grad;
}

double element3d::integrate_scalar(const std::string &scalar_field_name) const {
    constexpr double CM_TO_MICRON = 1e-12;
    double           volume       = this->get_volume() * CM_TO_MICRON;
    double           integral     = m_vertices[0]->get_scalar_data(scalar_field_name) + m_vertices[1]->get_scalar_data(scalar_field_name) +
                      m_vertices[2]->get_scalar_data(scalar_field_name) + m_vertices[3]->get_scalar_data(scalar_field_name);
    constexpr double NUM_VERTICES = 4.0;
    return integral * volume / NUM_VERTICES;
}

vector3 element3d::integrate_vector(const std::string &vector_field_name) const {
    double  volume   = this->get_volume();
    vector3 integral = m_vertices[0]->get_vector_data(vector_field_name) + m_vertices[1]->get_vector_data(vector_field_name) +
                       m_vertices[2]->get_vector_data(vector_field_name) + m_vertices[3]->get_vector_data(vector_field_name);
    constexpr double NUM_VERTICES = 4.0;
    return integral * volume / NUM_VERTICES;
}

std::map<std::shared_ptr<element>, vector3> element3d::compute_element_line_intersection(const vector3 &point_A,
                                                                                         const vector3 &point_B) const {
    std::map<std::shared_ptr<element>, vector3> map_position_face_intersections;
    auto                                        list_faces = get_faces_as_element2d();
    for (const auto &face : list_faces) {
        auto intersection_result = face.compute_line_triangle_intersection_3d(point_A, point_B);
        if (intersection_result.has_value()) {
            std::shared_ptr<element2d> face_intersection       = std::make_shared<element2d>(face);
            map_position_face_intersections[face_intersection] = intersection_result.value();
        }
    }
    return map_position_face_intersections;
}

/**
 * @brief Draw a random point inside the tetrahedron.
 *
 * It is done by using the barycentic coordinates.
 * TO BE CHECKED, not sure the distribution is uniform.
 *
 * @return vector3
 */
vector3 element3d::draw_uniform_random_point_inside_element() const {
    std::mt19937                           random_generator{std::random_device()()};
    std::uniform_real_distribution<double> uniform_unit_distribution{std::uniform_real_distribution<double>{0.0, 1.0}};

    double lambda_random_vtxA = uniform_unit_distribution(random_generator);
    double lambda_random_vtxB = uniform_unit_distribution(random_generator);
    double lambda_random_vtxC = uniform_unit_distribution(random_generator);
    double lambda_vtxD        = uniform_unit_distribution(random_generator);
    double lambda_sum         = lambda_random_vtxA + lambda_random_vtxB + lambda_random_vtxC + lambda_vtxD;
    lambda_random_vtxA /= lambda_sum;
    lambda_random_vtxB /= lambda_sum;
    lambda_random_vtxC /= lambda_sum;
    lambda_vtxD /= lambda_sum;

    vector3 random_point = lambda_random_vtxA * *m_vertices[0] + lambda_random_vtxB * *m_vertices[1] + lambda_random_vtxC * *m_vertices[2] +
                           lambda_vtxD * *m_vertices[3];

    return random_point;
}

vector3 element3d::draw_uniform_random_point_inside_element(std::minstd_rand &random_generator) const {
    std::uniform_real_distribution<double> uniform_unit_distribution{std::uniform_real_distribution<double>{0.0, 1.0}};

    double lambda_random_vtxA = uniform_unit_distribution(random_generator);
    double lambda_random_vtxB = uniform_unit_distribution(random_generator);
    double lambda_random_vtxC = uniform_unit_distribution(random_generator);
    double lambda_vtxD        = uniform_unit_distribution(random_generator);
    double lambda_sum         = lambda_random_vtxA + lambda_random_vtxB + lambda_random_vtxC + lambda_vtxD;
    lambda_random_vtxA /= lambda_sum;
    lambda_random_vtxB /= lambda_sum;
    lambda_random_vtxC /= lambda_sum;
    lambda_vtxD /= lambda_sum;

    vector3 random_point = lambda_random_vtxA * *m_vertices[0] + lambda_random_vtxB * *m_vertices[1] + lambda_random_vtxC * *m_vertices[2] +
                           lambda_vtxD * *m_vertices[3];

    return random_point;
}


}  // namespace mesh

}  // namespace uepm