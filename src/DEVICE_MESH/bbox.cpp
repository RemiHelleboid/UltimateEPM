/**
 * @file bbox.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-29
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "bbox.hpp"

#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "element.hpp"

namespace uepm {

namespace mesh {

std::vector<bbox> bbox::split_2d_box_in_quadrants() const {
    const vector3 box_center = get_center();
    const bbox    bottom_left_box(m_x_min, box_center.x(), m_y_min, box_center.y(), 0.0, 0.0);
    const bbox    bottom_right_box(box_center.x(), m_x_max, m_y_min, box_center.y(), 0.0, 0.0);
    const bbox    top_right_box(box_center.x(), m_x_max, box_center.y(), m_y_max, 0.0, 0.0);
    const bbox    top_left_box(m_x_min, box_center.x(), box_center.y(), m_y_max, 0.0, 0.0);
    return {bottom_left_box, bottom_right_box, top_right_box, top_left_box};
}

std::vector<bbox> bbox::split_3d_box_in_octants() const {
    const vector3 box_center = get_center();
    const bbox    bottom_left_front_box(m_x_min, box_center.x(), m_y_min, box_center.y(), m_z_min, box_center.z());
    const bbox    bottom_right_front_box(box_center.x(), m_x_max, m_y_min, box_center.y(), m_z_min, box_center.z());
    const bbox    top_right_front_box(box_center.x(), m_x_max, box_center.y(), m_y_max, m_z_min, box_center.z());
    const bbox    top_left_front_box(m_x_min, box_center.x(), box_center.y(), m_y_max, m_z_min, box_center.z());
    const bbox    bottom_left_back_box(m_x_min, box_center.x(), m_y_min, box_center.y(), box_center.z(), m_z_max);
    const bbox    bottom_right_back_box(box_center.x(), m_x_max, m_y_min, box_center.y(), box_center.z(), m_z_max);
    const bbox    top_right_back_box(box_center.x(), m_x_max, box_center.y(), m_y_max, box_center.z(), m_z_max);
    const bbox    top_left_back_box(m_x_min, box_center.x(), box_center.y(), m_y_max, box_center.z(), m_z_max);

    return {bottom_left_front_box,
            bottom_right_front_box,
            top_right_front_box,
            top_left_front_box,
            bottom_left_back_box,
            bottom_right_back_box,
            top_right_back_box,
            top_left_back_box};
}

bool bbox::is_overlapping(const bbox& second_box) const {
    const bool noOverlap = this->m_x_min > second_box.m_x_max || second_box.m_x_min > this->m_x_max || this->m_y_min > second_box.m_y_max ||
                           second_box.m_y_min > this->m_y_max || this->m_z_min > second_box.m_z_max || second_box.m_z_min > this->m_z_max;
    return !noOverlap;
}

bool bbox::is_overlapping_2d(const bbox& second_box) const {
    const bool noOverlap = this->m_x_min > second_box.m_x_max || second_box.m_x_min > this->m_x_max || this->m_y_min > second_box.m_y_max ||
                           second_box.m_y_min > this->m_y_max;
    return !noOverlap;
}

bool bbox::is_box_overlapping(const element& my_element) const {
    const bool im_overlapping = this->is_overlapping(my_element.get_bounding_box());
    return im_overlapping;
}

bool bbox::is_box_overlapping_2d(const element& my_element) const {
    const bool im_overlapping = this->is_overlapping_2d(my_element.get_bounding_box());
    return im_overlapping;
}

bool bbox::is_inside(const vector3& location) const {
    return (location.x() <= m_x_max && location.x() >= m_x_min && location.y() <= m_y_max && location.y() >= m_y_min &&
            location.z() <= m_z_max && location.z() >= m_z_min);
}

bool bbox::is_inside_2d(const vector3& location) const {
    return (location.x() <= m_x_max && location.x() >= m_x_min && location.y() <= m_y_max && location.y() >= m_y_min);
}

bool bbox::is_overlapping_triangle(const element& triangle) const {
    const auto& list_p_vertices = triangle.get_vertices();
    // If one vertices lies into the box, the the triangle overlapps the box.
    if (is_inside(*list_p_vertices[0]) || is_inside(*list_p_vertices[1]) || is_inside(*list_p_vertices[2])) {
        return true;
    }
    return (this->is_box_overlapping_2d(triangle));
}

bool bbox::is_overlapping_tetra(const element& tetra) const {
    const auto& list_p_vertices = tetra.get_vertices();
    // If one vertices lies into the box, the the tetra overlapps the box.
    if (is_inside(*list_p_vertices[0]) || is_inside(*list_p_vertices[1]) || is_inside(*list_p_vertices[2]) ||
        is_inside(*list_p_vertices[3])) {
        return true;
    }
    return (this->is_overlapping(tetra.get_bounding_box()));
}

bool bbox::is_overlapping_circle_2d(const vector3& center, double radius) const {
    const double x_min          = m_x_min - radius;
    const double x_max          = m_x_max + radius;
    const double y_min          = m_y_min - radius;
    const double y_max          = m_y_max + radius;
    const double x_center       = center.x();
    const double y_center       = center.y();
    const bool   im_overlapping = (x_center - x_min) * (x_center - x_max) <= 0 && (y_center - y_min) * (y_center - y_max) <= 0;
    return im_overlapping;
}

bool bbox::is_overlapping_sphere_3d(const vector3& center, double radius) const {
    const double x_min          = m_x_min - radius;
    const double x_max          = m_x_max + radius;
    const double y_min          = m_y_min - radius;
    const double y_max          = m_y_max + radius;
    const double z_min          = m_z_min - radius;
    const double z_max          = m_z_max + radius;
    const double x_center       = center.x();
    const double y_center       = center.y();
    const double z_center       = center.z();
    const bool   im_overlapping = (x_center - x_min) * (x_center - x_max) <= 0 && (y_center - y_min) * (y_center - y_max) <= 0 &&
                                  (z_center - z_min) * (z_center - z_max) <= 0;
    return im_overlapping;
}

std::ostream& operator<<(std::ostream& os, const bbox& my_box) {
    os << my_box.get_x_min() << "," << my_box.get_x_max() << "," << my_box.get_y_min() << "," << my_box.get_y_max() << ","
       << my_box.get_z_min() << "," << my_box.get_z_max();
    return os;
}

std::vector<vector3> bbox::generate_mesh_grid_2d(std::size_t N_x, std::size_t N_y) const {
    std::vector<vector3> MeshGrid;
    std::cout << "Number of point on the grid : " << N_x * N_y << "\n";
    MeshGrid.reserve(N_x * N_y);
    const double dx = N_x == 0 ? 0 : (m_x_max - m_x_min) / static_cast<double>(N_x - 1);
    const double dy = N_y == 0 ? 0 : (m_y_max - m_y_min) / static_cast<double>(N_y - 1);
    for (std::size_t index_x = 0; index_x < N_x; ++index_x) {
        double x_coord = m_x_min + index_x * dx;
        for (std::size_t index_y = 0; index_y < N_y; ++index_y) {
            double y_coord = m_y_min + index_y * dy;
            MeshGrid.push_back(vector3(x_coord, y_coord, 0.0));
        }
    }
    return MeshGrid;
}

std::vector<vector3> bbox::generate_mesh_grid_3d(std::size_t N_x, std::size_t N_y, std::size_t N_z) const {
    std::vector<vector3> MeshGrid;
    MeshGrid.reserve(N_x * N_y * N_z);
    const double dx = N_x <= 1 ? 0 : (m_x_max - m_x_min) / static_cast<double>(N_x - 1);
    const double dy = N_y <= 1 ? 0 : (m_y_max - m_y_min) / static_cast<double>(N_y - 1);
    const double dz = N_z <= 1 ? 0 : (m_z_max - m_z_min) / static_cast<double>(N_z - 1);
    for (std::size_t index_x = 0; index_x < N_x; ++index_x) {
        double x_coord = m_x_min + dx * index_x;
        for (std::size_t index_y = 0; index_y < N_y; ++index_y) {
            double y_coord = m_y_min + dy * index_y;
            for (std::size_t index_z = 0; index_z < N_z; ++index_z) {
                double z_coord = m_z_min + dz * index_z;
                MeshGrid.push_back(vector3(x_coord, y_coord, z_coord));
            }
        }
    }
    return MeshGrid;
}

std::vector<vector3> bbox::generate_inner_mesh_grid_3d(std::size_t N_x, std::size_t N_y, std::size_t N_z) const {
    constexpr double     one_half = 1.0 / 2.0;
    std::vector<vector3> MeshGrid;
    MeshGrid.reserve(N_x * N_y * N_z);
    const double inner_x = (m_x_max - m_x_min) * (1 - 1.0 / static_cast<double>(N_x - 1));
    const double inner_y = (m_y_max - m_y_min) * (1 - 1.0 / static_cast<double>(N_y - 1));
    const double inner_z = (m_z_max - m_z_min) * (1 - 1.0 / static_cast<double>(N_z - 1));
    const double dx      = inner_x / static_cast<double>(N_x - 1);
    const double dy      = inner_y / static_cast<double>(N_y - 1);
    const double dz      = inner_z / static_cast<double>(N_z - 1);
    for (std::size_t index_x = 0; index_x < N_x; ++index_x) {
        double x_coord = m_x_min + one_half * dx + dx * index_x;
        for (std::size_t index_y = 0; index_y < N_y; ++index_y) {
            double y_coord = m_y_min + one_half * dy + dy * index_y;
            for (std::size_t index_z = 0; index_z < N_z; ++index_z) {
                double z_coord = m_z_min + one_half * dz + dz * index_z;
                MeshGrid.push_back(vector3(x_coord, y_coord, z_coord));
            }
        }
    }
    return MeshGrid;
}

std::vector<vector3> bbox::generate_inner_mesh_grid_2d(std::size_t N_x, std::size_t N_y) const {
    constexpr double     one_half = 1.0 / 2.0;
    std::vector<vector3> MeshGrid;
    MeshGrid.reserve(N_x * N_y);
    const double inner_x = (m_x_max - m_x_min) * (1 - 1.0 / static_cast<double>(N_x - 1));
    const double inner_y = (m_y_max - m_y_min) * (1 - 1.0 / static_cast<double>(N_y - 1));
    const double dx      = inner_x / static_cast<double>(N_x - 1);
    const double dy      = inner_y / static_cast<double>(N_y - 1);
    for (std::size_t index_x = 0; index_x < N_x; ++index_x) {
        double x_coord = m_x_min + one_half * dx + dx * index_x;
        for (std::size_t index_y = 0; index_y < N_y; ++index_y) {
            double y_coord = m_y_min + one_half * dy + dy * index_y;
            MeshGrid.push_back(vector3(x_coord, y_coord, 0.0));
        }
    }
    return MeshGrid;
}

/**
 * @brief Find intersection between a segment and the bounding box
 * The segment MUST have one point inside the bounding box and one outside.
 *
 * @param line_start
 * @param line_end
 * @return box_line_intersect_result = std::optional<std::pair<vector3, box_face>>
 */
box_line_intersect_result bbox::find_box_line_intersection(const vector3& line_start, const vector3& line_end) const {
    // Check X planes
    double t_x_min = (m_x_min - line_start.x()) / (line_end.x() - line_start.x());
    double t_x_max = (m_x_max - line_start.x()) / (line_end.x() - line_start.x());
    // Check Y planes
    double t_y_min = (m_y_min - line_start.y()) / (line_end.y() - line_start.y());
    double t_y_max = (m_y_max - line_start.y()) / (line_end.y() - line_start.y());
    // Check Z planes
    double t_z_min = (m_z_min - line_start.z()) / (line_end.z() - line_start.z());
    double t_z_max = (m_z_max - line_start.z()) / (line_end.z() - line_start.z());

    vector3 intersection_point_x_min = line_start + t_x_min * (line_end - line_start);
    if (intersection_point_x_min.y() >= m_y_min && intersection_point_x_min.y() <= m_y_max && intersection_point_x_min.z() >= m_z_min &&
        intersection_point_x_min.z() <= m_z_max) {
        return std::make_pair(intersection_point_x_min, box_face::x_min);
    }
    vector3 intersection_point_x_max = line_start + t_x_max * (line_end - line_start);
    if (intersection_point_x_max.y() >= m_y_min && intersection_point_x_max.y() <= m_y_max && intersection_point_x_max.z() >= m_z_min &&
        intersection_point_x_max.z() <= m_z_max) {
        return std::make_pair(intersection_point_x_max, box_face::x_max);
    }
    vector3 intersection_point_y_min = line_start + t_y_min * (line_end - line_start);
    if (intersection_point_y_min.x() >= m_x_min && intersection_point_y_min.x() <= m_x_max && intersection_point_y_min.z() >= m_z_min &&
        intersection_point_y_min.z() <= m_z_max) {
        return std::make_pair(intersection_point_y_min, box_face::y_min);
    }
    vector3 intersection_point_y_max = line_start + t_y_max * (line_end - line_start);
    if (intersection_point_y_max.x() >= m_x_min && intersection_point_y_max.x() <= m_x_max && intersection_point_y_max.z() >= m_z_min &&
        intersection_point_y_max.z() <= m_z_max) {
        return std::make_pair(intersection_point_y_max, box_face::y_max);
    }
    vector3 intersection_point_z_min = line_start + t_z_min * (line_end - line_start);
    if (intersection_point_z_min.x() >= m_x_min && intersection_point_z_min.x() <= m_x_max && intersection_point_z_min.y() >= m_y_min &&
        intersection_point_z_min.y() <= m_y_max) {
        return std::make_pair(intersection_point_z_min, box_face::z_min);
    }
    vector3 intersection_point_z_max = line_start + t_z_max * (line_end - line_start);
    if (intersection_point_z_max.x() >= m_x_min && intersection_point_z_max.x() <= m_x_max && intersection_point_z_max.y() >= m_y_min &&
        intersection_point_z_max.y() <= m_y_max) {
        return std::make_pair(intersection_point_z_max, box_face::z_max);
    }
    return std::nullopt;
}

/**
 * @brief Draw a random point inside the triangle.
 *
 * It is done by using the barycentic coordinates
 *
 * @return vector3
 */
vector3 bbox::draw_uniform_random_point_inside_box() const {
    std::mt19937                           random_generator{std::random_device()()};
    std::uniform_real_distribution<double> uniform_unit_distribution{std::uniform_real_distribution<double>{0.0, 1.0}};

    double x_random = uniform_unit_distribution(random_generator);
    double y_random = uniform_unit_distribution(random_generator);
    double z_random = uniform_unit_distribution(random_generator);

    x_random = m_x_min + (m_x_max - m_x_min) * x_random;
    y_random = m_y_min + (m_y_max - m_y_min) * y_random;
    z_random = m_z_min + (m_z_max - m_z_min) * z_random;

    return vector3{x_random, y_random, z_random};
}

}  // namespace mesh

}  // namespace uepm