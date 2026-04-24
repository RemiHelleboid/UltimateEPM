/**
 * @file grid_mesh_interpolator.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-08-02
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "vector.hpp"
#include "vertex.hpp"

namespace uepm {

namespace mesh {

class quadrangle_2d {
 private:
    std::array<vector3, 4> m_vertices;
    std::array<double, 4>  m_scalar_values = {0.0, 0.0, 0.0, 0.0};
    std::array<vector3, 4> m_vector_values {vector3{}, vector3{}, vector3{}, vector3{}};

 public:
    quadrangle_2d() = delete;
    quadrangle_2d(const vector3 &v0, const vector3 &v1, const vector3 &v2, const vector3 &v3) : m_vertices{v0, v1, v2, v3} {}
    explicit quadrangle_2d(const std::array<vector3, 4> &vertices) : m_vertices{vertices} {}

    quadrangle_2d(const quadrangle_2d &other) = default;
    quadrangle_2d(quadrangle_2d &&other)      = default;
    quadrangle_2d &operator=(const quadrangle_2d &other) = default;
    quadrangle_2d &operator=(quadrangle_2d &&other) = default;
    ~quadrangle_2d()                                = default;

    void set_scalar_field_values(const std::array<double, 4> &values) { m_scalar_values = values; }
    void set_scalar_field_values(double v0, double v1, double v2, double v3) { m_scalar_values = {v0, v1, v2, v3}; }

    void set_vector_field_values(const std::array<vector3, 4> &values) { m_vector_values = values; }
    void set_vector_field_values(const vector3 &v0, const vector3 &v1, const vector3 &v2, const vector3 &v3) {
        m_vector_values = {v0, v1, v2, v3};
    }

    bool check_im_a_rectangle() const {
        const double length_0 = (m_vertices[1] - m_vertices[0]).norm();
        const double length_1 = (m_vertices[3] - m_vertices[2]).norm();
        const double height_0 = (m_vertices[1] - m_vertices[2]).norm();
        const double height_1 = (m_vertices[3] - m_vertices[0]).norm();
        return (length_0 == length_1) && (height_0 == height_1);
    }

    double get_area() const {
        const double length = (m_vertices[1] - m_vertices[0]).norm();
        const double height = (m_vertices[1] - m_vertices[2]).norm();
        return length * height;
    }

    const std::array<vector3, 4> &get_vertices() const { return m_vertices; }
    const vector3 &               get_vertex(int i) const { return m_vertices[i]; }

    double                       get_scalar_value(int i) const { return m_scalar_values[i]; }
    const std::array<double, 4> &get_scalar_values() const { return m_scalar_values; }

    const vector3 &               get_vector_value(int i) const { return m_vector_values[i]; }
    const std::array<vector3, 4> &get_vector_values() const { return m_vector_values; }

    const vector3 &get_nearest_vertex(const vector3 &point) const {
        double min_distance = std::numeric_limits<double>::max();
        int    min_index    = 0;
        for (int i = 0; i < 4; ++i) {
            const double distance = (m_vertices[i] - point).norm();
            if (distance < min_distance) {
                min_distance = distance;
                min_index    = i;
            }
        }
        return m_vertices[min_index];
    }

    int get_nearest_vertex_index(const vector3 &point) const {
        double min_distance = std::numeric_limits<double>::max();
        int    min_index    = 0;
        for (int i = 0; i < 4; ++i) {
            const double distance = (m_vertices[i] - point).norm();
            if (distance < min_distance) {
                min_distance = distance;
                min_index    = i;
            }
        }
        return min_index;
    }

    double get_length() const { return (m_vertices[1] - m_vertices[0]).norm(); }
    double get_height() const { return (m_vertices[3] - m_vertices[0]).norm(); }

    /**
     * @brief Nearest neighbor interpolation within a quadrangle for a scalar field.
     *
     * @param location
     * @return double
     */
    double nearest_neighbor_scalar_interpolation(const vector3 &location) const {
        int index_nearest_vertex = get_nearest_vertex_index(location);
        return m_scalar_values[index_nearest_vertex];
    }

    /**
     * @brief Nearest neighbor interpolation within a quadrangle for a vector field.
     *
     * @param location
     * @return vector3
     */
    vector3 nearest_neighbor_vector_interpolation(const vector3 &location) const {
        int index_nearest_vertex = get_nearest_vertex_index(location);
        return m_vector_values[index_nearest_vertex];
    }

    /**
     * @brief Bilinear interpolation within a quadrangle for a scalar field.
     *
     * https://fr.wikipedia.org/wiki/Interpolation_bilin%C3%A9aire
     *
     * Defined with values at vertices as argument because it will be reused for vector field interpolation.
     *
     * @param location
     * @param v0
     * @param v1
     * @param v2
     * @param v3
     * @return double
     */
    double bilinear_interpolation_scalar_at_location(const vector3 &location, double v0, double v1, double v2, double v3) const {
        const double dx = location.x() - m_vertices[0].x();
        const double dy = location.y() - m_vertices[0].y();

        const double delta_x = m_vertices[1].x() - m_vertices[0].x();
        const double delta_y = m_vertices[3].y() - m_vertices[0].y();

        const double delta_fx  = v1 - v0;
        const double delta_fy  = v3 - v0;
        const double delta_fxy = v0 + v2 - v1 - v3;

        const double f_xy = delta_fx * (dx / delta_x) + delta_fy * (dy / delta_y) + delta_fxy * (dx * dy / (delta_x * delta_y)) + v0;

        return f_xy;
    }

    /**
     * @brief Bilinear interpolation of the vector field at a given location.
     *
     * This is not the most efficient way to do it, but it is simple and works (for now).
     * Indeed, many things are recalculated (delta_x,y dx,y).
     *
     * @param location
     * @return double
     */
    double bilinear_interpolation_scalar_at_location(const vector3 &location) {
        return bilinear_interpolation_scalar_at_location(location,
                                                         m_scalar_values[0],
                                                         m_scalar_values[1],
                                                         m_scalar_values[2],
                                                         m_scalar_values[3]);
    }

    vector3 bilinear_interpolation_vector_at_location(const vector3 &location) {
        const double inter_x = bilinear_interpolation_scalar_at_location(location,
                                                                         m_vector_values[0].x(),
                                                                         m_vector_values[1].x(),
                                                                         m_vector_values[2].x(),
                                                                         m_vector_values[3].x());
        const double inter_y = bilinear_interpolation_scalar_at_location(location,
                                                                         m_vector_values[0].y(),
                                                                         m_vector_values[1].y(),
                                                                         m_vector_values[2].y(),
                                                                         m_vector_values[3].y());
        const double inter_z = bilinear_interpolation_scalar_at_location(location,
                                                                         m_vector_values[0].z(),
                                                                         m_vector_values[1].z(),
                                                                         m_vector_values[2].z(),
                                                                         m_vector_values[3].z());
        return vector3(inter_x, inter_y, inter_z);
    }
};

class quadrangle_3d {
 private:
    std::array<vector3, 8> m_vertices;
    std::array<double, 8>  m_scalar_values = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                              0.0, 0.0};
    std::array<vector3, 8> m_vector_values {vector3{}, vector3{}, vector3{}, vector3{}, vector3{}, vector3{}, vector3{}, vector3{}};

 public:
    quadrangle_3d() = default;
    quadrangle_3d(const vector3 &v0,
                  const vector3 &v1,
                  const vector3 &v2,
                  const vector3 &v3,
                  const vector3 &v4,
                  const vector3 &v5,
                  const vector3 &v6,
                  const vector3 &v7)
        : m_vertices{v0, v1, v2, v3, v4, v5, v6, v7} {}
    explicit quadrangle_3d(const std::array<vector3, 8> &vertices) : m_vertices{vertices} {}
    quadrangle_3d(const quadrangle_3d &other) = default;
    quadrangle_3d(quadrangle_3d &&other)      = default;
    quadrangle_3d &operator=(const quadrangle_3d &other) = default;
    quadrangle_3d &operator=(quadrangle_3d &&other) = default;
    ~quadrangle_3d()                                = default;
    void set_scalar_field_values(const std::array<double, 8> &values) { m_scalar_values = values; }
    void set_scalar_field_values(double v0, double v1, double v2, double v3, double v4, double v5, double v6, double v7) {
        m_scalar_values = {v0, v1, v2, v3, v4, v5, v6, v7};
    }
    void set_vector_field_values(const std::array<vector3, 8> &values) { m_vector_values = values; }
    void set_vector_field_values(const vector3 &v0,
                                 const vector3 &v1,
                                 const vector3 &v2,
                                 const vector3 &v3,
                                 const vector3 &v4,
                                 const vector3 &v5,
                                 const vector3 &v6,
                                 const vector3 &v7) {
        m_vector_values = {v0, v1, v2, v3, v4, v5, v6, v7};
    }

    vector3 get_nearest_vertex_index(const vector3 &location) const {
        double  min_distance = std::numeric_limits<double>::max();
        vector3 nearest_vertex;
        for (const auto &vertex : m_vertices) {
            const double distance = (vertex - location).norm();
            if (distance < min_distance) {
                min_distance   = distance;
                nearest_vertex = vertex;
            }
        }
        return nearest_vertex;
    }

    std::size_t get_nearest_vertex_index(const vector3 &location, double &distance) const {
        distance                         = std::numeric_limits<double>::max();
        std::size_t nearest_vertex_index = 0;
        for (std::size_t i = 0; i < m_vertices.size(); ++i) {
            const double current_distance = (m_vertices[i] - location).norm();
            if (current_distance < distance) {
                distance             = current_distance;
                nearest_vertex_index = i;
            }
        }
        return nearest_vertex_index;
    }

    double nearest_neighbor_scalar_interpolation(const vector3 &location) const {
        double min_distance  = std::numeric_limits<double>::max();
        double nearest_value = 0.0;
        for (std::size_t i = 0; i < m_vertices.size(); ++i) {
            const double distance = (m_vertices[i] - location).norm();
            if (distance < min_distance) {
                min_distance  = distance;
                nearest_value = m_scalar_values[i];
            }
        }
        return nearest_value;
    }

    vector3 nearest_neighbor_vector_interpolation(const vector3 &location) const {
        double  min_distance = std::numeric_limits<double>::max();
        vector3 nearest_value{0.0, 0.0, 0.0};
        for (std::size_t i = 0; i < m_vertices.size(); ++i) {
            const double distance = (m_vertices[i] - location).norm();
            if (distance < min_distance) {
                min_distance  = distance;
                nearest_value = m_vector_values[i];
            }
        }
        return nearest_value;
    }

    /**
     * @brief Trilinear (3D) interpolation of the scalar field at a given location.
     *
     * The interpolation is done on the 8 vertices of the quadrangle.
     * The values a_0, a_x, a_y, a_z, a_xy, a_xz, a_yz, a_xyz are computed from a python script using sympy.
     * (It is stored in the file "trilinear_interpolation.py" in the scripts folder).
     *
     * @param location
     * @param v0
     * @param v1
     * @param v2
     * @param v3
     * @param v4
     * @param v5
     * @param v6
     * @param v7
     * @return double
     */
    double trilinear_interpolation_scalar_at_location(const vector3 &location,
                                                      double         v0,
                                                      double         v1,
                                                      double         v2,
                                                      double         v3,
                                                      double         v4,
                                                      double         v5,
                                                      double         v6,
                                                      double         v7) const {
        const double x1 = m_vertices[0].x();
        const double x2 = m_vertices[1].x();
        const double y1 = m_vertices[0].y();
        const double y2 = m_vertices[3].y();
        const double z1 = m_vertices[0].z();
        const double z2 = m_vertices[4].z();

        const double minus_volume = ((x1 - x2) * (y1 - y2) * (z1 - z2));

        const double a_0 = (-v0 * x2 * y2 * z2 + v1 * x1 * y2 * z2 - v2 * x1 * y1 * z2 + v3 * x2 * y1 * z2 + v4 * x2 * y2 * z1 -
                            v5 * x1 * y2 * z1 + v6 * x1 * y1 * z1 - v7 * x2 * y1 * z1) /
                           minus_volume;
        const double a_x =
            -(-v0 * y2 * z2 + v1 * y2 * z2 - v2 * y1 * z2 + v3 * y1 * z2 + v4 * y2 * z1 - v5 * y2 * z1 + v6 * y1 * z1 - v7 * y1 * z1) /
            minus_volume;
        const double a_y =
            (v0 * x2 * z2 - v1 * x1 * z2 + v2 * x1 * z2 - v3 * x2 * z2 - v4 * x2 * z1 + v5 * x1 * z1 - v6 * x1 * z1 + v7 * x2 * z1) /
            minus_volume;
        const double a_z =
            (v0 * x2 * y2 - v1 * x1 * y2 + v2 * x1 * y1 - v3 * x2 * y1 - v4 * x2 * y2 + v5 * x1 * y2 - v6 * x1 * y1 + v7 * x2 * y1) /
            minus_volume;
        const double a_xy  = (-v0 * z2 + v1 * z2 - v2 * z2 + v3 * z2 + v4 * z1 - v5 * z1 + v6 * z1 - v7 * z1) / minus_volume;
        const double a_xz  = -(v0 * y2 - v1 * y2 + v2 * y1 - v3 * y1 - v4 * y2 + v5 * y2 - v6 * y1 + v7 * y1) / minus_volume;
        const double a_yz  = (-v0 * x2 + v1 * x1 - v2 * x1 + v3 * x2 + v4 * x2 - v5 * x1 + v6 * x1 - v7 * x2) / minus_volume;
        const double a_xyz = (v0 - v1 + v2 - v3 - v4 + v5 - v6 + v7) / minus_volume;

        const double f_xyz = a_0 + a_x * location.x() + a_y * location.y() + a_z * location.z() + a_xy * location.x() * location.y() +
                             a_xz * location.x() * location.z() + a_yz * location.y() * location.z() +
                             a_xyz * location.x() * location.y() * location.z();

        return f_xyz;
    }

    double trilinear_interpolation_scalar_at_location(const vector3 &location) const {
        return trilinear_interpolation_scalar_at_location(location,
                                                          m_scalar_values[0],
                                                          m_scalar_values[1],
                                                          m_scalar_values[2],
                                                          m_scalar_values[3],
                                                          m_scalar_values[4],
                                                          m_scalar_values[5],
                                                          m_scalar_values[6],
                                                          m_scalar_values[7]);
    }

    vector3 trilinear_interpolation_vector_at_location(const vector3 &location) {
        const double interp_Vx = trilinear_interpolation_scalar_at_location(location,
                                                                            m_vector_values[0].x(),
                                                                            m_vector_values[1].x(),
                                                                            m_vector_values[2].x(),
                                                                            m_vector_values[3].x(),
                                                                            m_vector_values[4].x(),
                                                                            m_vector_values[5].x(),
                                                                            m_vector_values[6].x(),
                                                                            m_vector_values[7].x());

        const double interp_Vy = trilinear_interpolation_scalar_at_location(location,
                                                                            m_vector_values[0].y(),
                                                                            m_vector_values[1].y(),
                                                                            m_vector_values[2].y(),
                                                                            m_vector_values[3].y(),
                                                                            m_vector_values[4].y(),
                                                                            m_vector_values[5].y(),
                                                                            m_vector_values[6].y(),
                                                                            m_vector_values[7].y());

        const double interp_Vz = trilinear_interpolation_scalar_at_location(location,
                                                                            m_vector_values[0].z(),
                                                                            m_vector_values[1].z(),
                                                                            m_vector_values[2].z(),
                                                                            m_vector_values[3].z(),
                                                                            m_vector_values[4].z(),
                                                                            m_vector_values[5].z(),
                                                                            m_vector_values[6].z(),
                                                                            m_vector_values[7].z());

        return vector3(interp_Vx, interp_Vy, interp_Vz);
    }
};

}  // namespace mesh

}  // namespace uepm