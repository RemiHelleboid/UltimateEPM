/**
 * @file sphere.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-08-15
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <vector>

#include "bbox.hpp"
#include "vector_mesh.hpp"

namespace uepm::mesh {

class sphere {
 private:
    vector3 m_center;
    double  m_radius;

 public:
    sphere() = default;
    sphere(const vector3 &center, double radius) : m_center(center), m_radius(radius) {}
    sphere(const sphere &other)                = default;
    sphere(sphere &&other) noexcept            = default;
    sphere &operator=(const sphere &other)     = default;
    sphere &operator=(sphere &&other) noexcept = default;
    ~sphere()                                  = default;

    vector3 center() const { return m_center; }
    double  radius() const { return m_radius; }
    void    set_center(const vector3 &center) { m_center = center; }
    void    set_radius(double radius) { m_radius = radius; }
    bool    is_inside(const vector3 &point) const { return (point - m_center).norm() < m_radius; }
    bool    is_overlapping(const bbox &box) const {
        // Check if the sphere is overlapping with the bounding box
        double x_min = box.get_x_min();
        double x_max = box.get_x_max();
        double y_min = box.get_y_min();
        double y_max = box.get_y_max();
        double z_min = box.get_z_min();
        double z_max = box.get_z_max();

        // Find the closest point on the bounding box to the center of the sphere
        double closest_x = std::max(x_min, std::min(m_center.x(), x_max));
        double closest_y = std::max(y_min, std::min(m_center.y(), y_max));
        double closest_z = std::max(z_min, std::min(m_center.z(), z_max));

        // Calculate the distance from the center of the sphere to the closest point
        vector3 closest_point(closest_x, closest_y, closest_z);
        double  distance_squared = (closest_point - m_center).norm_squared();

        // The sphere overlaps with the bounding box if the distance is less than or equal to the radius squared
        return distance_squared <= (m_radius * m_radius);
    }
};

}  // namespace uepm::mesh