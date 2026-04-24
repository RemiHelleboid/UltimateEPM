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

#include <memory>
#include <vector>

#include "bbox.hpp"
#include "vector3.hpp"

namespace uepm {

class sphere {
 private:
    vector3 m_center;
    double  m_radius;

 public:
    sphere(){};
    sphere(const vector3 &center, double radius);
    sphere(const sphere &other);
    sphere(sphere &&other) noexcept;
    sphere &operator=(const sphere &other);
    sphere &operator=(sphere &&other) noexcept;
    ~sphere()       = default;

    vector3 center() const;
    double  radius() const;
    void    set_center(const vector3 &center) { m_center = center; }
    void    set_radius(double radius) { m_radius = radius; }
    bool    is_inside(const vector3 &point) const { return (point - m_center).norm() < m_radius; }
    bool    is_overlapping(const bbox &box) const bool is_overlapping(const sphere &other) const {
        return (other.m_center - m_center).norm() < m_radius + other.m_radius;
    }
};

}  // namespace uepm