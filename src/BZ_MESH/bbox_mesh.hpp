/**
 * @file box_bz.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-08-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <memory>
#include <random>
#include <vector>

#include "vector.hpp"

namespace uepm::mesh_bz {

/**
 * @brief Class representing a box in 3D space.
 * The box is of the form AABB (Axis-Aligned Bounding Box)
 *
 *
 */
class bbox_mesh {
 private:
    double m_x_min{0.0};
    double m_x_max{0.0};
    double m_y_min{0.0};
    double m_y_max{0.0};
    double m_z_min{0.0};
    double m_z_max{0.0};

 public:
    bbox_mesh(){};
    bbox_mesh(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max)
        : m_x_min(x_min),
          m_x_max(x_max),
          m_y_min(y_min),
          m_y_max(y_max),
          m_z_min(z_min),
          m_z_max(z_max) {
        if (!check_order()) {
            throw std::invalid_argument("Order of bbox_mesh corner is bad.");
        }
    }

    bbox_mesh(const vector3 &bottom_left_front_corner, const vector3 &up_right_back_corner)
        : m_x_min(bottom_left_front_corner.x()),
          m_x_max(up_right_back_corner.x()),
          m_y_min(bottom_left_front_corner.y()),
          m_y_max(up_right_back_corner.y()),
          m_z_min(bottom_left_front_corner.z()),
          m_z_max(up_right_back_corner.z()) {
        if (!check_order()) {
            throw std::invalid_argument("Order of bbox_mesh corner is bad.");
        }
    }

    bool check_order() const { return (m_x_min <= m_x_max) && (m_y_min <= m_y_max) && (m_z_min <= m_z_max); }

    double get_x_min() const { return m_x_min; }
    double get_x_max() const { return m_x_max; }
    double get_y_min() const { return m_y_min; }
    double get_y_max() const { return m_y_max; }
    double get_z_min() const { return m_z_min; }
    double get_z_max() const { return m_z_max; }

    double get_x_size() const { return m_x_max - m_x_min; }
    double get_y_size() const { return m_y_max - m_y_min; }
    double get_z_size() const { return m_z_max - m_z_min; }

    double get_diagonal_size() const {
        return sqrt(get_x_size() * get_x_size() + get_y_size() * get_y_size() + get_z_size() * get_z_size());
    }

    double get_volume() const {
        const double surface = fabs(m_x_max - m_x_min) * fabs(m_y_max - m_y_min) * fabs(m_z_max - m_z_min);
        return surface;
    }

    vector3 get_center() const {
        constexpr double one_half = 1.0 / 2.0;
        return {one_half * (m_x_min + m_x_max), one_half * (m_y_min + m_y_max), one_half * (m_z_min + m_z_max)};
    }

    bool is_inside(const vector3 &location) const {
        return (location.x() > m_x_min) && (location.x() < m_x_max) && (location.y() > m_y_min) && (location.y() < m_y_max) &&
               (location.z() > m_z_min) && (location.z() < m_z_max);
    }

    /**
     * @brief Split the box into 8 equal sub-boxes.
     *
     * @return std::array<bbox_mesh, 8>
     */
    std::array<bbox_mesh, 8> split_3d_box_in_octants() const {
        std::array<bbox_mesh, 8> octants;
        const vector3         box_center = get_center();

        octants[0] = bbox_mesh(m_x_min, box_center.x(), m_y_min, box_center.y(), m_z_min, box_center.z());
        octants[1] = bbox_mesh(box_center.x(), m_x_max, m_y_min, box_center.y(), m_z_min, box_center.z());
        octants[2] = bbox_mesh(box_center.x(), m_x_max, box_center.y(), m_y_max, m_z_min, box_center.z());
        octants[3] = bbox_mesh(m_x_min, box_center.x(), box_center.y(), m_y_max, m_z_min, box_center.z());
        octants[4] = bbox_mesh(m_x_min, box_center.x(), m_y_min, box_center.y(), box_center.z(), m_z_max);
        octants[5] = bbox_mesh(box_center.x(), m_x_max, m_y_min, box_center.y(), box_center.z(), m_z_max);
        octants[6] = bbox_mesh(box_center.x(), m_x_max, box_center.y(), m_y_max, box_center.z(), m_z_max);
        octants[7] = bbox_mesh(m_x_min, box_center.x(), box_center.y(), m_y_max, box_center.z(), m_z_max);
        return octants;
    }

    void dilate(double factor) {
        m_x_min *= factor;
        m_x_max *= factor;
        m_y_min *= factor;
        m_y_max *= factor;
        m_z_min *= factor;
        m_z_max *= factor;
    }

    void translate(const vector3 &translation) {
        m_x_min += translation.x();
        m_x_max += translation.x();
        m_y_min += translation.y();
        m_y_max += translation.y();
        m_z_min += translation.z();
        m_z_max += translation.z();
    }

    bool is_overlapping(const bbox_mesh &second_box) const {
        const bool noOverlap = this->m_x_min > second_box.m_x_max || second_box.m_x_min > this->m_x_max ||
                               this->m_y_min > second_box.m_y_max || second_box.m_y_min > this->m_y_max ||
                               this->m_z_min > second_box.m_z_max || second_box.m_z_min > this->m_z_max;
        return !noOverlap;
    }

    friend std::ostream &operator<<(std::ostream &os, const bbox_mesh &my_box) {
        os << "bbox_mesh: x_min = " << my_box.m_x_min << " x_max = " << my_box.m_x_max << " y_min = " << my_box.m_y_min
           << " y_max = " << my_box.m_y_max << " z_min = " << my_box.m_z_min << " z_max = " << my_box.m_z_max;
        return os;
    }
};

}  // namespace uepm::mesh_bz