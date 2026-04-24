#pragma once
#include <array>
#include <optional>
#include <random>
#include <utility>
#include <vector>

#include "vector.hpp"

namespace uepm {

namespace mesh {

class element;

enum class box_face { x_min, x_max, y_min, y_max, z_min, z_max };

using box_line_intersect_result = std::optional<std::pair<vector3, box_face>>;

class bbox {
 private:
    double m_x_min{0.0};
    double m_x_max{0.0};
    double m_y_min{0.0};
    double m_y_max{0.0};
    double m_z_min{0.0};
    double m_z_max{0.0};

 public:
    bbox() = default;
    bbox(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max)
        : m_x_min(x_min),
          m_x_max(x_max),
          m_y_min(y_min),
          m_y_max(y_max),
          m_z_min(z_min),
          m_z_max(z_max) {
        if (!check_order()) {
            throw std::invalid_argument("Order of bbox corner is bad.");
        }
    }

    bbox(const vector3 &bottom_left_front_corner, const vector3 &up_right_back_corner)
        : m_x_min(bottom_left_front_corner.x()),
          m_x_max(up_right_back_corner.x()),
          m_y_min(bottom_left_front_corner.y()),
          m_y_max(up_right_back_corner.y()),
          m_z_min(bottom_left_front_corner.z()),
          m_z_max(up_right_back_corner.z()) {
        if (!check_order()) {
            throw std::invalid_argument("Order of bbox corner is bad.");
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

    double get_surface() const { return fabs(m_x_max - m_x_min) * fabs(m_y_max - m_y_min); }

    double get_volume() const { return fabs(m_x_max - m_x_min) * fabs(m_y_max - m_y_min) * fabs(m_z_max - m_z_min); }

    vector3 get_center() const {
        const double one_half = 0.5;
        return vector3{one_half * (m_x_min + m_x_max), one_half * (m_y_min + m_y_max), one_half * (m_z_min + m_z_max)};
    }

    bool              is_inside(const vector3 &location) const;
    bool              is_inside_2d(const vector3 &location) const;
    std::vector<bbox> split_2d_box_in_quadrants() const;
    std::vector<bbox> split_3d_box_in_octants() const;

    void dilate(double factor) {
        m_x_min *= factor;
        m_x_max *= factor;
        m_y_min *= factor;
        m_y_max *= factor;
        m_z_min *= factor;
        m_z_max *= factor;
    }

    bool is_overlapping(const bbox &second_box) const;
    bool is_overlapping_2d(const bbox &second_box) const;
    bool is_overlapping_triangle(const element &triangle) const;
    bool is_overlapping_tetra(const element &tetra) const;
    bool is_box_overlapping(const element &my_element) const;
    bool is_box_overlapping_2d(const element &my_element) const;
    bool is_overlapping_circle_2d(const vector3 &center, double radius) const;
    bool is_overlapping_sphere_3d(const vector3 &center, double radius) const;

    box_line_intersect_result find_box_line_intersection(const vector3 &line_start, const vector3 &line_end) const;

    std::vector<vector3> generate_mesh_grid_2d(std::size_t N_x, std::size_t N_y) const;
    std::vector<vector3> generate_mesh_grid_3d(std::size_t N_x, std::size_t N_y, std::size_t N_z) const;
    std::vector<vector3> generate_inner_mesh_grid_2d(std::size_t N_x, std::size_t N_y) const;
    std::vector<vector3> generate_inner_mesh_grid_3d(std::size_t N_x, std::size_t N_y, std::size_t N_z) const;

    vector3              draw_uniform_random_point_inside_box() const;
    friend std::ostream &operator<<(std::ostream &os, const bbox &my_box);

    template <std::uniform_random_bit_generator Generator>
    vector3 draw_uniform_random_point_inside_box(Generator &random_generator) const {
        std::uniform_real_distribution<double> uniform_unit_distribution(0.0, 1.0);

        double x_random = uniform_unit_distribution(random_generator);
        double y_random = uniform_unit_distribution(random_generator);
        double z_random = uniform_unit_distribution(random_generator);

        x_random = m_x_min + (m_x_max - m_x_min) * x_random;
        y_random = m_y_min + (m_y_max - m_y_min) * y_random;
        z_random = m_z_min + (m_z_max - m_z_min) * z_random;

        return vector3{x_random, y_random, z_random};
    }
};

}  // namespace mesh

}  // namespace uepm