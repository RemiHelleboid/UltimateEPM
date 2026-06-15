/**
 * @file vector.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief  Common vector class header.
 * @version 0.1
 * @date 2026-05-18
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "float_comparison.hpp"

namespace uepm::common {

enum class permutation_type { XY, XZ, YZ, XYZ, YZX, ZXY };

class vector3 {
 private:
    double m_x;
    double m_y;
    double m_z;

 public:
    constexpr vector3() noexcept : m_x(0.0), m_y(0.0), m_z(0.0) {}
    constexpr vector3(double x, double y) noexcept : m_x(x), m_y(y), m_z(0.0) {}
    constexpr vector3(double x, double y, double z) noexcept : m_x(x), m_y(y), m_z(z) {}

    explicit constexpr vector3(const std::array<double, 2> &array) noexcept : m_x(array[0]), m_y(array[1]), m_z(0.0) {}

    explicit constexpr vector3(const std::array<double, 3> &array) noexcept
        : m_x(array[0]),
          m_y(array[1]),
          m_z(array[2]) {}

    explicit vector3(const std::vector<double> &values) {
        if (values.size() < 3) {
            throw std::invalid_argument("vector3 requires at least 3 values");
        }

        m_x = values[0];
        m_y = values[1];
        m_z = values[2];
    }

    vector3(const vector3 &)                = default;
    vector3(vector3 &&) noexcept            = default;
    vector3 &operator=(const vector3 &)     = default;
    vector3 &operator=(vector3 &&) noexcept = default;
    ~vector3()                              = default;

    constexpr double x() const noexcept { return m_x; }
    constexpr double y() const noexcept { return m_y; }
    constexpr double z() const noexcept { return m_z; }

    constexpr void set_x(double x) noexcept { m_x = x; }
    constexpr void set_y(double y) noexcept { m_y = y; }
    constexpr void set_z(double z) noexcept { m_z = z; }

    constexpr void set_coordinates(double x, double y, double z) noexcept {
        m_x = x;
        m_y = y;
        m_z = z;
    }

    constexpr vector3 to_2d() const noexcept { return vector3{m_x, m_y, 0.0}; }

    constexpr void to_2d_inplace() noexcept { m_z = 0.0; }

    double norm() const noexcept { return std::sqrt(norm_squared()); }

    constexpr double norm_squared() const noexcept { return m_x * m_x + m_y * m_y + m_z * m_z; }

    void re_normalize() noexcept {
        const double v_norm = norm();

        if (v_norm == 0.0) {
            return;
        }

        m_x /= v_norm;
        m_y /= v_norm;
        m_z /= v_norm;
    }

    constexpr void apply_reflection(int x_reflection, int y_reflection, int z_reflection) noexcept {
        m_x *= static_cast<double>(x_reflection);
        m_y *= static_cast<double>(y_reflection);
        m_z *= static_cast<double>(z_reflection);
    }

    void apply_permutation(permutation_type permutation) noexcept {
        switch (permutation) {
            case permutation_type::XY:
                std::swap(m_x, m_y);
                break;

            case permutation_type::XZ:
                std::swap(m_x, m_z);
                break;

            case permutation_type::YZ:
                std::swap(m_y, m_z);
                break;

            case permutation_type::XYZ:
                break;

            case permutation_type::YZX: {
                const double old_x = m_x;
                const double old_y = m_y;
                const double old_z = m_z;

                m_x = old_y;
                m_y = old_z;
                m_z = old_x;
                break;
            }

            case permutation_type::ZXY: {
                const double old_x = m_x;
                const double old_y = m_y;
                const double old_z = m_z;

                m_x = old_z;
                m_y = old_x;
                m_z = old_y;
                break;
            }
        }
    }

    constexpr double dot(const vector3 &rhs) const noexcept { return m_x * rhs.m_x + m_y * rhs.m_y + m_z * rhs.m_z; }

    constexpr vector3 &operator+=(const vector3 &rhs) noexcept {
        m_x += rhs.m_x;
        m_y += rhs.m_y;
        m_z += rhs.m_z;
        return *this;
    }

    constexpr vector3 &operator-=(const vector3 &rhs) noexcept {
        m_x -= rhs.m_x;
        m_y -= rhs.m_y;
        m_z -= rhs.m_z;
        return *this;
    }

    constexpr vector3 &operator*=(double lambda) noexcept {
        m_x *= lambda;
        m_y *= lambda;
        m_z *= lambda;
        return *this;
    }

    constexpr vector3 &operator/=(double lambda) noexcept {
        m_x /= lambda;
        m_y /= lambda;
        m_z /= lambda;
        return *this;
    }

    friend constexpr vector3 operator+(const vector3 &lhs, const vector3 &rhs) noexcept {
        vector3 result = lhs;
        result += rhs;
        return result;
    }

    friend constexpr vector3 operator-(const vector3 &lhs, const vector3 &rhs) noexcept {
        vector3 result = lhs;
        result -= rhs;
        return result;
    }

    friend constexpr vector3 operator*(const vector3 &lhs, double lambda) noexcept {
        vector3 result = lhs;
        result *= lambda;
        return result;
    }

    friend constexpr vector3 operator*(double lambda, const vector3 &rhs) noexcept {
        vector3 result = rhs;
        result *= lambda;
        return result;
    }

    friend constexpr vector3 operator/(const vector3 &lhs, double lambda) noexcept {
        vector3 result = lhs;
        result /= lambda;
        return result;
    }

    friend constexpr double dot(const vector3 &lhs, const vector3 &rhs) noexcept { return lhs.dot(rhs); }

    friend constexpr vector3 middle(const vector3 &lhs, const vector3 &rhs) noexcept {
        constexpr double one_half = 0.5;

        return vector3{
            one_half * (lhs.m_x + rhs.m_x),
            one_half * (lhs.m_y + rhs.m_y),
            one_half * (lhs.m_z + rhs.m_z),
        };
    }

    friend constexpr vector3 point_pair_to_vector(const vector3 &lhs, const vector3 &rhs) noexcept {
        return vector3{
            rhs.m_x - lhs.m_x,
            rhs.m_y - lhs.m_y,
            rhs.m_z - lhs.m_z,
        };
    }

    friend double distance(const vector3 &lhs, const vector3 &rhs) noexcept {
        return point_pair_to_vector(lhs, rhs).norm();
    }

    friend constexpr vector3 cross_product(const vector3 &lhs, const vector3 &rhs) noexcept {
        return vector3{
            lhs.m_y * rhs.m_z - lhs.m_z * rhs.m_y,
            lhs.m_z * rhs.m_x - lhs.m_x * rhs.m_z,
            lhs.m_x * rhs.m_y - lhs.m_y * rhs.m_x,
        };
    }

    friend constexpr double scalar_triple_product(const vector3 &v1, const vector3 &v2, const vector3 &v3) noexcept {
        return v1.dot(cross_product(v2, v3));
    }

    friend double compute_cos_angle(const vector3 &lhs, const vector3 &rhs) noexcept {
        constexpr double epsilon = 1.0e-13;

        const double norm_product = lhs.norm() * rhs.norm();

        if (norm_product < epsilon) {
            return 1.0;
        }

        return lhs.dot(rhs) / norm_product;
    }

    friend bool is_point_between_two_others(const vector3 &a,
                                            const vector3 &b,
                                            const vector3 &point,
                                            double         epsilon = 1.0e-9) noexcept {
        const double d_ab  = distance(a, b);
        const double d_sum = distance(a, point) + distance(point, b);

        return d_ab <= d_sum + epsilon && d_ab >= d_sum - epsilon;
    }

    friend constexpr double double_scalar_product_2d(const vector3 &lhs, const vector3 &rhs) noexcept {
        return lhs.m_x * rhs.m_y - lhs.m_y * rhs.m_x;
    }

    friend std::optional<vector3> compute_line_line_intersection(const vector3 &a,
                                                                 const vector3 &b,
                                                                 const vector3 &c,
                                                                 const vector3 &d) noexcept {
        constexpr double epsilon = 1.0e-14;

        const vector3 r   = b - a;
        const vector3 s   = d - c;
        const vector3 cma = c - a;

        const double rxs = double_scalar_product_2d(r, s);

        if (std::fabs(rxs) < epsilon) {
            return std::nullopt;
        }

        const double t = double_scalar_product_2d(cma, s) / rxs;
        const double u = double_scalar_product_2d(cma, r) / rxs;

        if (t >= -epsilon && t <= 1.0 + epsilon && u >= -epsilon && u <= 1.0 + epsilon) {
            return a + t * r;
        }

        return std::nullopt;
    }

    friend std::ostream &operator<<(std::ostream &os, const vector3 &vect) {
        os << vect.m_x << ',' << vect.m_y << ',' << vect.m_z;
        return os;
    }

    friend bool operator==(const vector3 &lhs, const vector3 &rhs) noexcept {
        return utils::doubles_are_equal(lhs.m_x, rhs.m_x) && utils::doubles_are_equal(lhs.m_y, rhs.m_y) &&
               utils::doubles_are_equal(lhs.m_z, rhs.m_z);
    }

    friend bool operator!=(const vector3 &lhs, const vector3 &rhs) noexcept { return !(lhs == rhs); }
};

static_assert(std::is_trivially_copyable_v<vector3>);
static_assert(std::is_standard_layout_v<vector3>);

}  // namespace uepm::common

namespace uepm::mesh_bz {
using vector3 = uepm::common::vector3;
}

namespace uepm::mesh {
using vector3 = uepm::common::vector3;
}