/**
 * @file vector.hpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief Vector class header.
 * @version 0.1
 * @date 2021-07-10
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <optional>
#include <vector>

namespace bz_mesh {

enum class permutaion_type { XY, XZ, YZ, XYZ, YZX, ZXY };

class alignas(32) vector3 {
 private:
    double m_x;
    double m_y;
    double m_z;

 public:
    vector3() : m_x(0.0), m_y(0.0), m_z(0.0) {}
    vector3(double x, double y) : m_x(x), m_y(y), m_z(0.0) {}
    vector3(double x, double y, double z) : m_x(x), m_y(y), m_z(z) {}
    vector3(const vector3 &)            = default;
    vector3 &operator=(const vector3 &) = default;
    vector3(vector3 &&)                 = default;
    vector3 &operator=(vector3 &&)      = default;
    ~vector3()                          = default;

    double x() const { return m_x; }
    double y() const { return m_y; }
    double z() const { return m_z; }

    void set_x(double x) { m_x = x; }
    void set_y(double y) { m_y = y; }
    void set_z(double z) { m_z = z; }
    void set_coordinates(double x, double y, double z) {
        m_x = x;
        m_y = y;
        m_z = z;
    }

    vector3 to_2d() const { return vector3{m_x, m_y, 0.0}; }

    double norm() const { return std::sqrt(m_x * m_x + m_y * m_y + m_z * m_z); }
    double norm_squared() const { return (m_x * m_x + m_y * m_y + m_z * m_z); }

    void re_normalize() {
        const double v_norm = norm();
        if (v_norm == 0.0) return;  // avoid division by zero
        m_x /= v_norm;
        m_y /= v_norm;
        m_z /= v_norm;
    }

    /**
     * @brief Apply a reflection on the vector.
     * (x_reflection, y_reflection, z_reflection) should typically be in {-1, 1}.
     */
    void apply_reflection(int x_reflection, int y_reflection, int z_reflection) {
        m_x = x_reflection * m_x;
        m_y = y_reflection * m_y;
        m_z = z_reflection * m_z;
    }

    /**
     * @brief Apply a permutation on the vector (6 possibilities).
     * XY, XZ, YZ are simple swaps; XYZ is identity; YZX and ZXY are 3-cycles.
     */
    void apply_permutation(permutaion_type permutation) {
        switch (permutation) {
            case permutaion_type::XY:
                std::swap(m_x, m_y);
                break;
            case permutaion_type::XZ:
                std::swap(m_x, m_z);
                break;
            case permutaion_type::YZ:
                std::swap(m_y, m_z);
                break;
            case permutaion_type::XYZ:
                // identity
                break;
            case permutaion_type::YZX: {
                // (x, y, z) -> (y, z, x)
                const double ox = m_x, oy = m_y, oz = m_z;
                m_x = oy;
                m_y = oz;
                m_z = ox;
                break;
            }
            case permutaion_type::ZXY: {
                // (x, y, z) -> (z, x, y)
                const double ox = m_x, oy = m_y, oz = m_z;
                m_x = oz;
                m_y = ox;
                m_z = oy;
                break;
            }
        }
    }

    double dot(const vector3 &rhs) const { return (m_x * rhs.m_x + m_y * rhs.m_y + m_z * rhs.m_z); }

    friend double dot(const vector3 &lhs, const vector3 &rhs) { return (lhs.m_x * rhs.m_x + lhs.m_y * rhs.m_y + lhs.m_z * rhs.m_z); }

    friend vector3 middle(const vector3 &lhs, const vector3 &rhs) {
        constexpr double one_half = 1.0 / 2.0;
        return vector3{one_half * (rhs.m_x + lhs.m_x), one_half * (rhs.m_y + lhs.m_y), one_half * (rhs.m_z + lhs.m_z)};
    }

    friend vector3 point_pair_to_vector(const vector3 &lhs, const vector3 &rhs) {
        return vector3{rhs.m_x - lhs.m_x, rhs.m_y - lhs.m_y, rhs.m_z - lhs.m_z};
    }

    friend double distance(const vector3 &lhs, const vector3 &rhs) { return point_pair_to_vector(lhs, rhs).norm(); }

    friend inline vector3 cross_product(const vector3 &lhs, const vector3 &rhs) {
        return vector3(lhs.m_y * rhs.m_z - lhs.m_z * rhs.m_y, lhs.m_z * rhs.m_x - lhs.m_x * rhs.m_z, lhs.m_x * rhs.m_y - lhs.m_y * rhs.m_x);
    }

    friend inline double scalar_triple_product(const vector3 &v1, const vector3 &v2, const vector3 &v3) {
        return v1.dot(cross_product(v2, v3));
    }

    friend inline double compute_cos_angle(const vector3 &V1, const vector3 &V2) {
        const double dot_product  = V1.dot(V2);
        const double norm_product = V1.norm() * V2.norm();
        return (norm_product < 1.0e-13) ? 1.0 : dot_product / norm_product;
    }

    friend inline bool is_point_between_two_others(const vector3 &A, const vector3 &B, const vector3 &P, double epsilon = 1e-9) {
        const double d_AB  = distance(A, B);
        const double d_sum = distance(A, P) + distance(P, B);
        return (d_AB <= d_sum + epsilon && d_AB >= d_sum - epsilon);
    }

    friend double double_scalar_product_2d(const vector3 &v1, const vector3 &v2) {
        // 2D determinant using (x,y) components
        return v1.m_x * v2.m_y - v1.m_y * v2.m_x;
    }

    /**
     * @brief Compute the intersection between two 2D segments [A B] and [C D] (using x,y only).
     * Returns std::nullopt if parallel/colinear or if intersection lies outside either segment.
     */
    friend std::optional<vector3> compute_line_line_intersection(const vector3 &A, const vector3 &B, const vector3 &C, const vector3 &D) {
        constexpr double eps = 1e-14;
        const vector3    r   = B - A;  // segment 1
        const vector3    s   = D - C;  // segment 2
        const vector3    qmp = C - A;  // from A to C

        const double rxs = double_scalar_product_2d(r, s);
        if (std::fabs(rxs) < eps) {
            // Parallel (or colinear) in 2D
            return {};
        }

        // Barycentric parameters:
        // A + t*r  intersects  C + u*s
        const double t = double_scalar_product_2d(qmp, s) / rxs;
        const double u = double_scalar_product_2d(qmp, r) / rxs;

        if (t >= -eps && t <= 1.0 + eps && u >= -eps && u <= 1.0 + eps) {
            const vector3 P = A + t * r;
            return {P};
        }
        return {};
    }

    vector3 &operator+=(const vector3 &lhs) {
        m_x += lhs.m_x;
        m_y += lhs.m_y;
        m_z += lhs.m_z;
        return *this;
    }

    vector3 &operator-=(const vector3 &lhs) {
        m_x -= lhs.m_x;
        m_y -= lhs.m_y;
        m_z -= lhs.m_z;
        return *this;
    }

    vector3 &operator*=(const double lambda) {
        m_x *= lambda;
        m_y *= lambda;
        m_z *= lambda;
        return *this;
    }

    vector3 &operator/=(const double lambda) {
        m_x /= lambda;
        m_y /= lambda;
        m_z /= lambda;
        return *this;
    }

    friend vector3 operator+(const vector3 &lhs, const vector3 &rhs) {
        vector3 v = lhs;
        v += rhs;
        return v;
    }

    friend vector3 operator-(const vector3 &lhs, const vector3 &rhs) {
        vector3 v = lhs;
        v -= rhs;
        return v;
    }

    friend vector3 operator*(const vector3 &lhs, const double lambda) {
        vector3 v = lhs;
        v *= lambda;
        return v;
    }

    friend vector3 operator*(const double lambda, const vector3 &lhs) {
        vector3 v = lhs;
        v *= lambda;
        return v;
    }

    friend vector3 operator/(const vector3 &lhs, const double lambda) {
        vector3 v = lhs;
        v /= lambda;
        return v;
    }

    friend std::ostream &operator<<(std::ostream &os, const vector3 &vect) {
        os << vect.m_x << ',' << vect.m_y << ',' << vect.m_z;
        return os;
    }
};

static_assert(std::is_trivially_copyable_v<vector3>, "vector3 must be trivially copyable");
static_assert(std::is_standard_layout_v<vector3>, "vector3 must be standard-layout");

}  // namespace bz_mesh