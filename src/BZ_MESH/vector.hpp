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

#include <cmath>
#include <iostream>
#include <optional>
#include <vector>
namespace bz_mesh {

enum class permutaion_type { XY, XZ, YZ, XYZ, YZX, ZXY };

class vector3 {
 private:
    double m_x;
    double m_y;
    double m_z;

 public:
    vector3() : m_x(0u), m_y(0u), m_z(0u) {}
    vector3(double x, double y) : m_x(x), m_y(y), m_z(0u) {}
    vector3(double x, double y, double z) : m_x(x), m_y(y), m_z(z) {}

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

    double norm() const { return sqrt(m_x * m_x + m_y * m_y + m_z * m_z); }

    void re_normalize() {
        const double v_norm = norm();
        m_x /= v_norm;
        m_y /= v_norm;
        m_z /= v_norm;
    }

    /**
     * @brief Apply a reflection on the vector.
     * For example, apply_reflection(1, -1, 1) will apply the reflexion with respect to the y axis.
     * There is 2^3 = 8 possible reflexions.
     *
     * @param x_reflection
     * @param y_reflection
     * @param z_reflection
     */
    void apply_reflection(int x_reflection, int y_reflection, int z_reflection) {
        m_x = x_reflection * m_x;
        m_y = y_reflection * m_y;
        m_z = z_reflection * m_z;
    }

    /**
     * @brief Apply a permutation on the vector.
     * There is 6 possible permutations.
     *
     * @param permutation
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
                break;
            case permutaion_type::YZX:
                std::swap(m_y, m_z);
                std::swap(m_x, m_y);
                break;
            case permutaion_type::ZXY:
                std::swap(m_z, m_y);
                std::swap(m_x, m_y);
                break;
        }
    }

    double dot(const vector3 &vector_rhs) const { return (m_x * vector_rhs.m_x + m_y * vector_rhs.m_y + m_z * vector_rhs.m_z); }

    friend vector3 middle(const vector3 &vector_lhs, const vector3 &vector_rhs) {
        constexpr double one_half = 1.0 / 2.0;
        return vector3{one_half * (vector_rhs.m_x + vector_lhs.m_x),
                       one_half * (vector_rhs.m_y + vector_lhs.m_y),
                       one_half * (vector_rhs.m_z + vector_lhs.m_z)};
    }

    friend vector3 point_pair_to_vector(const vector3 &vector_lhs, const vector3 &vector_rhs) {
        return vector3{vector_rhs.m_x - vector_lhs.m_x, vector_rhs.m_y - vector_lhs.m_y, vector_rhs.m_z - vector_lhs.m_z};
    }

    friend double distance(const vector3 &vector_lhs, const vector3 &vector_rhs) {
        return point_pair_to_vector(vector_lhs, vector_rhs).norm();
    }

    friend inline vector3 cross_product(const vector3 &vector_lhs, const vector3 &vector_rhs) {
        return vector3(vector_lhs.m_y * vector_rhs.m_z - vector_lhs.m_z * vector_rhs.m_y,
                       vector_lhs.m_z * vector_rhs.m_x - vector_lhs.m_x * vector_rhs.m_z,
                       vector_lhs.m_x * vector_rhs.m_y - vector_lhs.m_y * vector_rhs.m_x);
    }

    friend inline double scalar_triple_product(const vector3 &vector_1, const vector3 &vector_2, const vector3 &vector_3) {
        return vector_1.dot(cross_product(vector_2, vector_3));
    }

    /**
     * @brief Check if a point is on the segment [A, B]
     *
     * @param point_A
     * @param point_B
     * @param point_to_test
     * @return true
     * @return false
     */
    friend inline bool is_point_between_two_others(const vector3 &point_A,
                                                   const vector3 &point_B,
                                                   const vector3 &point_to_test,
                                                   double         epsilon = 1e-9) {
        const double d_AB  = distance(point_A, point_B);
        const double d_sum = distance(point_A, point_to_test) + distance(point_to_test, point_B);
        return (d_AB <= d_sum + epsilon && d_AB >= d_sum - epsilon);
    }

    /**
     * @brief Compute a particular type of scalar product.
     * It is mainly used for the computation of segment-segment intersection.
     * It is equivalent to the determinant of the matrix made by using the tow vectors as columns.
     *
     * @param vector_1
     * @param vector_2
     */
    friend double double_scalar_product_2d(const vector3 &vector_1, const vector3 &vector_2) {
        return vector_1.m_x * vector_2.m_y - vector_1.m_y * vector_2.m_x;
    }

    /**
     * @brief Compute the intersection between two segments [vector_A vector_B] and [vector_C vector_D].
     *
     */
    friend std::optional<vector3> compute_line_line_intersection(const vector3 &vector_A,
                                                                 const vector3 &vector_B,
                                                                 const vector3 &vector_C,
                                                                 const vector3 &vector_D) {
        constexpr double epsilon_intersection   = 1e-14;  // Maybe move that as some kind of global variable.
        const vector3    segment_1              = vector_B - vector_A;
        const vector3    segment_2              = vector_D - vector_C;
        const vector3    starting_point_segment = vector_D - vector_A;

        const double segment_12_scalar_cross_product = double_scalar_product_2d(segment_1, segment_2);
        if (fabs(segment_12_scalar_cross_product) < epsilon_intersection) {
            // std::cout << "Line are colinear : "<< segment_12_scalar_cross_product <<"\n";
            return {};
        }

        const double barycentric_intersection_segment_1 =
            double_scalar_product_2d(starting_point_segment, segment_2) / segment_12_scalar_cross_product;
        const double barycentric_intersection_segment_2 =
            double_scalar_product_2d(starting_point_segment, segment_2) / segment_12_scalar_cross_product;

        if (barycentric_intersection_segment_1 >= -epsilon_intersection && barycentric_intersection_segment_1 <= 1 + epsilon_intersection &&
            barycentric_intersection_segment_2 >= -epsilon_intersection && barycentric_intersection_segment_2 <= 1 + epsilon_intersection) {
            return {vector_A + barycentric_intersection_segment_1 * segment_1};
        }
        return {};
    }

    vector3 &operator+=(const vector3 &vector_lhs) {
        m_x += vector_lhs.m_x;
        m_y += vector_lhs.m_y;
        m_z += vector_lhs.m_z;
        return *this;
    }

    vector3 &operator-=(const vector3 &vector_lhs) {
        m_x -= vector_lhs.m_x;
        m_y -= vector_lhs.m_y;
        m_z -= vector_lhs.m_z;
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

    friend vector3 operator+(const vector3 &vector_lhs, const vector3 &vector_rhs) {
        vector3 SumVector = vector_lhs;
        SumVector += vector_rhs;
        return SumVector;
    }

    friend vector3 operator-(const vector3 &vector_lhs, const vector3 &vector_rhs) {
        vector3 DiffVector = vector_lhs;
        DiffVector -= vector_rhs;
        return DiffVector;
    }

    friend vector3 operator*(const vector3 &vector_lhs, const double lambda) {
        vector3 MultVector = vector_lhs;
        MultVector *= lambda;
        return MultVector;
    }

    friend vector3 operator*(const double lambda, const vector3 &vector_lhs) {
        vector3 MultVector = vector_lhs;
        MultVector *= lambda;
        return MultVector;
    }

    friend vector3 operator/(const vector3 &vector_lhs, const double lambda) {
        vector3 DivVector = vector_lhs;
        DivVector /= lambda;
        return DivVector;
    }

    friend std::ostream &operator<<(std::ostream &os, const vector3 &vect) {
        os << vect.m_x << ',' << vect.m_y << ',' << vect.m_z;
        return os;
    }
};

}  // namespace bz_mesh
