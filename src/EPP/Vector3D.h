#pragma once

#include <cmath>
#include <iostream>
#include <tuple>

template <typename T>
class Vector3D {
 public:
    T X;
    T Y;
    T Z;

    Vector3D();
    Vector3D(const T& v);
    template <typename O>
    Vector3D(const Vector3D<O>& other);
    explicit Vector3D(const T& x, const T& y, const T& z);

    template <typename O>
    Vector3D& operator=(const Vector3D<O>& other);

    const Vector3D& operator+() const;
    Vector3D        operator-() const;

    template <typename O>
    Vector3D operator+(const Vector3D<O>& other) const;
    template <typename O>
    Vector3D operator-(const Vector3D<O>& other) const;
    template <typename O>
    T operator*(const Vector3D<O>& other) const;
    template <typename O>
    Vector3D operator%(const Vector3D<O>& other) const;

    Vector3D operator*(T s) const;
    Vector3D operator/(T s) const;

    template <typename O>
    Vector3D& operator+=(const Vector3D<O>& other);
    template <typename O>
    Vector3D& operator-=(const Vector3D<O>& other);
    template <typename O>
    T operator*=(const Vector3D<O>& other);
    template <typename O>
    Vector3D& operator%=(const Vector3D<O>& other);

    Vector3D& operator*=(T s);
    Vector3D& operator/=(T s);

    double   Length() const;
    Vector3D Normalize() const;

    template <typename O, typename A>
    Vector3D RotateAround(const Vector3D<O>& other, A angle) const;
    template <typename O, typename A>
    Vector3D RotateTowards(const Vector3D<O>& other, A angle) const;

    double getTheta() const {
        double cosTheta = Z / Length();
        if (std::isnan(cosTheta) || std::isinf(cosTheta) || cosTheta > 1. || cosTheta < -1) cosTheta = (cosTheta < 0 ? -1 : 1);

        return std::acos(cosTheta);
    }

    double getPhi() const { return atan2(Y, X); }

    /**
     * @brief Compute the cosinus of the angle between two vectors.
     *
     * @param V1
     * @param V2
     * @return double
     */
    friend double compte_cos_angle(const Vector3D& V1, const Vector3D& V2) {
        double dot_product  = V1 * V2;
        double norm_product = V1.Length() * V2.Length();
        return (norm_product < 1.0e-13) ? 1.0 : dot_product / (V1.Length() * V2.Length());
    }

    friend std::ostream& operator<<(std::ostream& os, const Vector3D& v) {
        os << "(" << v.X << ", " << v.Y << ", " << v.Z << ")";
        return os;
    }
};

template <typename T>
Vector3D<T> operator*(T o, const Vector3D<T>& t) {
    return t * o;
}

template <typename T>
Vector3D<T> operator/(T o, const Vector3D<T>& t) {
    return t / o;
}

template <typename T>
Vector3D<T> cross_product(const Vector3D<T>& v1, const Vector3D<T>& v2) {
    return Vector3D<T>(v1.Y * v2.Z - v1.Z * v2.Y, v1.Z * v2.X - v1.X * v2.Z, v1.X * v2.Y - v1.Y * v2.Z);
}

template <typename T>
bool operator==(const Vector3D<T>& f, const Vector3D<T>& t) {
    return f.X == t.X && f.Y == t.Y && f.Z == t.Z;
}


template <typename T>
bool operator<(const Vector3D<T>& lhs, const Vector3D<T>& rhs) {
    return std::tie(lhs.X, lhs.Y, lhs.Z) < std::tie(rhs.X, rhs.Y, rhs.Z);
}

#ifndef _VECTOR_3D_IMPL
#include "Vector3D.inl"
#endif