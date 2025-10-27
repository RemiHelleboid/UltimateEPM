/**
 * @file numerical_helper.hpp
 *
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief Numerical helper functions for various calculations.
 * @version 0.1
 * @date 2025-10-24
 *
 *
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

namespace uepm::numerical {

/**
 * @brief Generate linearly spaced points between start and end.
 *
 * @tparam T Numeric type (e.g., double, float).
 * @param start Starting value.
 * @param end Ending value.
 * @param num_points Number of points to generate.
 * @return std::vector<T> Vector of linearly spaced points.
 */
template < typename T>
std::vector<T> linspace(T start, T end, std::size_t num_points) {
    std::vector<T> points;
    if (num_points == 0) {
        return points;
    }
    if (num_points == 1) {
        points.push_back(start);
        return points;
    }
    T step = (end - start) / static_cast<T>(num_points - 1);
    for (std::size_t i = 0; i < num_points; ++i) {
        points.push_back(start + i * step);
    }
    return points;
}

/**
 * @brief Compute the differences between consecutive elements in a vector.
 * 
 * @tparam T 
 * @param input 
 * @return std::vector<T> 
 */
template <typename T>
std::vector<T> diff(const std::vector<T>& input) {
    std::vector<T> differences;
    if (input.size() < 2) {
        return differences;
    }
    differences.reserve(input.size() - 1);
    for (std::size_t i = 1; i < input.size(); ++i) {
        differences.push_back(input[i] - input[i - 1]);
    }
    return differences;
}

/**
 * @brief Returns the indices that would sort the input vector.
 *
 * @tparam T Numeric type of the input vector.
 * @param v Input vector to be sorted.
 * @return std::vector<std::size_t> Indices that would sort the vector.
 */
template <class T>
std::vector<std::size_t> argsort(const std::vector<T>& v) {
    std::vector<std::size_t> p(v.size());
    std::iota(p.begin(), p.end(), 0);
    std::stable_sort(p.begin(), p.end(), [&](std::size_t i, std::size_t j) { return v[i] < v[j]; });
    return p;
}

}  // namespace uepm::numerical