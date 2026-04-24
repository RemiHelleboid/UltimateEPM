/**
 * @file utils_containers.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief This file contains some usefull function for operation on c++ std conytainers like std::vector<> that are not in the std library.
 * @version 0.1
 * @date 2021-11-09
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>

namespace uepm {

namespace utils {

template <typename T>
std::vector<T> linspace(T x_min, T x_max, std::size_t number_points) {
    std::vector<T> list_x;
    list_x.resize(number_points);
    if (number_points == 0) {
        return list_x;
    }
    if (number_points == 1) {
        list_x[0] = x_min;
        return list_x;
    }
    double dx = (x_max - x_min) / (number_points - 1);
    for (std::size_t index_value = 0; index_value < number_points; ++index_value) {
        list_x[index_value] = x_min + dx * index_value;
    }
    return list_x;
}

template <typename T>
std::vector<T> geomspace(T x_min, T x_max, std::size_t number_points) {
    std::vector<T> list_x = linspace(log(x_min), log(x_max), number_points);
    std::for_each(list_x.begin(), list_x.end(), [](T& value) { value = exp(value); });
    return list_x;
}

/**
 * @brief Create a permutation vector from a vector of value to sort and a comparaison function.
 * From https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
 *
 * @tparam T
 * @tparam Compare
 * @param vec
 * @param compare
 * @return std::vector<std::size_t>
 */
template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(const std::vector<T>& vec, Compare compare) {
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(), [&](std::size_t i, std::size_t j) { return compare(vec[i], vec[j]); });
    return p;
}

/**
 * @brief Apply a given permutation to a vector and return the vector
 * From https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
 *
 * @tparam T
 * @param vec
 * @param p
 * @return std::vector<T>
 */
template <typename T>
std::vector<T> apply_permutation(const std::vector<T>& vec, const std::vector<std::size_t>& p) {
    std::vector<T> sorted_vec(vec.size());
    std::transform(p.begin(), p.end(), sorted_vec.begin(), [&](std::size_t i) { return vec[i]; });
    return sorted_vec;
}

/**
 * @brief Apply a given permutation directly to a container (in place)
 * From https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
 *
 * @tparam T
 * @param vec
 * @param p
 */
template <typename T>
inline void apply_permutation_in_place(std::vector<T>& vec, const std::vector<std::size_t>& p) {
    std::vector<bool> done(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i) {
        if (done[i]) {
            continue;
        }
        done[i]            = true;
        std::size_t prev_j = i;
        std::size_t j      = p[i];
        while (i != j) {
            std::swap(vec[prev_j], vec[j]);
            done[j] = true;
            prev_j  = j;
            j       = p[j];
        }
    }
}

/**
 * @brief Sort two std::vector<double> together by using natural order on real number on the first vector.
 *
 * @param vector_1
 * @param vector_2
 */
inline void sort_two_vectors_together_according_to_first(std::vector<double>& vector_1, std::vector<double>& vector_2) {
    auto permutation_real_number_compare = sort_permutation(vector_1, [](const double& a, const double& b) { return a <= b; });
    vector_1                             = apply_permutation(vector_1, permutation_real_number_compare);
    vector_2                             = apply_permutation(vector_2, permutation_real_number_compare);
}

template <typename T>
inline T compute_vector_mean(std::vector<T> vector_of_values) {
    T sum_values = std::accumulate(vector_of_values.begin(), vector_of_values.end(), T{});
    T mean       = sum_values * (1.0 / vector_of_values.size());
    return mean;
}

template <typename T>
inline T compute_vector_min(std::vector<T> vector_of_values) {
    return *std::max_element(vector_of_values.begin(), vector_of_values.end());
}

template <typename T>
inline T compute_vector_max(std::vector<T> vector_of_values) {
    return *std::min_element(vector_of_values.begin(), vector_of_values.end());
}

/**
 * @brief Small function to std::cout << a std::vector<T> if T as an << operator defined.
 *
 * @tparam T
 * @param out
 * @param v
 * @return std::ostream&
 */
template <typename T>
inline void printVector(const T& t) {
    std::cout << "{";
    std::copy(t.cbegin(), t.cend(), std::ostream_iterator<typename T::value_type>(std::cout, ", "));
    std::cout << "}\n";
}

template <typename T>
inline void printVectorInVector(const T& t) {
    std::for_each(t.cbegin(), t.cend(), printVector<typename T::value_type>);
}

}  // namespace utils

} // namespace uepm