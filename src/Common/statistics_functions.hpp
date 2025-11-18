/**
 * @file statistics_functions.hpp
 * @brief Super-simple statistical helpers (sum, mean, variance, stdev, min, max, median)
 * @date 2025-11-06
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace uepm::statistics {

// Sum
inline double sum(const std::vector<double>& v) { return std::accumulate(v.begin(), v.end(), 0.0); }

// Mean
inline double mean(const std::vector<double>& v) {
    if (v.empty()) {
        throw std::invalid_argument("mean: empty vector");
    }
    return sum(v) / static_cast<double>(v.size());
}

// Mean
inline double weighted_mean(const std::vector<double>& v, const std::vector<double>& weights) {
    if (v.empty() || weights.empty() || v.size() != weights.size()) {
        throw std::invalid_argument("weighted_mean: invalid vectors");
    }
    double weighted_sum = 0.0;
    double total_weight = 0.0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        weighted_sum += v[i] * weights[i];
        total_weight += weights[i];
    }
    return weighted_sum / total_weight;
}

// Variance (sample=true -> Bessel corrected)
inline double variance(const std::vector<double>& v, bool sample = true) {
    const std::size_t n = v.size();
    if (n == 0) {
        throw std::invalid_argument("variance: empty vector");
    }
    if (sample && n < 2) {
        throw std::invalid_argument("variance: need at least 2 values for sample variance");
    }

    const double m   = mean(v);
    double       acc = 0.0;
    for (double x : v) {
        const double d = x - m;
        acc += d * d;
    }
    return acc / static_cast<double>(sample ? (n - 1) : n);
}

// Standard deviation
inline double stdev(const std::vector<double>& v, bool sample = true) { return std::sqrt(variance(v, sample)); }

// Min / Max
inline double min(const std::vector<double>& v) {
    if (v.empty()) {
        throw std::invalid_argument("min: empty vector");
    }
    return *std::min_element(v.begin(), v.end());
}

inline double max(const std::vector<double>& v) {
    if (v.empty()) {
        throw std::invalid_argument("max: empty vector");
    }
    return *std::max_element(v.begin(), v.end());
}

}  // namespace uepm::statistics
