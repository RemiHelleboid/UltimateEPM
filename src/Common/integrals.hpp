/**
 * @file integrals.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief Numerical integration utilities.
 * @version 0.1
 * @date 2025-10-12
 *
 *
 */

#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <span>
#include <type_traits>
#include <vector>

namespace uepm::integrate {

template <typename T>
concept Real = std::is_floating_point_v<T>;

template <Real T>
inline bool is_strictly_increasing(std::span<const T> x) noexcept {
    for (std::size_t i = 1; i < x.size(); ++i) {
        if (!(x[i] > x[i - 1])) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Linear interpolation on a segment [x0,y0]-[x1,y1] at xq.
 *
 * @tparam T
 * @param x0
 * @param y0
 * @param x1
 * @param y1
 * @param xq
 * @return constexpr T
 */
template <Real T>
constexpr T lerp_on_segment(T x0, T y0, T x1, T y1, T xq) noexcept {
    const T dx = x1 - x0;
    if (dx == T(0)) {
        return y0;
    }
    const T t = (xq - x0) / dx;
    // return y0 + t * (y1 - y0);
    return std::fma(t, (y1 - y0), y0);
}

/**
 * @brief Segment index for value q in array x.
 * 
 * @tparam T 
 * @param x 
 * @param q 
 * @return std::size_t 
 */
template <Real T>
inline std::size_t segment_index(std::span<const T> x, T q) {
    const auto it = std::lower_bound(x.begin(), x.end(), q);
    if (it == x.begin()) {
        return 0;
    }
    if (it == x.end()) {
        return x.size() - 2;
    }
    std::size_t idx = static_cast<std::size_t>(it - x.begin());
    if (*it == q && idx > 0) {
        --idx;  // exact hit → take left segment
    }
    return std::min<std::size_t>(idx, x.size() - 2);
}

/**
 * @brief Compute the trapezoidal integral of a function given its samples.
 *
 * @tparam T
 * @param x
 * @param y
 * @return T
 */
template <Real T>
inline T trapz(std::span<const T> x, std::span<const T> y) {
    assert(x.size() == y.size());
    const std::size_t n = x.size();
    if (n < 2) {
        return T(0);
    }
    assert(is_strictly_increasing(x) && "x must be strictly increasing.");

    T s = T(0);
    for (std::size_t i = 1; i < n; ++i) {
        const T dx = x[i] - x[i - 1];
        s += dx * (y[i] + y[i - 1]) / T(2);
    }
    return s;
}

/**
 * @brief Compute the cumulative trapezoidal integral of a function given its samples.
 *
 * @tparam T
 * @param x
 * @param y
 * @param initial
 * @return std::vector<T>
 */
template <Real T>
inline std::vector<T> cumtrapz(std::span<const T> x, std::span<const T> y, T initial = T(0)) {
    assert(x.size() == y.size());
    const std::size_t n = x.size();
    std::vector<T>    out(n, initial);
    if (n < 2) {
        return out;
    }

    assert(is_strictly_increasing(x) && "x must be strictly increasing.");

    for (std::size_t i = 1; i < n; ++i) {
        const T dx = x[i] - x[i - 1];
        out[i]     = out[i - 1] + dx * (y[i] + y[i - 1]) / T(2);
    }
    return out;
}

/**
 * @brief Sa
 *
 * @tparam T
 * @param x
 * @param y
 * @param a
 * @param b
 * @return T
 */
template <Real T>
inline T integrate(std::span<const T> x, std::span<const T> y, T a, T b) {
    assert(x.size() == y.size());
    assert(a <= b);
    const std::size_t n = x.size();
    if (n < 2 || a == b) {
        return T(0);
    }

    assert(is_strictly_increasing(x) && "x must be strictly increasing.");

    // Clip to data support (no extrapolation surprises).
    const T lo = std::max(a, x.front());
    const T hi = std::min(b, x.back());
    if (hi <= lo) {
        return T(0);
    }

    const std::size_t il = segment_index(x, lo);
    const std::size_t ih = segment_index(x, hi);

    const T yl = lerp_on_segment(x[il], y[il], x[il + 1], y[il + 1], lo);
    const T yh = lerp_on_segment(x[ih], y[ih], x[ih + 1], y[ih + 1], hi);

    T s = T(0);

    if (il == ih) {
        // Both bounds inside one panel.
        return (hi - lo) * (yl + yh) / T(2);
    }

    // Left partial panel [lo, x[il+1]]
    s += (x[il + 1] - lo) * (yl + y[il + 1]) / T(2);

    // Full interior panels
    for (std::size_t i = il + 1; i < ih; ++i) {
        s += (x[i + 1] - x[i]) * (y[i + 1] + y[i]) / T(2);
    }

    // Right partial panel [x[ih], hi]
    s += (hi - x[ih]) * (y[ih] + yh) / T(2);

    return s;
}

template <Real T>
inline T trapz(const std::vector<T>& x, const std::vector<T>& y) {
    return trapz<T>(std::span<const T>(x), std::span<const T>(y));
}
template <Real T>
inline std::vector<T> cumtrapz(const std::vector<T>& x, const std::vector<T>& y, T initial = T(0)) {
    return cumtrapz<T>(std::span<const T>(x), std::span<const T>(y), initial);
}
template <Real T>
inline T integrate(const std::vector<T>& x, const std::vector<T>& y, T a, T b) {
    return integrate<T>(std::span<const T>(x), std::span<const T>(y), a, b);
}
template <Real T>
inline T simpson_uniform(const std::vector<T>& x, const std::vector<T>& y, T tol = T(1e-12)) {
    return simpson_uniform<T>(std::span<const T>(x), std::span<const T>(y), tol);
}

}  // namespace uepm::integrate
