/**
 * @file bessel_func.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-08-10
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <cmath>
#include <iostream>

constexpr double pi = M_PI;

inline double J_neg_v(double v, double x) {
    return std::cos(-v * pi) * std::cyl_bessel_j(-v, x) - std::sin(-v * pi) * std::cyl_neumann(-v, x);
}

inline double generalized_bessel_cylindrical(double nu, double x) {
    if (nu >= 0) {
        return std::cyl_bessel_j(nu, x);
    } else {
        return J_neg_v(nu, x);
    }
}

inline double generalized_bessel(double nu, double x) {
    constexpr double epsilon = 1e-14;
    if (x < 0) {
        throw std::invalid_argument("Bessel function x argument must be positive. (x = " + std::to_string(x) + ")");
    }
    if (x < epsilon && nu == 0) {
        return 1.0;
    }
    if (x < epsilon && nu == 1) {
        return 0.0;
    }
    if (x < epsilon && nu >= 2) {
        return 0.0;
    }

    return sqrt(M_PI / (2.0 * x)) * generalized_bessel_cylindrical(nu + 0.5, x);
}