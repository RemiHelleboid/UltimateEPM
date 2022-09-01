/**
 * @file Constants.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-07-26
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <cmath>

namespace EmpiricalPseudopotential {

namespace Constants {

constexpr double Ryd_to_eV      = 13.6;
constexpr double h              = 6.62606957e-34;
constexpr double h_bar          = 6.62606957e-34 / (2 * M_PI);
constexpr double m0             = 9.10938356e-31;
constexpr double q              = 1.6e-19;
constexpr double eps_zero       = 8.854187e-12;
constexpr double k_b            = 1.38064852e-23;
constexpr double ab             = 5.291e-11;
constexpr double angstrom_to_m  = 1e-10;
constexpr double ryd_to_hartree = 1.0 / 2.0;

}  // namespace Constants

}  // namespace EmpiricalPseudopotential