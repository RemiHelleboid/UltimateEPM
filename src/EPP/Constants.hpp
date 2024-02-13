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
constexpr double h_bar_eV       = 6.582119514e-16;
constexpr double m0             = 9.10938356e-31;
constexpr double q              = 1.602176634e-19;
constexpr double eps_zero       = 8.854187e-12;
constexpr double k_b            = 1.38064852e-23;
constexpr double k_b_eV         = 8.617333262145e-5;
constexpr double bohr_radius    = 5.291e-11;
constexpr double angstrom_to_m  = 1e-10;
constexpr double ryd_to_hartree = 1.0 / 2.0;
constexpr double hartree_to_eV  = 27.21138602;
constexpr double hartree_to_ryd = 2.0;
constexpr double hartree_to_J   = 4.3597447222071e-18;
constexpr double eV_to_J        = 1.602176634e-19;
constexpr double eV_to_ryd      = 1.0 / 13.6056980659;
constexpr double eV_to_hartree  = 1.0 / 27.21138602;
constexpr double eV_to_kg_m2    = 1.783e-36;
constexpr double eV_to_cm_1     = 8065.54429;
constexpr double eV_to_mol      = 9.64853399e4;
constexpr double pi             = M_PI;

}  // namespace Constants

}  // namespace EmpiricalPseudopotential