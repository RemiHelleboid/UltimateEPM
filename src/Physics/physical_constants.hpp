/**
 * @file physical_constants.h
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-05
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once

#include <cmath>

namespace uepm {

namespace physic {
namespace constant {

/**
 * @brief Vacuum absolute permittivity [F.m^-1]
 *
 */
constexpr double vacuum_permittivity = 8.8541878128e-12;

/**
 * @brief Elementary electronic charge [C]
 *
 */
constexpr double elementary_charge = 1.602176634e-19;

/**
 * @brief Boltzmann thermodynamic constant k_B [J.K^-1]
 *
 */
constexpr double boltzmann_constant_SI = 1.380649e-23;

/**
 * @brief Boltzmann thermodynamic constant k_B [eV.K^-1]
 *
 */
constexpr double boltzmann_constant_eV = boltzmann_constant_SI / elementary_charge;  // in eV.K^-1

/**
 * @brief Speed of light in vacuum [m.s^-1]
 *
 */
constexpr double speed_of_light = 299792458.0;

/**
 * @brief Planck constant [J.s]
 *
 */
constexpr double planck_constant = 6.62607015e-34;

/**
 * @brief Reduced Planck constant [J.s]
 *
 */
constexpr double reduced_planck_constant = planck_constant / (2.0 * M_PI);

}  //  namespace constant
}  // namespace physic
} // namespace uepm