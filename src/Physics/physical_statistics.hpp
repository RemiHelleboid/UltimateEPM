/**
 * @file physiscal_statistics.hpp
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
#include <vector>

#include "physical_constants.hpp"

namespace uepm {

namespace physic {

inline double fermi_dirac_statistics(double energy, double total_chemical_energy, double temperature) {
    double exponent = (energy - total_chemical_energy) / (constant::boltzmann_constant_SI * temperature);
    double n_i      = 1.0 / (1 + exp(exponent));
    return n_i;
}

inline double boson_einstein_statistics(double frequency, double temperature) {
    double exponent = (constant::planck_constant * frequency) / (constant::boltzmann_constant_SI * temperature);
    double n_i      = 1.0 / (exp(exponent) - 1);
    return n_i;
}




}  // namespace physic

}  // namespace uepm