/**
 * @file physics_function.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-10-21
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <functional>
#include <memory>
#include <optional>
#include <random>
#include <vector>

#include "physical_constants.hpp"

namespace uepm {

namespace physic {

inline double debye_length(double net_doping, double temperature, double absolute_permittivity) {
    return std::sqrt(absolute_permittivity * constant::boltzmann_constant_SI * temperature /
                     (constant::elementary_charge * constant::elementary_charge * net_doping));
}

inline double fermi_dirac_distribution(double energy, double temperature) {
    return 1.0 / (std::exp(energy / (constant::boltzmann_constant_SI * temperature)) + 1.0);
}   


}  // namespace physic

}  // namespace uepm