/**
 * @file valley_model.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-04-12
 * 
 * 
 */

 
#include "valley_model.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace uepm::amc {

valley_model::vector3 valley_model::draw_random_k_valley_at_energy(double energy_eV, std::mt19937_64& rng) const {
    if (energy_eV < 0.0) {
        throw std::invalid_argument("energy must be non-negative");
    }

    std::uniform_real_distribution<double> unif01(0.0, 1.0);

    const double uz  = 2.0 * unif01(rng) - 1.0;
    const double phi = 2.0 * uepm::constants::pi * unif01(rng);
    const double rxy = std::sqrt(std::max(0.0, 1.0 - uz * uz));

    const vector3 direction{rxy * std::cos(phi), rxy * std::sin(phi), uz};

    return k_valley_from_energy_direction(energy_eV, direction);
}

}  // namespace uepm::amc