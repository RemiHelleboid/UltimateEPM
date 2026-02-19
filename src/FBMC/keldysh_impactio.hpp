/**
 * @file keldysh_impactio.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief Impact ionization rates using Keldysh formula.
 * @version 0.1
 * @date 2026-02-17
 * 
 * 
 */

#pragma once

#include <cmath>

namespace uepm::fbmc {

/**
 * @brief Keldysh impact ionization model.
 * Filled with typical parameters for Si [Kamakura, 1994], but can be adapted for other materials.
 * 
 */
struct KeldyshImpactIonization {
    double m_P0 = 1.0e11;  // Pre-exponential factor (1/s)
    double m_alpha = 4.6;  // Exponent
    double m_E_threshold = 1.1;  // Threshold energy for impact ionization (eV)

    KeldyshImpactIonization() = default;
    KeldyshImpactIonization(double P0, double alpha, double E_threshold) : m_P0(P0), m_alpha(alpha), m_E_threshold(E_threshold) {}

    double compute_rate(double energy_eV) const {
        if (energy_eV < m_E_threshold) {
            return 0.0;  // No ionization below threshold
        }
        double excess_energy = energy_eV - m_E_threshold;
        return m_P0 * std::pow(excess_energy, m_alpha);
    }
};

}  // namespace uepm::fbmc