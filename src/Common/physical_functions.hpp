/**
 * @file physical_functions.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2025-10-12
 *
 *
 */

#pragma once

#include <cmath>
#include <numbers>

#include "physical_constants.hpp"

namespace uepm {

namespace physics {

inline double fermi_dirac_distribution(double energy_eV, double fermi_level_eV, double temperature_K) {
    const double kT = uepm::constants::k_b_eV * temperature_K;

    // Handle T <= 0 K as the T→0 limit (Heaviside step).
    if (!(kT > 0.0)) {
        if (energy_eV < fermi_level_eV) {
            return 1.0;
        }
        if (energy_eV > fermi_level_eV) {
            return 0.0;
        }
        return 0.5;  // at E = Ef, take the symmetric limit
    }

    const double x = (energy_eV - fermi_level_eV) / kT;

    // In double precision, |x| >= 40 puts f within ~1e-17 of 0 or 1.
    constexpr double X_CUTOFF = 40.0;
    if (x >= X_CUTOFF) {
        return 0.0;
    }
    if (x <= -X_CUTOFF) {
        return 1.0;
    }

    // Stable logistic: avoid large exp(x) when x > 0.
    if (x > 0.0) {
        const double emx = std::exp(-x);
        return emx / (1.0 + emx);
    } else {
        const double ex = std::exp(x);
        return 1.0 / (1.0 + ex);
    }
}

inline double d_de_fermi_dirac_dE(double energy_eV, double fermi_level_eV, double temperature_K) {
    const double kT = uepm::constants::k_b_eV * temperature_K;

    // T <= 0: derivative is a Dirac delta in theory; return 0 numerically.
    if (!(kT > 0.0)) {
        return 0.0;
    }

    const double x = (energy_eV - fermi_level_eV) / kT;

    constexpr double X_CUTOFF = 40.0;
    if (x >= X_CUTOFF || x <= -X_CUTOFF) {
        return 0.0;
    }

    // Compute f stably, then use df/dE = -(1/kT) * f * (1 - f).
    double f;
    if (x > 0.0) {
        const double emx = std::exp(-x);
        f                = emx / (1.0 + emx);
    } else {
        const double ex = std::exp(x);
        f               = 1.0 / (1.0 + ex);
    }
    return -(f * (1.0 - f)) / kT;
}

inline double bose_einstein_distribution(double energy_eV, double temperature_K) {
    // N0 = 1 / (exp(E / kT) - 1)
    const double x = energy_eV / (uepm::constants::k_b_eV * temperature_K);
    return 1.0 / std::expm1(x);  // stable for small x
}

}  // namespace PhysicalFunctions

}  // namespace uepm