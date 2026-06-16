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
#include <stdexcept>

#include "yaml-cpp/yaml.h"

namespace uepm::fbmc {

/**
 * @brief Keldysh impact ionization model.
 * Filled with typical parameters for Si [Kamakura, 1994], but can be adapted for other materials.
 *
 */
struct KeldyshImpactIonization {
    static constexpr double default_P0_s_1             = 1.0e11;
    static constexpr double default_alpha              = 4.6;
    static constexpr double default_energy_threshold_eV = 1.1;

    double m_P0          = default_P0_s_1;              // Pre-exponential factor (1/s)
    double m_alpha       = default_alpha;               // Exponent
    double m_E_threshold = default_energy_threshold_eV;  // Threshold energy for impact ionization (eV)

    KeldyshImpactIonization() = default;
    KeldyshImpactIonization(double P0, double alpha, double E_threshold)
        : m_P0(P0),
          m_alpha(alpha),
          m_E_threshold(E_threshold) {}

    void load_from_yaml(const YAML::Node& node) {
        m_P0          = node["P0"].as<double>();
        m_alpha       = node["alpha"].as<double>();
        m_E_threshold = node["energy_threshold"].as<double>();
        validate();
    }

    double compute_rate(double energy_eV) const {
        validate();
        if (!std::isfinite(energy_eV)) {
            throw std::invalid_argument("KeldyshImpactIonization::compute_rate: energy must be finite");
        }
        if (energy_eV < m_E_threshold) {
            return 0.0;  // No ionization below threshold
        }
        double excess_energy = energy_eV - m_E_threshold;
        return m_P0 * std::pow(excess_energy, m_alpha);
    }

    void validate() const {
        if (!(m_P0 >= 0.0) || !std::isfinite(m_P0)) {
            throw std::invalid_argument("KeldyshImpactIonization: P0 must be finite and non-negative");
        }
        if (!(m_alpha >= 0.0) || !std::isfinite(m_alpha)) {
            throw std::invalid_argument("KeldyshImpactIonization: alpha must be finite and non-negative");
        }
        if (!(m_E_threshold >= 0.0) || !std::isfinite(m_E_threshold)) {
            throw std::invalid_argument("KeldyshImpactIonization: threshold energy must be finite and non-negative");
        }
    }
};

}  // namespace uepm::fbmc
