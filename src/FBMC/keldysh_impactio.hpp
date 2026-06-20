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
#include <limits>
#include <stdexcept>
#include <string>

#include "materials.hpp"
#include "yaml-cpp/yaml.h"

namespace uepm::fbmc {

/**
 * @brief Keldysh impact ionization model.
 * Parameters are loaded from the selected material impact-ionization YAML profile.
 */
struct KeldyshImpactIonization {
    double m_P0          = std::numeric_limits<double>::quiet_NaN();  // Pre-exponential factor (1/s)
    double m_alpha       = std::numeric_limits<double>::quiet_NaN();  // Exponent
    double m_E_threshold = std::numeric_limits<double>::quiet_NaN();  // Threshold energy (eV)

    KeldyshImpactIonization() = default;
    KeldyshImpactIonization(double P0, double alpha, double E_threshold)
        : m_P0(P0),
          m_alpha(alpha),
          m_E_threshold(E_threshold) {}

    void load_from_yaml(const YAML::Node& node) {
        if (!node || !node.IsMap()) {
            throw std::invalid_argument("Keldysh impact-ionization configuration must be a map");
        }
        if (!node["P0"] || !node["alpha"] || !node["energy_threshold"]) {
            throw std::invalid_argument(
                "Keldysh impact-ionization configuration requires P0, alpha, and energy_threshold");
        }
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

inline KeldyshImpactIonization load_keldysh_impact_ionization(const uepm::physics::material_repository& repository,
                                                              const std::string&                        material_symbol,
                                                              const std::string&                        parameter_set) {
    const auto filename = repository.parameter_file(material_symbol, "impact_ionization", parameter_set);
    const auto config   = YAML::LoadFile(filename.string());
    if (!config || !config.IsMap() || !config["material"] || !config["model"] || !config["parameter_set"] ||
        config["material"].as<std::string>() != material_symbol ||
        config["model"].as<std::string>() != "impact_ionization" ||
        config["parameter_set"].as<std::string>() != parameter_set) {
        throw std::runtime_error("Invalid impact-ionization parameter file '" + filename.string() + "'.");
    }

    KeldyshImpactIonization result;
    try {
        result.load_from_yaml(config);
    } catch (const std::exception& error) {
        throw std::runtime_error("Invalid Keldysh impact-ionization configuration in '" + filename.string() +
                                 "': " + error.what());
    }
    return result;
}

}  // namespace uepm::fbmc
