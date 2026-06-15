/**
 * @file admc_mobility.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-15
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "admc_mobility.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

#include "yaml-cpp/yaml.h"

namespace uepm::ADMC {
namespace {

constexpr double cm2_to_m2           = 1.0e-4;
constexpr double V_per_m_to_V_per_cm = 1.0e-2;

void require_positive_finite(double value, const char* name) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(std::string(name) + " must be positive and finite");
    }
}

void validate_carrier_parameters(const carrier_mobility_parameters& parameters) {
    require_positive_finite(parameters.low_field.mu_min_300_cm2_per_V_s, "minimum mobility");
    require_positive_finite(parameters.low_field.mu_dop_300_cm2_per_V_s, "doping mobility");
    require_positive_finite(parameters.low_field.reference_doping_300_cm_3, "reference doping");
    require_positive_finite(parameters.low_field.doping_exponent_300, "doping exponent");
    require_positive_finite(parameters.high_field.saturation_velocity_300_cm_per_s, "saturation velocity");
    require_positive_finite(parameters.high_field.beta_300, "Canali beta");
    const double exponents[] = {
        parameters.low_field.mu_min_temperature_exponent,
        parameters.low_field.mu_dop_temperature_exponent,
        parameters.low_field.reference_doping_temperature_exponent,
        parameters.low_field.doping_exponent_temperature_exponent,
        parameters.high_field.saturation_velocity_temperature_exponent,
        parameters.high_field.beta_temperature_exponent,
    };
    for (const double exponent : exponents) {
        if (!std::isfinite(exponent)) {
            throw std::invalid_argument("mobility temperature exponents must be finite");
        }
    }
}

double required_double(const YAML::Node& node, const char* key, const std::string& filename) {
    const auto value = node[key];
    if (!value) {
        throw std::runtime_error("ADMC mobility parameter file '" + filename + "' does not define '" + key + "'.");
    }
    return value.as<double>();
}

carrier_mobility_parameters parse_carrier_parameters(const YAML::Node&  node,
                                                     const std::string& carrier,
                                                     const std::string& filename) {
    if (!node || !node.IsMap()) {
        throw std::runtime_error("ADMC mobility parameter file '" + filename + "' does not define '" + carrier + "'.");
    }
    const auto low_field  = node["low_field"];
    const auto high_field = node["high_field"];
    if (!low_field || !low_field.IsMap() || !high_field || !high_field.IsMap()) {
        throw std::runtime_error("ADMC mobility parameter file '" + filename +
                                 "' must define low_field and high_field parameters for '" + carrier + "'.");
    }

    return {
        .low_field =
            {
                .mu_min_300_cm2_per_V_s      = required_double(low_field, "mu_min_300_cm2_per_V_s", filename),
                .mu_min_temperature_exponent = required_double(low_field, "mu_min_temperature_exponent", filename),
                .mu_dop_300_cm2_per_V_s      = required_double(low_field, "mu_dop_300_cm2_per_V_s", filename),
                .mu_dop_temperature_exponent = required_double(low_field, "mu_dop_temperature_exponent", filename),
                .reference_doping_300_cm_3   = required_double(low_field, "reference_doping_300_cm_3", filename),
                .reference_doping_temperature_exponent =
                    required_double(low_field, "reference_doping_temperature_exponent", filename),
                .doping_exponent_300 = required_double(low_field, "doping_exponent_300", filename),
                .doping_exponent_temperature_exponent =
                    required_double(low_field, "doping_exponent_temperature_exponent", filename),
            },
        .high_field =
            {
                .saturation_velocity_300_cm_per_s =
                    required_double(high_field, "saturation_velocity_300_cm_per_s", filename),
                .saturation_velocity_temperature_exponent =
                    required_double(high_field, "saturation_velocity_temperature_exponent", filename),
                .beta_300                  = required_double(high_field, "beta_300", filename),
                .beta_temperature_exponent = required_double(high_field, "beta_temperature_exponent", filename),
            },
    };
}

}  // namespace

void silicon_mobility_parameters::validate() const {
    require_positive_finite(reference_temperature_K, "reference temperature");
    validate_carrier_parameters(electron);
    validate_carrier_parameters(hole);
}

silicon_mobility_parameters load_silicon_mobility_parameters(const uepm::physics::material_repository& repository,
                                                             const std::string&                        parameter_set) {
    const auto material = repository.load_material("Si");
    if (material.id != uepm::physics::material_id::silicon) {
        throw std::runtime_error("ADMC silicon mobility parameters require material 'Si'.");
    }

    const auto filename = repository.parameter_file("Si", "admc", parameter_set);
    const auto config   = YAML::LoadFile(filename.string());
    if (!config.IsMap() || !config["schema_version"] || config["schema_version"].as<int>() != 1 ||
        !config["material"] || config["material"].as<std::string>() != "Si" || !config["model"] ||
        config["model"].as<std::string>() != "arora_canali_mobility" || !config["parameter_set"] ||
        config["parameter_set"].as<std::string>() != parameter_set) {
        throw std::runtime_error("Invalid ADMC mobility parameter file '" + filename.string() + "'.");
    }

    silicon_mobility_parameters parameters;
    parameters.reference_temperature_K = required_double(config, "reference_temperature_K", filename.string());
    parameters.electron                = parse_carrier_parameters(config["electron"], "electron", filename.string());
    parameters.hole                    = parse_carrier_parameters(config["hole"], "hole", filename.string());
    parameters.validate();
    return parameters;
}

silicon_mobility_parameters load_silicon_mobility_parameters(const std::string& parameter_set) {
    return load_silicon_mobility_parameters(uepm::physics::material_repository{}, parameter_set);
}

silicon_arora_canali_mobility::silicon_arora_canali_mobility()
    : silicon_arora_canali_mobility(load_silicon_mobility_parameters()) {}

silicon_arora_canali_mobility::silicon_arora_canali_mobility(silicon_mobility_parameters parameters)
    : m_parameters(std::move(parameters)) {
    m_parameters.validate();
}

const carrier_mobility_parameters& silicon_arora_canali_mobility::parameters_for(carrier_type type) const {
    switch (type) {
        case carrier_type::electron:
            return m_parameters.electron;
        case carrier_type::hole:
            return m_parameters.hole;
    }
    throw std::invalid_argument("invalid ADMC carrier type");
}

double silicon_arora_canali_mobility::low_field_mobility_m2_per_V_s(carrier_type type,
                                                                    double       temperature_K,
                                                                    double       doping_concentration_cm_3) const {
    require_positive_finite(temperature_K, "temperature");
    if (!std::isfinite(doping_concentration_cm_3)) {
        throw std::invalid_argument("doping concentration must be finite");
    }

    const auto&  p                    = parameters_for(type).low_field;
    const double relative_temperature = temperature_K / m_parameters.reference_temperature_K;
    const double mu_min_cm2_per_V_s =
        p.mu_min_300_cm2_per_V_s * std::pow(relative_temperature, p.mu_min_temperature_exponent);
    const double mu_dop_cm2_per_V_s =
        p.mu_dop_300_cm2_per_V_s * std::pow(relative_temperature, p.mu_dop_temperature_exponent);
    const double reference_doping_cm_3 =
        p.reference_doping_300_cm_3 * std::pow(relative_temperature, p.reference_doping_temperature_exponent);
    const double doping_exponent =
        p.doping_exponent_300 * std::pow(relative_temperature, p.doping_exponent_temperature_exponent);
    const double doping_ratio = std::abs(doping_concentration_cm_3) / reference_doping_cm_3;

    return cm2_to_m2 * (mu_min_cm2_per_V_s + mu_dop_cm2_per_V_s / (1.0 + std::pow(doping_ratio, doping_exponent)));
}

double silicon_arora_canali_mobility::mobility_m2_per_V_s(carrier_type type,
                                                          double       temperature_K,
                                                          double       doping_concentration_cm_3,
                                                          double       electric_field_V_per_m) const {
    if (!std::isfinite(electric_field_V_per_m) || electric_field_V_per_m < 0.0) {
        throw std::invalid_argument("electric field magnitude must be finite and non-negative");
    }

    const auto&  p                    = parameters_for(type).high_field;
    const double relative_temperature = temperature_K / m_parameters.reference_temperature_K;
    const double low_field_mobility   = low_field_mobility_m2_per_V_s(type, temperature_K, doping_concentration_cm_3);
    const double low_field_mobility_cm2_per_V_s = low_field_mobility / cm2_to_m2;
    const double saturation_velocity_cm_per_s =
        p.saturation_velocity_300_cm_per_s * std::pow(relative_temperature, p.saturation_velocity_temperature_exponent);
    const double beta            = p.beta_300 * std::pow(relative_temperature, p.beta_temperature_exponent);
    const double field_V_per_cm  = electric_field_V_per_m * V_per_m_to_V_per_cm;
    const double velocity_ratio  = low_field_mobility_cm2_per_V_s * field_V_per_cm / saturation_velocity_cm_per_s;
    const double field_reduction = std::pow(1.0 + std::pow(velocity_ratio, beta), 1.0 / beta);

    return low_field_mobility / field_reduction;
}

}  // namespace uepm::ADMC
