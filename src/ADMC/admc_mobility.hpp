/**
 * @file admc_mobility.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-15
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <string>

#include "admc_carrier.hpp"
#include "materials.hpp"

namespace uepm::ADMC {

struct arora_mobility_parameters {
    double mu_min_300_cm2_per_V_s;
    double mu_min_temperature_exponent;
    double mu_dop_300_cm2_per_V_s;
    double mu_dop_temperature_exponent;
    double reference_doping_300_cm_3;
    double reference_doping_temperature_exponent;
    double doping_exponent_300;
    double doping_exponent_temperature_exponent;
};

struct canali_high_field_parameters {
    double saturation_velocity_300_cm_per_s;
    double saturation_velocity_temperature_exponent;
    double beta_300;
    double beta_temperature_exponent;
};

struct carrier_mobility_parameters {
    arora_mobility_parameters    low_field;
    canali_high_field_parameters high_field;
};

struct silicon_mobility_parameters {
    double                      reference_temperature_K = 300.0;
    carrier_mobility_parameters electron;
    carrier_mobility_parameters hole;

    void validate() const;
};

silicon_mobility_parameters load_silicon_mobility_parameters(const uepm::physics::material_repository& repository,
                                                             const std::string& parameter_set = "arora-canali");
silicon_mobility_parameters load_silicon_mobility_parameters(const std::string& parameter_set = "arora-canali");

class silicon_arora_canali_mobility {
 public:
    silicon_arora_canali_mobility();
    explicit silicon_arora_canali_mobility(silicon_mobility_parameters parameters);

    double low_field_mobility_m2_per_V_s(carrier_type type,
                                         double       temperature_K,
                                         double       doping_concentration_cm_3) const;
    double mobility_m2_per_V_s(carrier_type type,
                               double       temperature_K,
                               double       doping_concentration_cm_3,
                               double       electric_field_V_per_m) const;

    const silicon_mobility_parameters& parameters() const noexcept { return m_parameters; }

 private:
    const carrier_mobility_parameters& parameters_for(carrier_type type) const;

    silicon_mobility_parameters m_parameters;
};

}  // namespace uepm::ADMC
