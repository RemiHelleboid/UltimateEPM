/**
 * @file admc_transport.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-15
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <utility>

#include "admc_carrier.hpp"
#include "admc_mobility.hpp"

namespace uepm::ADMC {

struct admc_local_environment {
    vector3 electric_field_V_per_m{};
    double  doping_concentration_cm_3 = 0.0;
    double  lattice_temperature_K     = 300.0;

    void validate() const;
};

double einstein_diffusion_m2_per_s(double mobility_m2_per_V_s, double temperature_K);

class admc_transport_kernel {
 public:
    admc_transport_kernel() = default;
    explicit admc_transport_kernel(silicon_arora_canali_mobility mobility_model)
        : m_mobility_model(std::move(mobility_model)) {}

    void initialize_particle_state(admc_particle& particle, const admc_local_environment& environment) const;
    void step(admc_particle&                particle,
              const admc_local_environment& environment,
              double                        time_step_s,
              const vector3&                standard_normal_draw) const;

    const silicon_arora_canali_mobility& mobility_model() const noexcept { return m_mobility_model; }

 private:
    silicon_arora_canali_mobility m_mobility_model;
};

}  // namespace uepm::ADMC
