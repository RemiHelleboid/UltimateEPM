/**
 * @file admc_transport.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-06-15
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#include "admc_transport.hpp"

#include <cmath>
#include <stdexcept>

#include "physical_constants.hpp"

namespace uepm::ADMC {

void admc_local_environment::validate() const {
    if (!std::isfinite(electric_field_V_per_m.x()) || !std::isfinite(electric_field_V_per_m.y()) ||
        !std::isfinite(electric_field_V_per_m.z())) {
        throw std::invalid_argument("electric field components must be finite");
    }
    if (!std::isfinite(doping_concentration_cm_3)) {
        throw std::invalid_argument("doping concentration must be finite");
    }
    if (!std::isfinite(lattice_temperature_K) || lattice_temperature_K <= 0.0) {
        throw std::invalid_argument("lattice temperature must be positive and finite");
    }
}

double einstein_diffusion_m2_per_s(double mobility_m2_per_V_s, double temperature_K) {
    if (!std::isfinite(mobility_m2_per_V_s) || mobility_m2_per_V_s < 0.0) {
        throw std::invalid_argument("mobility must be finite and non-negative");
    }
    if (!std::isfinite(temperature_K) || temperature_K <= 0.0) {
        throw std::invalid_argument("temperature must be positive and finite");
    }
    return mobility_m2_per_V_s * uepm::constants::k_B * temperature_K / uepm::constants::q_e;
}

void admc_transport_kernel::step(admc_particle&                 particle,
                                 const admc_local_environment& environment,
                                 double                         time_step_s,
                                 const vector3&                 standard_normal_draw) const {
    environment.validate();
    if (!std::isfinite(time_step_s) || time_step_s <= 0.0) {
        throw std::invalid_argument("ADMC time step must be positive and finite");
    }
    if (!std::isfinite(standard_normal_draw.x()) || !std::isfinite(standard_normal_draw.y()) ||
        !std::isfinite(standard_normal_draw.z())) {
        throw std::invalid_argument("normal draw components must be finite");
    }

    auto& state = particle.state();
    state.previous_position_m        = state.position_m;
    state.electric_field_V_per_m     = environment.electric_field_V_per_m;
    state.doping_concentration_cm_3  = environment.doping_concentration_cm_3;
    state.lattice_temperature_K      = environment.lattice_temperature_K;
    state.mobility_m2_per_V_s =
        m_mobility_model.mobility_m2_per_V_s(particle.type(),
                                             environment.lattice_temperature_K,
                                             environment.doping_concentration_cm_3,
                                             environment.electric_field_V_per_m.norm());
    state.diffusion_m2_per_s =
        einstein_diffusion_m2_per_s(state.mobility_m2_per_V_s, environment.lattice_temperature_K);
    state.drift_velocity_m_per_s =
        carrier_charge_sign(particle.type()) * state.mobility_m2_per_V_s * environment.electric_field_V_per_m;

    const double diffusion_sigma_m = std::sqrt(2.0 * state.diffusion_m2_per_s * time_step_s);
    const vector3 displacement_m =
        time_step_s * state.drift_velocity_m_per_s + diffusion_sigma_m * standard_normal_draw;
    state.position_m += displacement_m;
    state.total_velocity_m_per_s = displacement_m / time_step_s;
    state.time_s += time_step_s;
}

}  // namespace uepm::ADMC
