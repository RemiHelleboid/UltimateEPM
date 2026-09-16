/**
 * @file mmmc_particle_transfer.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-07-10
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "mmmc_particle_transfer.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "unit_conversion.hpp"

namespace uepm::MMMC {
namespace {

ADMC::vector3 to_meters(const mesh::vector3& position_um) {
    return {position_um.x() * units::micron_to_meter,
            position_um.y() * units::micron_to_meter,
            position_um.z() * units::micron_to_meter};
}

mesh::vector3 to_microns(const ADMC::vector3& position_m) {
    return {position_m.x() * units::meter_to_micron,
            position_m.y() * units::meter_to_micron,
            position_m.z() * units::meter_to_micron};
}

void bias_pbmc_velocity_to_admc_drift(PBMC::pbmc_particle&         particle,
                                      const ADMC::vector3&         drift_velocity_m_per_s,
                                      PBMC::pbmc_transport_kernel& transport) {
    if (!std::isfinite(drift_velocity_m_per_s.x()) || !std::isfinite(drift_velocity_m_per_s.y()) ||
        !std::isfinite(drift_velocity_m_per_s.z())) {
        throw std::invalid_argument("ADMC drift velocity must be finite for MMMC transfer.");
    }

    const double drift_speed_m_per_s = drift_velocity_m_per_s.norm();
    if (!std::isfinite(drift_speed_m_per_s)) {
        throw std::invalid_argument("ADMC drift speed must be finite for MMMC transfer.");
    }
    if (drift_speed_m_per_s == 0.0) {
        return;
    }

    // ADMC drift velocity is an ensemble mean, whereas PBMC velocity is a
    // microscopic sample. Keep the thermally initialized PBMC state with
    // probability 1-p and align it with the drift with probability p. Since
    // the unbiased thermal sample has zero mean, p = |v_drift|/|v_aligned|
    // reproduces v_drift in expectation whenever the sampled microscopic
    // speed can represent it. Very low-energy samples saturate at p = 1.
    const PBMC::particle_state unbiased_state = particle.state();
    transport.set_particle_velocity_direction_preserving_energy(particle, drift_velocity_m_per_s);

    const double aligned_speed_m_per_s = particle.state().velocity.norm();
    if (!std::isfinite(aligned_speed_m_per_s) || aligned_speed_m_per_s <= 0.0) {
        particle.state() = unbiased_state;
        return;
    }

    const double alignment_probability = std::min(1.0, drift_speed_m_per_s / aligned_speed_m_per_s);
    if (transport.uniform01() >= alignment_probability) {
        particle.state() = unbiased_state;
    }
}

}  // namespace

ADMC::carrier_type to_admc_carrier_type(PBMC::particle_type type) {
    switch (type) {
        case PBMC::particle_type::electron:
            return ADMC::carrier_type::electron;
        case PBMC::particle_type::hole:
            return ADMC::carrier_type::hole;
    }
    throw std::invalid_argument("invalid PBMC particle type for MMMC transfer.");
}

PBMC::particle_type to_pbmc_particle_type(ADMC::carrier_type type) {
    switch (type) {
        case ADMC::carrier_type::electron:
            return PBMC::particle_type::electron;
        case ADMC::carrier_type::hole:
            return PBMC::particle_type::hole;
    }
    throw std::invalid_argument("invalid ADMC carrier type for MMMC transfer.");
}

ADMC::device_admc_particle convert_pbmc_to_admc(const PBMC::pbmc_particle& particle, std::size_t new_index) {
    const auto& pbmc_state = particle.state();

    ADMC::device_admc_particle result{
        .particle =
            ADMC::admc_particle(new_index, to_admc_carrier_type(particle.type()), to_meters(pbmc_state.position)),
        .containing_element = pbmc_state.m_containing_element,
        .weight             = particle.weight(),
        .crossed_contact    = pbmc_state.m_crossed_contact,
    };

    auto& admc_state                     = result.particle.state();
    admc_state.time_s                    = pbmc_state.time;
    admc_state.previous_position_m       = to_meters(pbmc_state.previous_position);
    admc_state.electric_field_V_per_m    = pbmc_state.electric_field * units::electric_field_V_per_cm_to_V_per_m;
    admc_state.drift_velocity_m_per_s    = pbmc_state.velocity;
    admc_state.total_velocity_m_per_s    = pbmc_state.velocity;
    admc_state.doping_concentration_cm_3 = pbmc_state.doping_concentration_cm_3;
    admc_state.lattice_temperature_K     = pbmc_state.lattice_temperature_K;

    return result;
}

PBMC::pbmc_particle convert_admc_to_pbmc(const ADMC::device_admc_particle& particle,
                                         std::size_t                       new_index,
                                         PBMC::pbmc_transport_kernel&      transport) {
    const auto& admc_state = particle.particle.state();

    PBMC::particle_state pbmc_state{};
    pbmc_state.time              = admc_state.time_s;
    pbmc_state.position          = to_microns(admc_state.position_m);
    pbmc_state.previous_position = to_microns(admc_state.previous_position_m);
    pbmc_state.velocity          = admc_state.drift_velocity_m_per_s;
    pbmc_state.electric_field    = admc_state.electric_field_V_per_m / units::electric_field_V_per_cm_to_V_per_m;
    pbmc_state.doping_concentration_cm_3   = admc_state.doping_concentration_cm_3;
    pbmc_state.impurity_concentration_cm_3 = std::abs(admc_state.doping_concentration_cm_3);
    pbmc_state.lattice_temperature_K       = admc_state.lattice_temperature_K;
    pbmc_state.m_containing_element        = particle.containing_element;
    pbmc_state.m_crossed_contact           = particle.crossed_contact;
    pbmc_state.valley_index                = 0;

    PBMC::pbmc_particle result(new_index, to_pbmc_particle_type(particle.particle.type()), pbmc_state, particle.weight);
    transport.initialize_particle_state(result, admc_state.lattice_temperature_K);
    bias_pbmc_velocity_to_admc_drift(result, admc_state.drift_velocity_m_per_s, transport);
    return result;
}

}  // namespace uepm::MMMC
