/**
 * @file pbmc_self_consistent_device_simulation_base.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-09
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "pbmc_self_consistent_device_simulation_base.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace uepm::PBMC {

void options_self_consistent_device_pbmc_common::validate() const {
    if (m_poisson_frequency == 0) {
        throw std::invalid_argument("Poisson frequency must be positive.");
    }
    if (!std::isfinite(m_anode_voltage)) {
        throw std::invalid_argument("Anode voltage must be finite.");
    }
    if (!std::isfinite(m_cathode_voltage)) {
        throw std::invalid_argument("Cathode voltage must be finite.");
    }
    if (!std::isfinite(m_contact_injection_particle_weight) || m_contact_injection_particle_weight <= 0.0) {
        throw std::invalid_argument("Contact injection particle weight must be positive.");
    }
    if (!std::isfinite(m_initial_particle_weight) || m_initial_particle_weight <= 0.0) {
        throw std::invalid_argument("Initial particle weight must be positive.");
    }
    if (!std::isfinite(m_ramo_current_to_quench_current_sign) || m_ramo_current_to_quench_current_sign == 0.0) {
        throw std::invalid_argument("Ramo current to quench current sign must be finite and non-zero.");
    }
    if (!std::isfinite(m_avalanche_voltage_drop_threshold_V) || m_avalanche_voltage_drop_threshold_V <= 0.0) {
        throw std::invalid_argument("Avalanche voltage-drop threshold must be positive and finite.");
    }
    if (!std::isfinite(m_quench_high_field_threshold_V_per_cm) || m_quench_high_field_threshold_V_per_cm <= 0.0) {
        throw std::invalid_argument("Quench high-field threshold must be positive and finite.");
    }
    if (!std::isfinite(m_quench_quiet_time_s) || m_quench_quiet_time_s <= 0.0) {
        throw std::invalid_argument("Quench quiet time must be positive and finite.");
    }
}

self_consistent_device_pbmc_simulation_base::self_consistent_device_pbmc_simulation_base(
    const device::device&                            simulation_device,
    const options_device_PBMC&                        simulation_options,
    const options_self_consistent_device_pbmc_common& common_options,
    int                                              seed_random_generator)
    : device_pbmc_simulation(simulation_device, simulation_options, seed_random_generator),
      m_common_options(common_options),
      m_quench_circuit(common_options.m_passive_quench_circuit),
      m_avalanche_detector(common_options.m_avalanche_voltage_drop_threshold_V),
      m_successful_quench_detector(common_options.m_quench_quiet_time_s) {
    validate_common_self_consistent_options();
}

self_consistent_device_pbmc_simulation_base::self_consistent_device_pbmc_simulation_base(
    const device::device&                            simulation_device,
    const options_device_PBMC&                        simulation_options,
    const options_self_consistent_device_pbmc_common& common_options,
    const mesh::vector3&                             starting_position,
    std::size_t                                      number_electrons_start,
    std::size_t                                      number_holes_start,
    int                                              seed_random_generator)
    : device_pbmc_simulation(simulation_device,
                            simulation_options,
                            starting_position,
                            number_electrons_start,
                            number_holes_start,
                            seed_random_generator),
      m_common_options(common_options),
      m_quench_circuit(common_options.m_passive_quench_circuit),
      m_avalanche_detector(common_options.m_avalanche_voltage_drop_threshold_V),
      m_successful_quench_detector(common_options.m_quench_quiet_time_s) {
    validate_common_self_consistent_options();
}

void self_consistent_device_pbmc_simulation_base::validate_common_self_consistent_options() const {
    m_common_options.validate();
}

const options_self_consistent_device_pbmc_common& self_consistent_device_pbmc_simulation_base::common_options() const {
    return m_common_options;
}

std::size_t self_consistent_device_pbmc_simulation_base::poisson_frequency() const {
    return m_common_options.m_poisson_frequency;
}

double self_consistent_device_pbmc_simulation_base::anode_voltage_for_poisson() const {
    if (m_quench_circuit.is_enabled() && m_common_options.m_quench_biased_contact == quench_biased_contact::anode) {
        return m_quench_circuit.device_voltage_V();
    }
    return m_common_options.m_anode_voltage;
}

double self_consistent_device_pbmc_simulation_base::cathode_voltage_for_poisson() const {
    if (m_quench_circuit.is_enabled() && m_common_options.m_quench_biased_contact == quench_biased_contact::cathode) {
        return m_quench_circuit.device_voltage_V();
    }
    return m_common_options.m_cathode_voltage;
}

double self_consistent_device_pbmc_simulation_base::device_bias_voltage_for_history() const {
    return cathode_voltage_for_poisson() - anode_voltage_for_poisson();
}

void self_consistent_device_pbmc_simulation_base::advance_quench_circuit(double averaged_ramo_current_A,
                                                                        double dt_s,
                                                                        double sample_time_s) {
    if (!m_quench_circuit.is_enabled()) {
        return;
    }
    const double circuit_current_A = m_common_options.m_ramo_current_to_quench_current_sign * averaged_ramo_current_A;
    m_quench_circuit.advance(circuit_current_A, dt_s);
    m_avalanche_detector.update(true, sample_time_s, m_quench_circuit.voltage_drop_V());
}

void self_consistent_device_pbmc_simulation_base::update_successful_quench_detection(
    double      sample_time_s,
    std::size_t impact_events_before_step) {
    const bool had_impact_ionization_event =
        m_simulation_history.m_impact_ionization_positions.size() > impact_events_before_step;
    const bool has_high_field_particle =
        max_particle_electric_field_V_per_cm() >= m_common_options.m_quench_high_field_threshold_V_per_cm;
    m_successful_quench_detector.update(m_avalanche_detector.state().m_detected,
                                        sample_time_s,
                                        has_high_field_particle,
                                        had_impact_ionization_event);
}

double self_consistent_device_pbmc_simulation_base::max_particle_electric_field_V_per_cm() const {
    double max_field_V_per_cm = 0.0;
    for (const auto& particle : m_list_particles) {
        max_field_V_per_cm = std::max(max_field_V_per_cm, particle->state().electric_field.norm());
    }
    return max_field_V_per_cm;
}

bool self_consistent_device_pbmc_simulation_base::quench_enabled() const { return m_quench_circuit.is_enabled(); }

double self_consistent_device_pbmc_simulation_base::quench_supply_voltage_for_history() const {
    if (!m_quench_circuit.is_enabled()) {
        return 0.0;
    }
    return m_quench_circuit.options().m_bias_voltage_V;
}

double self_consistent_device_pbmc_simulation_base::quench_node_voltage_for_history() const {
    if (!m_quench_circuit.is_enabled()) {
        return 0.0;
    }
    return m_quench_circuit.device_voltage_V();
}

double self_consistent_device_pbmc_simulation_base::quench_device_current_for_history() const {
    if (!m_quench_circuit.is_enabled()) {
        return 0.0;
    }
    return m_quench_circuit.device_current_A();
}

double self_consistent_device_pbmc_simulation_base::quench_resistor_current_for_history() const {
    if (!m_quench_circuit.is_enabled()) {
        return 0.0;
    }
    return m_quench_circuit.resistor_current_A();
}

double self_consistent_device_pbmc_simulation_base::quench_voltage_drop_for_history() const {
    if (!m_quench_circuit.is_enabled()) {
        return 0.0;
    }
    return m_quench_circuit.voltage_drop_V();
}

const avalanche_detection_state& self_consistent_device_pbmc_simulation_base::avalanche_detection() const {
    return m_avalanche_detector.state();
}

const successful_quench_detection_state& self_consistent_device_pbmc_simulation_base::successful_quench_detection()
    const {
    return m_successful_quench_detector.state();
}

}  // namespace uepm::PBMC
