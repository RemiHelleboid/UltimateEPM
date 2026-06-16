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

#include <fmt/core.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "physical_constants.hpp"

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
    if (!std::isfinite(m_intrinsic_concentration_cm_3) || m_intrinsic_concentration_cm_3 <= 0.0) {
        throw std::invalid_argument("Intrinsic concentration must be positive and finite.");
    }
    if (!std::isfinite(m_built_in_contact_voltage_scale)) {
        throw std::invalid_argument("Built-in contact voltage scale must be finite.");
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

namespace {

struct contact_doping_summary {
    double donor_cm_3    = 0.0;
    double acceptor_cm_3 = 0.0;
};

contact_doping_summary summarize_contact_doping(
    const std::vector<std::shared_ptr<mesh::element>>& contact_elements) {
    double donor_integral    = 0.0;
    double acceptor_integral = 0.0;
    double measure_integral  = 0.0;

    for (const auto& element : contact_elements) {
        if (!element) {
            continue;
        }
        const double measure = std::abs(element->get_measure());
        if (measure <= 0.0) {
            continue;
        }
        const auto barycenter = element->get_barycenter();
        measure_integral += measure;
        donor_integral += measure * element->interpolate_scalar_at_location("DonorConcentration", barycenter);
        acceptor_integral += measure * element->interpolate_scalar_at_location("AcceptorConcentration", barycenter);
    }

    if (measure_integral <= 0.0) {
        throw std::runtime_error("Cannot estimate built-in potential from an empty contact element set.");
    }

    return {.donor_cm_3 = donor_integral / measure_integral, .acceptor_cm_3 = acceptor_integral / measure_integral};
}

double contact_equilibrium_voltage_offset_V(const contact_doping_summary& summary,
                                            double                         intrinsic_concentration_cm_3,
                                            double                         temperature_K) {
    const double donor_excess    = summary.donor_cm_3 - summary.acceptor_cm_3;
    const double acceptor_excess = summary.acceptor_cm_3 - summary.donor_cm_3;
    const double thermal_voltage = uepm::constants::k_B * temperature_K / uepm::constants::q_e;

    if (donor_excess > 0.0) {
        return thermal_voltage * std::log(donor_excess / intrinsic_concentration_cm_3);
    }
    if (acceptor_excess > 0.0) {
        return -thermal_voltage * std::log(acceptor_excess / intrinsic_concentration_cm_3);
    }
    return 0.0;
}

}  // namespace

void self_consistent_device_pbmc_simulation_base::update_built_in_contact_voltage_offsets(
    const std::vector<std::shared_ptr<mesh::element>>& anode_contact_elements,
    const std::vector<std::shared_ptr<mesh::element>>& cathode_contact_elements) {
    m_anode_built_in_voltage_offset_V   = 0.0;
    m_cathode_built_in_voltage_offset_V = 0.0;

    if (!m_common_options.m_enable_built_in_potential) {
        return;
    }

    const auto anode_doping   = summarize_contact_doping(anode_contact_elements);
    const auto cathode_doping = summarize_contact_doping(cathode_contact_elements);

    m_anode_built_in_voltage_offset_V =
        m_common_options.m_built_in_contact_voltage_scale *
        contact_equilibrium_voltage_offset_V(anode_doping,
                                             m_common_options.m_intrinsic_concentration_cm_3,
                                             m_simulation_options.m_lattice_temperature);
    m_cathode_built_in_voltage_offset_V =
        m_common_options.m_built_in_contact_voltage_scale *
        contact_equilibrium_voltage_offset_V(cathode_doping,
                                             m_common_options.m_intrinsic_concentration_cm_3,
                                             m_simulation_options.m_lattice_temperature);

    fmt::print("Built-in contact voltage offsets:\n");
    fmt::print("  intrinsic concentration: {:.6e} cm^-3\n", m_common_options.m_intrinsic_concentration_cm_3);
    fmt::print("  anode doping:   Nd={:.6e} cm^-3, Na={:.6e} cm^-3, offset={:.6e} V\n",
               anode_doping.donor_cm_3,
               anode_doping.acceptor_cm_3,
               m_anode_built_in_voltage_offset_V);
    fmt::print("  cathode doping: Nd={:.6e} cm^-3, Na={:.6e} cm^-3, offset={:.6e} V\n",
               cathode_doping.donor_cm_3,
               cathode_doping.acceptor_cm_3,
               m_cathode_built_in_voltage_offset_V);
    fmt::print("  built-in contact voltage difference cathode-anode: {:.6e} V\n",
               m_cathode_built_in_voltage_offset_V - m_anode_built_in_voltage_offset_V);
}

const options_self_consistent_device_pbmc_common& self_consistent_device_pbmc_simulation_base::common_options() const {
    return m_common_options;
}

std::size_t self_consistent_device_pbmc_simulation_base::poisson_frequency() const {
    return m_common_options.m_poisson_frequency;
}

double self_consistent_device_pbmc_simulation_base::anode_voltage_for_poisson() const {
    if (m_quench_circuit.is_enabled() && m_common_options.m_quench_biased_contact == quench_biased_contact::anode) {
        return m_quench_circuit.device_voltage_V() + m_anode_built_in_voltage_offset_V;
    }
    return m_common_options.m_anode_voltage + m_anode_built_in_voltage_offset_V;
}

double self_consistent_device_pbmc_simulation_base::cathode_voltage_for_poisson() const {
    if (m_quench_circuit.is_enabled() && m_common_options.m_quench_biased_contact == quench_biased_contact::cathode) {
        return m_quench_circuit.device_voltage_V() + m_cathode_built_in_voltage_offset_V;
    }
    return m_common_options.m_cathode_voltage + m_cathode_built_in_voltage_offset_V;
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
