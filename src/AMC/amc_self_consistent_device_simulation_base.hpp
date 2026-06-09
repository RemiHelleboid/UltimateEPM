/**
 * @file amc_self_consistent_device_simulation_base.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-09
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "amc_avalanche_detector.hpp"
#include "amc_quench_circuit.hpp"
#include "amc_quench_detector.hpp"
#include "device.hpp"
#include "device_amc_simulation.hpp"
#include "materials.hpp"
#include "vtkWriter.hpp"

namespace uepm::amc {

enum class quench_biased_contact { anode, cathode };

struct options_self_consistent_device_amc_common {
    bool        m_frozen_field_mode = false;
    std::size_t m_poisson_frequency = 10;

    double m_anode_voltage   = 0.0;
    double m_cathode_voltage = 0.0;

    double m_contact_injection_particle_weight = 2.0;

    bool   m_initialize_particles_from_doping = true;
    double m_initial_particle_weight          = 2.0;

    passive_quench_circuit_options m_passive_quench_circuit{};

    quench_biased_contact m_quench_biased_contact = quench_biased_contact::cathode;

    // Converts signed Ramo current into current drawn from the biased circuit node.
    double m_ramo_current_to_quench_current_sign   = 1.0;
    double m_background_ramo_current_A             = 0.0;
    double m_avalanche_voltage_drop_threshold_V = 1.0;
    double m_quench_high_field_threshold_V_per_cm = 1.0e5;
    double m_quench_quiet_time_s                   = 1.0e-11;

    void validate() const;
};

class self_consistent_device_amc_simulation_base : public device_amc_simulation {
 protected:
    options_self_consistent_device_amc_common m_common_options;
    passive_quench_circuit                    m_quench_circuit;
    voltage_drop_avalanche_detector           m_avalanche_detector;
    successful_quench_detector                m_successful_quench_detector;

    self_consistent_device_amc_simulation_base(const device::device&                            simulation_device,
                                               const options_device_amc&                        simulation_options,
                                               const options_self_consistent_device_amc_common& common_options,
                                               int seed_random_generator = 0);

    self_consistent_device_amc_simulation_base(const device::device&                            simulation_device,
                                               const options_device_amc&                        simulation_options,
                                               const options_self_consistent_device_amc_common& common_options,
                                               const mesh::vector3&                             starting_position,
                                               std::size_t                                      number_electrons_start,
                                               std::size_t                                      number_holes_start,
                                               int seed_random_generator = 0);

    void validate_common_self_consistent_options() const;

    const options_self_consistent_device_amc_common& common_options() const;

    std::size_t poisson_frequency() const;

    double anode_voltage_for_poisson() const;
    double cathode_voltage_for_poisson() const;
    double device_bias_voltage_for_history() const;

    void advance_quench_circuit(double averaged_ramo_current_A, double dt_s, double sample_time_s);
    void update_successful_quench_detection(double sample_time_s, std::size_t impact_events_before_step);
    double max_particle_electric_field_V_per_cm() const;

    bool   quench_enabled() const;
    double quench_supply_voltage_for_history() const;
    double quench_node_voltage_for_history() const;
    double quench_device_current_for_history() const;
    double quench_resistor_current_for_history() const;
    double quench_voltage_drop_for_history() const;

 public:
    [[nodiscard]] const avalanche_detection_state& avalanche_detection() const;
    [[nodiscard]] const successful_quench_detection_state& successful_quench_detection() const;
};

}  // namespace uepm::amc
