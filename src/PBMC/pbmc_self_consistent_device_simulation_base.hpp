/**
 * @file pbmc_self_consistent_device_simulation_base.hpp
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
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "device.hpp"
#include "device_pbmc_simulation.hpp"
#include "element.hpp"
#include "materials.hpp"
#include "pbmc_avalanche_detector.hpp"
#include "pbmc_quench_circuit.hpp"
#include "pbmc_quench_detector.hpp"
#include "nonlinear_poisson_solver.hpp"
#include "vtkWriter.hpp"

namespace uepm::PBMC {

struct scheduled_contact_voltage_event {
    double                        m_time_s = 0.0;
    std::map<std::string, double> m_contact_voltages_V;
};

struct options_self_consistent_device_pbmc_common {
    bool        m_frozen_field_mode = false;
    std::size_t m_poisson_frequency = 10;
    bool        m_nonlinear_steady_state_poisson = false;
    std::size_t m_nonlinear_poisson_warmup_steps = 10;
    uepm::fem::nonlinear_poisson_options m_nonlinear_poisson_options{};

    std::map<std::string, double>                m_contact_voltages_V;
    std::string                                  m_ramo_electrode;
    std::vector<scheduled_contact_voltage_event> m_contact_voltage_schedule;

    bool   m_enable_built_in_potential      = false;
    double m_built_in_contact_voltage_scale = 1.0;

    bool   m_enable_poisson_mixing                = false;
    double m_poisson_mixing_old_solution_fraction = 0.0;

    double m_contact_injection_particle_weight = 2.0;
    contact_injection_distribution m_contact_injection_distribution =
        contact_injection_distribution::velocity_weighted_maxwellian;

    bool        m_initialize_particles_from_doping = true;
    double      m_initial_particle_weight          = 2.0;
    std::string m_initial_particle_state_file;

    passive_quench_circuit_options m_passive_quench_circuit{};

    std::string m_quench_biased_contact = "cathode";

    // Converts signed Ramo current into current drawn from the biased circuit node.
    double m_ramo_current_to_quench_current_sign  = 1.0;
    double m_background_ramo_current_A            = 0.0;
    bool   m_auto_background_ramo_current         = false;
    double m_avalanche_voltage_drop_threshold_V   = 1.0;
    double m_quench_high_field_threshold_V_per_cm = 1.0e5;
    double m_quench_quiet_time_s                  = 1.0e-11;

    void validate() const;
};

double silicon_intrinsic_concentration_cm_3(double temperature_K);

class self_consistent_device_pbmc_simulation_base : public device_pbmc_simulation {
 protected:
    options_self_consistent_device_pbmc_common m_common_options;
    passive_quench_circuit                     m_quench_circuit;
    voltage_drop_avalanche_detector            m_avalanche_detector;
    successful_quench_detector                 m_successful_quench_detector;
    std::size_t                                m_next_contact_voltage_event_index = 0;

    std::map<std::string, double> m_built_in_contact_voltage_offsets_V;

    self_consistent_device_pbmc_simulation_base(const device::device&                             simulation_device,
                                                const options_device_PBMC&                        simulation_options,
                                                const options_self_consistent_device_pbmc_common& common_options,
                                                int seed_random_generator = 0);

    self_consistent_device_pbmc_simulation_base(const device::device&                             simulation_device,
                                                const options_device_PBMC&                        simulation_options,
                                                const options_self_consistent_device_pbmc_common& common_options,
                                                const mesh::vector3&                              starting_position,
                                                std::size_t number_electrons_start,
                                                std::size_t number_holes_start,
                                                int         seed_random_generator = 0);

    void validate_common_self_consistent_options() const;
    void update_built_in_contact_voltage_offset(const std::string&                                 contact_name,
                                                const std::vector<std::shared_ptr<mesh::element>>& contact_elements);
    bool is_transport_material_element(mesh::element& element);

    const options_self_consistent_device_pbmc_common& common_options() const;
    const std::map<std::string, double>&              contact_voltages_V() const;
    const std::string&                                ramo_electrode() const;

    std::size_t poisson_frequency() const;

    double contact_voltage_for_poisson(const std::string& contact_name) const;
    bool   apply_scheduled_contact_voltage_events(double time_s);
    double ramo_electrode_voltage_for_history() const;
    double reference_electrode_voltage_for_history() const;
    double device_bias_voltage_for_history() const;

    void   advance_quench_circuit(double averaged_ramo_current_A, double dt_s, double sample_time_s);
    void   update_successful_quench_detection(double sample_time_s, std::size_t impact_events_before_step);
    double max_particle_electric_field_V_per_cm() const;

    bool   quench_enabled() const;
    double quench_supply_voltage_for_history() const;
    double quench_node_voltage_for_history() const;
    double quench_device_current_for_history() const;
    double quench_resistor_current_for_history() const;
    double quench_voltage_drop_for_history() const;

 public:
    [[nodiscard]] const avalanche_detection_state&         avalanche_detection() const;
    [[nodiscard]] const successful_quench_detection_state& successful_quench_detection() const;
};

}  // namespace uepm::PBMC
