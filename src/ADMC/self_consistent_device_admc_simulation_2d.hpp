/**
 * @file self_consistent_device_admc_simulation_2d.hpp
 * @brief Self-consistent 2D ADMC device simulation.
 */

#pragma once

#include <cstddef>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "device_admc_simulation.hpp"
#include "materials.hpp"
#include "poisson_solver_2d.hpp"

namespace uepm::ADMC {

struct scheduled_contact_voltage_event {
    double m_time_s = 0.0;
    std::map<std::string, double> m_contact_voltages_V;
};

struct options_self_consistent_device_ADMC_common {
    std::size_t m_poisson_frequency = 10;

    std::map<std::string, double> m_contact_voltages_V;
    std::string                   m_ramo_electrode = "anode";
    std::vector<scheduled_contact_voltage_event> m_contact_voltage_schedule;

    bool   m_enable_built_in_potential      = false;
    double m_built_in_contact_voltage_scale = 1.0;

    bool   m_enable_poisson_mixing                = false;
    double m_poisson_mixing_old_solution_fraction = 0.0;

    bool   m_initialize_particles_from_doping = true;
    double m_initial_particle_weight          = 2.0;
    std::string m_initial_particle_state_file;
    double m_contact_injection_particle_weight = 2.0;

    void validate() const;
};

struct options_self_consistent_device_ADMC_2d {
    options_self_consistent_device_ADMC_common m_common{};

    double m_effective_depth_um = 1.0;

    void validate() const;
};

class self_consistent_device_admc_simulation_2d : public device_admc_simulation {
 public:
    self_consistent_device_admc_simulation_2d(const device::device& simulation_device,
                                              const options_device_ADMC& simulation_options,
                                              const options_self_consistent_device_ADMC_2d& self_consistent_options,
                                              const physics::material_database& material_database,
                                              const mesh::vector3& starting_position_um,
                                              std::size_t number_electrons_start,
                                              std::size_t number_holes_start,
                                              std::uint64_t random_seed = 5489u);

    void run_self_consistent_transport_simulation();

    void   reset_element_charges();
    void   add_particle_charges_to_elements();
    void   recompute_vertex_space_charge_from_element_charges(std::size_t accumulation_steps);
    double scale_integrated_2d_doping_to_carriers(double integrated_doping) const;

    const options_self_consistent_device_ADMC_2d& self_consistent_options() const noexcept {
        return m_self_consistent_options;
    }

 private:
    void validate_self_consistent_options() const;
    void initialize_contact_elements();
    void initialize_poisson_solver();
    void compute_unitary_potential();
    void update_self_consistent_potential(bool publish_mesh_functions = true);
    void initialize_particles_for_self_consistent_run();
    void update_built_in_contact_voltage_offset(
        const std::string& contact_name,
        const std::vector<std::shared_ptr<mesh::element>>& contact_elements);
    double contact_voltage_for_poisson(const std::string& contact_name) const;
    bool   apply_scheduled_contact_voltage_events(double time_s);

    void place_initial_charges_according_to_doping(double particle_weight);
    void add_charges_at_contacts(std::size_t poisson_frequency);
    void add_missing_contact_charge_to_poisson_reservoir(std::size_t accumulation_steps);
    options_self_consistent_device_ADMC_2d m_self_consistent_options;
    fem::poisson_solver_2d                m_poisson_solver;
    fem::EigenVector                      m_previous_poisson_solution;

    std::vector<std::size_t>                    m_list_element_contact;
    std::vector<std::shared_ptr<mesh::element>> m_list_element_contact_ptr;
    std::vector<double>                         m_list_element_contact_equilibrium_charge;
    std::map<std::string, double>               m_built_in_contact_voltage_offsets_V;
    std::size_t                                 m_next_contact_voltage_event_index = 0;

    std::minstd_rand m_contact_rng;
};

double silicon_intrinsic_concentration_admc_cm_3(double temperature_K);

}  // namespace uepm::ADMC
