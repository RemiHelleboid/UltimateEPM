/**
 * @file self_consistent_device_mmmc_simulation_2d.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-07-10
 * 
 * @copyright Copyright (c) 2026
 * 
 */

 
#pragma once

#include <cstddef>
#include <memory>
#include <random>
#include <vector>

#include "admc_transport.hpp"
#include "device_admc_simulation.hpp"
#include "mmmc_policy.hpp"
#include "pbmc_self_consistent_device_simulation_base.hpp"
#include "poisson_solver_2d.hpp"

namespace uepm::MMMC {

struct options_device_MMMC {
    PBMC::options_device_PBMC m_pbmc{};
    ADMC::options_device_ADMC m_admc{};

    void synchronize_from_pbmc();
    void validate() const;
};

struct options_self_consistent_device_MMMC_2d {
    PBMC::options_self_consistent_device_pbmc_common m_common{};
    bbox_transport_policy                            m_policy{};
    double                                           m_effective_depth_um = 1.0;

    void validate() const;
};

class self_consistent_device_mmmc_simulation_2d : public PBMC::self_consistent_device_pbmc_simulation_base {
 public:
    self_consistent_device_mmmc_simulation_2d(const device::device&                         simulation_device,
                                              const options_device_MMMC&                    simulation_options,
                                              const options_self_consistent_device_MMMC_2d& self_consistent_options,
                                              const physics::material_database&             material_database,
                                              int                                           seed_random_generator = 0);

    self_consistent_device_mmmc_simulation_2d(const device::device&                         simulation_device,
                                              const options_device_MMMC&                    simulation_options,
                                              const options_self_consistent_device_MMMC_2d& self_consistent_options,
                                              const physics::material_database&             material_database,
                                              const mesh::vector3&                          starting_position,
                                              std::size_t                                   number_electrons_start,
                                              std::size_t                                   number_holes_start,
                                              int                                           seed_random_generator = 0);

    void run_self_consistent_transport_simulation();

    [[nodiscard]] std::size_t get_total_number_electrons() const;
    [[nodiscard]] std::size_t get_total_number_holes() const;
    [[nodiscard]] std::size_t get_number_pbmc_particles() const noexcept;
    [[nodiscard]] std::size_t get_number_admc_particles() const noexcept;
    [[nodiscard]] std::size_t total_pbmc_to_admc_transfers() const noexcept;
    [[nodiscard]] std::size_t total_admc_to_pbmc_transfers() const noexcept;

    void export_current_mmmc_snapshot() const;
    void export_mmmc_particle_state_csv(const std::string& filename) const;

    void reset_element_charges();
    void add_particle_charges_to_elements();
    void recompute_vertex_space_charge_from_element_charges(std::size_t accumulation_steps);

 private:
    struct transfer_counters {
        std::size_t pbmc_to_admc = 0;
        std::size_t admc_to_pbmc = 0;
    };

    options_device_MMMC                     m_mmmc_options{};
    options_self_consistent_device_MMMC_2d  m_self_consistent_options{};
    ADMC::admc_transport_kernel             m_admc_transport{};
    std::vector<ADMC::device_admc_particle> m_admc_particles;
    fem::poisson_solver_2d                  m_poisson_solver;
    fem::EigenVector                        m_previous_poisson_solution;

    std::vector<std::size_t>                    m_list_element_contact;
    std::vector<std::shared_ptr<mesh::element>> m_list_element_contact_ptr;
    std::vector<std::size_t>                    m_list_element_contact_owner_index;
    std::vector<double>                         m_list_element_contact_equilibrium_charge;

    std::minstd_rand                 m_contact_rng;
    std::mt19937_64                  m_admc_rng;
    std::normal_distribution<double> m_standard_normal{0.0, 1.0};
    transfer_counters                m_last_transfer_counters{};
    transfer_counters                m_total_transfer_counters{};

    void validate_self_consistent_options() const;
    void initialize_contact_elements();
    void initialize_poisson_solver();
    void compute_unitary_potential();
    void initialize_particles_for_self_consistent_run();
    void place_initial_charges_according_to_doping(double particle_weight);
    void update_self_consistent_potential(bool publish_mesh_functions = true);

    double scale_integrated_2d_doping_to_carriers(double integrated_doping) const;
    double charge_deposition_factor(std::size_t accumulation_steps) const;
    double ramo_current_scale_factor() const override;
    double current_density_cell_volume_m3(const mesh::element& element) const override;

    void add_charges_at_contacts(std::size_t poisson_frequency);
    void add_missing_contact_charge_to_poisson_reservoir(std::size_t accumulation_steps);
    void apply_transport_policy();
    void advance_mmmc_particles_one_time_step();
    void advance_admc_particles_one_time_step(double dt_s);
    void update_admc_element_and_check_boundary(ADMC::device_admc_particle& particle);
    void remove_collected_admc_particles();
    void export_current_mmmc_particles_as_vtp(const std::string& directory) const;

    [[nodiscard]] ADMC::admc_local_environment local_admc_environment(const ADMC::device_admc_particle& particle) const;
    [[nodiscard]] ADMC::vector3                draw_standard_normal();
    [[nodiscard]] mesh::vector3                to_mesh_position_um(const ADMC::vector3& position_m) const;
    [[nodiscard]] ADMC::vector3                to_admc_position_m(const mesh::vector3& position_um) const;
    [[nodiscard]] ADMC::carrier_type           to_admc_type(PBMC::particle_type type) const;
    [[nodiscard]] PBMC::particle_type          to_pbmc_type(ADMC::carrier_type type) const;
    [[nodiscard]] PBMC::pbmc_transport_kernel& pbmc_transport_for(ADMC::carrier_type type);

    [[nodiscard]] std::pair<double, double> compute_admc_ramo_current() const;
    [[nodiscard]] double compute_admc_ramo_current_for_particle(const ADMC::device_admc_particle& particle) const;
    [[nodiscard]] double max_admc_particle_electric_field_V_per_cm() const;
};

}  // namespace uepm::MMMC
