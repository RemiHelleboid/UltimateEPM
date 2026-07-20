/**
 * @file pbmc_self_consistent_device_simulation_3d.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-26
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "device.hpp"
#include "device_pbmc_simulation.hpp"
#include "materials.hpp"
#include "pbmc_quench_circuit.hpp"
#include "pbmc_self_consistent_device_simulation_base.hpp"
#include "poisson_solver_3d.hpp"
#include "vtkWriter.hpp"

namespace uepm::PBMC {

struct options_self_consistent_device_pbmc_3d {
    options_self_consistent_device_pbmc_common m_common{};
};

class self_consistent_device_pbmc_simulation_3d : public self_consistent_device_pbmc_simulation_base {
 private:
    options_self_consistent_device_pbmc_3d m_self_consistent_options;
    uepm::fem::poisson_solver_3d           m_poisson_solver;
    uepm::fem::EigenVector                 m_previous_poisson_solution;

    std::vector<std::size_t>                    m_list_element_contact;
    std::vector<std::shared_ptr<mesh::element>> m_list_element_contact_ptr;
    std::vector<std::size_t>                    m_list_element_contact_owner_index;
    std::vector<double>                         m_list_element_contact_equilibrium_charge;
    std::vector<mesh::vector3>                  m_list_element_contact_inward_direction;
    std::vector<std::vector<mesh::vector3>>     m_list_element_contact_face_vertices;

    std::minstd_rand m_contact_rng;

    void validate_self_consistent_options() const;
    void initialize_poisson_solver();
    void compute_unitary_potential();
    void initialize_particles_for_self_consistent_run();

    void place_initial_charges_according_to_doping(double particle_weight = 1.0);
    void initialize_contact_elements();
    void add_charges_at_contacts(std::size_t poisson_frequency);
    void add_missing_contact_charge_to_poisson_reservoir(std::size_t accumulation_steps);
    void update_self_consistent_potential();

 public:
    self_consistent_device_pbmc_simulation_3d(const device::device&                         simulation_device,
                                              const options_device_PBMC&                    simulation_options,
                                              const options_self_consistent_device_pbmc_3d& self_consistent_options,
                                              const physics::material_database&             material_database,
                                              int                                           seed_random_generator = 0);

    self_consistent_device_pbmc_simulation_3d(const device::device&                         simulation_device,
                                              const options_device_PBMC&                    simulation_options,
                                              const options_self_consistent_device_pbmc_3d& self_consistent_options,
                                              const physics::material_database&             material_database,
                                              const mesh::vector3&                          starting_position,
                                              std::size_t                                   number_electrons_start,
                                              std::size_t                                   number_holes_start,
                                              int                                           seed_random_generator = 0);

    void run_self_consistent_transport_simulation();

    void add_particle_charges_to_elements();
    void reset_element_charges();
    void recompute_vertex_space_charge_from_element_charges(std::size_t accumulation_steps);

    void export_current_state();
};

}  // namespace uepm::PBMC
