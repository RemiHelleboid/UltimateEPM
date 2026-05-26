/**
 * @file amc_self_consistent_device_simulation_3d.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-26
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "amc_self_consistent_device_simulation_3d.hpp"

#include <fmt/core.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace uepm::amc {

void self_consistent_device_amc_simulation_3d::validate_self_consistent_options() const {
    if (m_dimension != 3) {
        throw std::invalid_argument("self_consistent_device_amc_simulation_3d requires a 3D device.");
    }

    if (m_self_consistent_options.m_poisson_frequency == 0) {
        throw std::invalid_argument("Poisson frequency must be positive.");
    }
}

/**
 * @brief We need to identify the elements adjacent to the contacts, and compute their equilibrium charge (integral of doping concentration)
 * to be able to properly update their charge during the simulation (so that they are not fixed at the equilibrium charge, but can exchange
 * charge with the device while still having a reference equilibrium charge for the poisson solver).
 *
 */
void self_consistent_device_amc_simulation_3d::initialize_contact_elements() {
    const auto list_contact_elem_anode = m_device.get_p_mesh()->get_idx_bulk_elements_adjacent_to_contact_region("anode");

    const auto list_contact_elem_cathode = m_device.get_p_mesh()->get_idx_bulk_elements_adjacent_to_contact_region("cathode");

    const auto list_bulk_elements = m_device.get_p_mesh()->get_list_bulk_element();

    m_list_element_contact.clear();
    m_list_element_contact_ptr.clear();
    m_list_element_contact_equilibrium_charge.clear();

    m_list_element_contact.reserve(list_contact_elem_anode.size() + list_contact_elem_cathode.size());
    m_list_element_contact_ptr.reserve(list_contact_elem_anode.size() + list_contact_elem_cathode.size());
    m_list_element_contact_equilibrium_charge.reserve(list_contact_elem_anode.size() + list_contact_elem_cathode.size());

    const auto add_contact_element = [&](std::size_t element_index) {
        if (element_index >= list_bulk_elements.size()) {
            throw std::runtime_error("Contact-adjacent bulk element index is out of range.");
        }

        auto element = list_bulk_elements[element_index];

        m_list_element_contact.push_back(element_index);
        m_list_element_contact_ptr.push_back(element);

        const double equilibrium_charge = element->integrate_scalar("DopingConcentration");

        m_list_element_contact_equilibrium_charge.push_back(equilibrium_charge);
    };

    for (const auto element_index : list_contact_elem_anode) {
        add_contact_element(element_index);
    }

    for (const auto element_index : list_contact_elem_cathode) {
        add_contact_element(element_index);
    }
}

/**
 * @brief We compute the charge to add at each contact-adjacent element based on the difference between the equilibrium charge (computed
 * from the doping concentration) and the current accumulated charge in the element, and we draw particles accordingly to add this charge at
 * the contacts.
 *
 * @param poisson_frequency
 */
void self_consistent_device_amc_simulation_3d::add_charges_at_contacts(std::size_t poisson_frequency) {
    if (poisson_frequency == 0) {
        throw std::invalid_argument("Poisson frequency must be positive.");
    }

    const std::size_t number_contact_elements = m_list_element_contact_ptr.size();

    std::vector<mesh::vector3> electron_positions;
    std::vector<mesh::vector3> hole_positions;

    std::vector<double> electron_charge_to_add(number_contact_elements, 0.0);
    std::vector<double> hole_charge_to_add(number_contact_elements, 0.0);

    double total_electron_charge_to_add = 0.0;
    double total_hole_charge_to_add     = 0.0;

    for (std::size_t i = 0; i < number_contact_elements; ++i) {
        auto& element = m_list_element_contact_ptr[i];

        const double element_charge             = element->get_n_charge() - element->get_p_charge();
        const double accumulated_element_charge = element_charge / static_cast<double>(poisson_frequency);
        const double equilibrium_charge         = m_list_element_contact_equilibrium_charge[i];
        const double charge_to_add              = equilibrium_charge - accumulated_element_charge;

        if (equilibrium_charge > 0.0) {
            electron_charge_to_add[i] = charge_to_add;
            total_electron_charge_to_add += charge_to_add;
        } else {
            hole_charge_to_add[i] = charge_to_add;
            total_hole_charge_to_add += charge_to_add;
        }
    }

    const auto draw_contact_index = [&](const std::vector<double>& charge) {
        std::vector<std::size_t> candidates;
        candidates.reserve(charge.size());

        for (std::size_t i = 0; i < charge.size(); ++i) {
            if (charge[i] != 0.0) {
                candidates.push_back(i);
            }
        }

        if (candidates.empty()) {
            return std::optional<std::size_t>{};
        }

        std::uniform_int_distribution<std::size_t> distribution(0, candidates.size() - 1);
        return std::optional<std::size_t>{candidates[distribution(m_contact_rng)]};
    };

    std::size_t number_electrons_to_place = static_cast<std::size_t>(std::max(0.0, total_electron_charge_to_add));
    std::size_t number_holes_to_place     = static_cast<std::size_t>(std::max(0.0, -total_hole_charge_to_add));

    while (electron_positions.size() < number_electrons_to_place) {
        const auto index = draw_contact_index(electron_charge_to_add);
        if (!index.has_value()) {
            break;
        }
        const std::size_t i = index.value();
        if (electron_charge_to_add[i] <= 0.0) {
            electron_charge_to_add[i] = 0.0;
            continue;
        }
        electron_positions.push_back(m_list_element_contact_ptr[i]->draw_uniform_random_point_inside_element(m_contact_rng));
        electron_charge_to_add[i] -= 1.0;
    }

    while (hole_positions.size() < number_holes_to_place) {
        const auto index = draw_contact_index(hole_charge_to_add);
        if (!index.has_value()) {
            break;
        }
        const std::size_t i = index.value();
        if (hole_charge_to_add[i] >= 0.0) {
            hole_charge_to_add[i] = 0.0;
            continue;
        }
        hole_positions.push_back(m_list_element_contact_ptr[i]->draw_uniform_random_point_inside_element(m_contact_rng));
        hole_charge_to_add[i] += 1.0;
    }

    add_particles_at_positions(electron_positions, particle_type::electron);
    add_particles_at_positions(hole_positions, particle_type::hole);
}

void self_consistent_device_amc_simulation_3d::compute_unitary_potential() {
    m_poisson_solver.compute_second_member(0.0);
    m_poisson_solver.apply_dirichlet_condition("anode", 1.0);
    m_poisson_solver.apply_dirichlet_condition("cathode", 0.0);
    m_poisson_solver.decompose_matrix();
    m_poisson_solver.solve_system();

    constexpr bool add_gradient = true;
    m_poisson_solver.add_solution_to_mesh_functions("RamoUnitaryPotential", add_gradient);
}

void self_consistent_device_amc_simulation_3d::initialize_poisson_solver() {
    m_poisson_solver.compute_stiffness_matrix();
    compute_unitary_potential();
}

self_consistent_device_amc_simulation_3d::self_consistent_device_amc_simulation_3d(
    const device::device&                        simulation_device,
    const options_device_amc&                    simulation_options,
    const options_self_consistent_device_amc_3d& self_consistent_options,
    const physic::material::list_materials&      list_materials,
    const std::string&                           simulation_name,
    int                                          seed_random_generator)
    : device_amc_simulation(simulation_device, simulation_options, simulation_name, seed_random_generator),
      m_self_consistent_options(self_consistent_options),
      m_poisson_solver(m_device.get_p_mesh(), m_device.get_p_mesh()->get_nb_vertices(), list_materials),
      m_contact_rng(seed_random_generator + 1) {
    validate_self_consistent_options();
    initialize_contact_elements();
    initialize_poisson_solver();
}

self_consistent_device_amc_simulation_3d::self_consistent_device_amc_simulation_3d(
    const device::device&                        simulation_device,
    const options_device_amc&                    simulation_options,
    const options_self_consistent_device_amc_3d& self_consistent_options,
    const physic::material::list_materials&      list_materials,
    const std::string&                           simulation_name,
    const mesh::vector3&                         starting_position,
    std::size_t                                  number_electrons_start,
    std::size_t                                  number_holes_start,
    int                                          seed_random_generator)
    : device_amc_simulation(simulation_device,
                            simulation_options,
                            simulation_name,
                            starting_position,
                            number_electrons_start,
                            number_holes_start,
                            seed_random_generator),
      m_self_consistent_options(self_consistent_options),
      m_poisson_solver(m_device.get_p_mesh(), m_device.get_p_mesh()->get_nb_vertices(), list_materials),
      m_contact_rng(seed_random_generator + 1) {
    validate_self_consistent_options();
    initialize_contact_elements();
    initialize_poisson_solver();
}

void self_consistent_device_amc_simulation_3d::add_particle_charges_to_elements() {
    for (const auto& p_particle : m_list_particles) {
        auto* element = p_particle->get_containing_element();
        if (element == nullptr) {
            continue;
        }
        if (p_particle->type() == particle_type::electron) {
            element->add_n_charge(p_particle->weight());
        } else {
            element->add_p_charge(p_particle->weight());
        }
    }
}

void self_consistent_device_amc_simulation_3d::reset_element_charges() {
    for (const auto& element : m_device.get_p_mesh()->get_list_bulk_element()) {
        element->reset_charge();
    }
}

void self_consistent_device_amc_simulation_3d::recompute_vertex_space_charge_from_element_charges(std::size_t accumulation_steps) {
    if (accumulation_steps == 0) {
        throw std::invalid_argument("accumulation_steps must be positive.");
    }
    const double factor = 1.0 / static_cast<double>(accumulation_steps);
    m_device.get_p_mesh()->convert_charge_on_element_into_charge_at_vtx(factor);
}

void self_consistent_device_amc_simulation_3d::update_self_consistent_potential() {
    m_poisson_solver.update_second_member();
    m_poisson_solver.apply_dirichlet_condition_second_member("anode", m_self_consistent_options.m_anode_voltage);
    m_poisson_solver.apply_dirichlet_condition_second_member("cathode", m_self_consistent_options.m_cathode_voltage);
    m_poisson_solver.solve_system();
    constexpr bool add_gradient = true;
    m_poisson_solver.add_solution_to_mesh_functions("DorySolution", add_gradient);
}

void self_consistent_device_amc_simulation_3d::run_self_consistent_transport_simulation() {
    const std::size_t total_iterations =
        static_cast<std::size_t>(std::ceil(m_simulation_options.m_t_max / m_simulation_options.m_time_step));

    fmt::print("START 3D SELF-CONSISTENT AMC SIMULATION\n");
    fmt::print("Total iterations: {}\n", total_iterations);
    fmt::print("Poisson frequency: {}\n", m_self_consistent_options.m_poisson_frequency);

    while (m_time <= m_simulation_options.m_t_max && !m_list_particles.empty()) {
        if (m_simulation_options.m_stop_simu_when_no_electron_remaining && get_number_electrons() == 0) {
            fmt::print("Stop: no electrons remaining in device.\n");
            break;
        }

        if (has_reached_avalanche()) {
            fmt::print("Stop: avalanche threshold reached.\n");
            break;
        }

        transport_particles_one_time_step();
        add_particle_charges_to_elements();

        const bool should_update_poisson = (m_iteration % m_self_consistent_options.m_poisson_frequency == 0) && (m_iteration != 0);

        if (should_update_poisson) {
            add_charges_at_contacts(m_self_consistent_options.m_poisson_frequency);
            add_particle_charges_to_elements();
            recompute_vertex_space_charge_from_element_charges(m_self_consistent_options.m_poisson_frequency + 1);
            update_self_consistent_potential();
            reset_element_charges();
        }

        
        m_time += m_simulation_options.m_time_step;
        ++m_iteration;
        
        m_simulation_history.add_data_to_history(m_time,
                                                 get_number_electrons(),
                                                 get_number_holes(),
                                                 m_simulation_history.m_impact_ionization_positions.size(),
                                                 m_anode_current,
                                                 m_cathode_current,
                                                 0.0);
    }

    fmt::print("END 3D SELF-CONSISTENT AMC SIMULATION\n");
}

}  // namespace uepm::amc