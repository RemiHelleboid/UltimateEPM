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

    if (common_options().m_poisson_frequency == 0) {
        throw std::invalid_argument("Poisson frequency must be positive.");
    }
}

/**
 * @brief We need to identify the elements adjacent to the contacts, and compute their equilibrium charge (integral of
 * doping concentration) to be able to properly update their charge during the simulation (so that they are not fixed at
 * the equilibrium charge, but can exchange charge with the device while still having a reference equilibrium charge for
 * the poisson solver).
 *
 */
void self_consistent_device_amc_simulation_3d::initialize_contact_elements() {
    const auto list_contact_elem_anode =
        m_device.get_p_mesh()->get_idx_bulk_elements_adjacent_to_contact_region("anode");

    const auto list_contact_elem_cathode =
        m_device.get_p_mesh()->get_idx_bulk_elements_adjacent_to_contact_region("cathode");

    const auto list_bulk_elements = m_device.get_p_mesh()->get_list_bulk_element();

    m_list_element_contact.clear();
    m_list_element_contact_ptr.clear();
    m_list_element_contact_equilibrium_charge.clear();

    m_list_element_contact.reserve(list_contact_elem_anode.size() + list_contact_elem_cathode.size());
    m_list_element_contact_ptr.reserve(list_contact_elem_anode.size() + list_contact_elem_cathode.size());
    m_list_element_contact_equilibrium_charge.reserve(list_contact_elem_anode.size() +
                                                      list_contact_elem_cathode.size());

    const auto add_contact_element = [&](std::size_t element_index) {
        if (element_index >= list_bulk_elements.size()) {
            throw std::runtime_error("Contact-adjacent bulk element index is out of range.");
        }

        auto element = list_bulk_elements[element_index];

        m_list_element_contact.push_back(element_index);
        m_list_element_contact_ptr.push_back(element);

        const double equilibrium_charge = element->integrate_scalar("DopingConcentration");

        m_list_element_contact_equilibrium_charge.push_back(equilibrium_charge);
        m_list_element_contact_residual_charge.push_back(0.0);
    };

    for (const auto element_index : list_contact_elem_anode) {
        add_contact_element(element_index);
    }

    for (const auto element_index : list_contact_elem_cathode) {
        add_contact_element(element_index);
    }
}

void self_consistent_device_amc_simulation_3d::place_initial_charges_according_to_doping(double particle_weight) {
    if (particle_weight <= 0.0) {
        throw std::invalid_argument("particle_weight must be positive.");
    }

    constexpr double micron3_to_cm3 = 1.0e-12;

    const std::string donor_field_name    = "DonorConcentration";
    const std::string acceptor_field_name = "AcceptorConcentration";

    auto* mesh = m_device.get_p_mesh();

    const double total_donor_charge = micron3_to_cm3 * mesh->integrate_over_mesh(donor_field_name);

    const double total_acceptor_charge = micron3_to_cm3 * mesh->integrate_over_mesh(acceptor_field_name);

    const std::size_t number_electrons = static_cast<std::size_t>(std::floor(total_donor_charge / particle_weight));

    const std::size_t number_holes = static_cast<std::size_t>(std::floor(total_acceptor_charge / particle_weight));

    fmt::print("Initial doping charge:\n");
    fmt::print("  donor carriers:    {:.6e}\n", total_donor_charge);
    fmt::print("  acceptor carriers: {:.6e}\n", total_acceptor_charge);
    fmt::print("Initial numerical particles:\n");
    fmt::print("  particle weight: {:.6e}\n", particle_weight);
    fmt::print("  electrons: {}\n", number_electrons);
    fmt::print("  holes:     {}\n", number_holes);

    if (number_electrons == 0 && number_holes == 0) {
        fmt::print("No initial doping particles created. Decrease particle weight.\n");
        return;
    }

    std::vector<mesh::vector3> electron_positions;
    std::vector<mesh::vector3> hole_positions;

    electron_positions.reserve(number_electrons);
    hole_positions.reserve(number_holes);

    std::uniform_real_distribution<double> uniform01(0.0, 1.0);

    const double max_donor_concentration = mesh->get_argmax_max_of_function(donor_field_name).second;

    const double max_acceptor_concentration = mesh->get_argmax_max_of_function(acceptor_field_name).second;

    if (number_electrons > 0 && max_donor_concentration <= 0.0) {
        throw std::runtime_error("Donor concentration maximum is non-positive.");
    }

    if (number_holes > 0 && max_acceptor_concentration <= 0.0) {
        throw std::runtime_error("Acceptor concentration maximum is non-positive.");
    }

    const mesh::bbox device_bbox = mesh->get_bounding_box();

    while (electron_positions.size() < number_electrons) {
        const mesh::vector3 position         = device_bbox.draw_uniform_random_point_inside_box(m_contact_rng);
        const double        donor_density    = mesh->interpolate_scalar_at_location(donor_field_name, position);
        const double        acceptor_density = mesh->interpolate_scalar_at_location(acceptor_field_name, position);
        if (acceptor_density > donor_density) {
            continue;
        }

        const double probability = donor_density / max_donor_concentration;

        if (uniform01(m_contact_rng) < probability) {
            electron_positions.push_back(position);
        }
    }

    while (hole_positions.size() < number_holes) {
        const mesh::vector3 position         = device_bbox.draw_uniform_random_point_inside_box(m_contact_rng);
        const double        acceptor_density = mesh->interpolate_scalar_at_location(acceptor_field_name, position);
        const double        donor_density    = mesh->interpolate_scalar_at_location(donor_field_name, position);

        if (donor_density > acceptor_density) {
            continue;
        }

        const double probability = acceptor_density / max_acceptor_concentration;
        if (uniform01(m_contact_rng) < probability) {
            hole_positions.push_back(position);
        }
    }

    m_list_particles.clear();
    m_list_particles.reserve(electron_positions.size() + hole_positions.size());

    add_particles_at_positions(electron_positions, particle_type::electron, particle_weight);

    add_particles_at_positions(hole_positions, particle_type::hole, particle_weight);

    reset_element_charges();
    add_particle_charges_to_elements();
    recompute_vertex_space_charge_from_element_charges(1);
    update_self_consistent_potential();
    reset_element_charges();

    fmt::print("Initial particles placed according to doping.\n");
}

void self_consistent_device_amc_simulation_3d::add_charges_at_contacts(std::size_t poisson_frequency) {
    if (poisson_frequency == 0) {
        throw std::invalid_argument("Poisson frequency must be positive.");
    }

    const double particle_weight = common_options().m_contact_injection_particle_weight;

    if (particle_weight <= 0.0) {
        throw std::invalid_argument("Contact injection particle weight must be positive.");
    }

    const std::size_t number_contact_elements = m_list_element_contact_ptr.size();

    std::vector<double> electron_charge_to_add(number_contact_elements, 0.0);
    std::vector<double> hole_charge_to_add(number_contact_elements, 0.0);

    double total_electron_charge_to_add = 0.0;
    double total_hole_charge_to_add     = 0.0;

    for (std::size_t i = 0; i < number_contact_elements; ++i) {
        auto& element = m_list_element_contact_ptr[i];

        const double element_charge          = element->get_n_charge() - element->get_p_charge();
        const double averaged_element_charge = element_charge / static_cast<double>(poisson_frequency);
        const double equilibrium_charge      = m_list_element_contact_equilibrium_charge[i];
        const double charge_to_add           = equilibrium_charge - averaged_element_charge;

        if (equilibrium_charge > 0.0 && charge_to_add > 0.0) {
            electron_charge_to_add[i] = charge_to_add;
            total_electron_charge_to_add += charge_to_add;
        } else if (equilibrium_charge < 0.0 && charge_to_add < 0.0) {
            hole_charge_to_add[i] = -charge_to_add;
            total_hole_charge_to_add += -charge_to_add;
        }
    }

    const std::size_t number_electrons_to_place =
        static_cast<std::size_t>(std::floor(total_electron_charge_to_add / particle_weight));

    const std::size_t number_holes_to_place =
        static_cast<std::size_t>(std::floor(total_hole_charge_to_add / particle_weight));

    std::vector<mesh::vector3> electron_positions;
    std::vector<mesh::vector3> hole_positions;

    electron_positions.reserve(number_electrons_to_place);
    hole_positions.reserve(number_holes_to_place);

    const auto has_capacity = [&]() {
        return m_list_particles.size() + electron_positions.size() + hole_positions.size() <
               m_simulation_options.m_max_number_particle;
    };

    std::uniform_int_distribution<std::size_t> contact_index_distribution(0, number_contact_elements - 1);

    while (electron_positions.size() < number_electrons_to_place && has_capacity()) {
        const std::size_t i = contact_index_distribution(m_contact_rng);

        if (electron_charge_to_add[i] <= 0.0) {
            continue;
        }

        electron_positions.push_back(
            m_list_element_contact_ptr[i]->draw_uniform_random_point_inside_element(m_contact_rng));

        electron_charge_to_add[i] -= particle_weight;
    }

    while (hole_positions.size() < number_holes_to_place && has_capacity()) {
        const std::size_t i = contact_index_distribution(m_contact_rng);

        if (hole_charge_to_add[i] <= 0.0) {
            continue;
        }

        hole_positions.push_back(
            m_list_element_contact_ptr[i]->draw_uniform_random_point_inside_element(m_contact_rng));

        hole_charge_to_add[i] -= particle_weight;
    }

    add_particles_at_positions(electron_positions, particle_type::electron, particle_weight);
    add_particles_at_positions(hole_positions, particle_type::hole, particle_weight);
}

void self_consistent_device_amc_simulation_3d::add_missing_contact_charge_to_poisson_reservoir(
    std::size_t accumulation_steps) {
    if (accumulation_steps == 0) {
        throw std::invalid_argument("accumulation_steps must be positive.");
    }

    const double accumulation_factor = static_cast<double>(accumulation_steps);

    for (std::size_t i = 0; i < m_list_element_contact_ptr.size(); ++i) {
        auto& element = m_list_element_contact_ptr[i];

        const double equilibrium_charge        = m_list_element_contact_equilibrium_charge[i];
        const double accumulated_mobile_charge = element->get_n_charge() - element->get_p_charge();
        const double target_accumulated_charge = equilibrium_charge * accumulation_factor;
        const double correction                = target_accumulated_charge - accumulated_mobile_charge;
        // fmt::print(
        //     "Contact element {}: equilibrium charge = {:.3e}, accumulated mobile charge = {:.3e}, target accumulated
        //     " "charge = {:.3e}, correction = {:.3e}\n", m_list_element_contact[i], equilibrium_charge,
        //     accumulated_mobile_charge,
        //     target_accumulated_charge,
        //     correction);

        if (correction > 0.0) {
            element->add_n_charge(correction);
        } else if (correction < 0.0) {
            element->add_p_charge(-correction);
        }
    }
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
    int                                          seed_random_generator)
    : self_consistent_device_amc_simulation_base(simulation_device,
                                                 simulation_options,
                                                 self_consistent_options.m_common,
                                                 seed_random_generator),
      m_self_consistent_options(self_consistent_options),
      m_poisson_solver(m_device.get_p_mesh(), m_device.get_p_mesh()->get_nb_vertices(), list_materials),
      m_contact_rng(seed_random_generator + 1) {
    validate_self_consistent_options();
    initialize_contact_elements();
    initialize_poisson_solver();
    place_initial_charges_according_to_doping();
}

self_consistent_device_amc_simulation_3d::self_consistent_device_amc_simulation_3d(
    const device::device&                        simulation_device,
    const options_device_amc&                    simulation_options,
    const options_self_consistent_device_amc_3d& self_consistent_options,
    const physic::material::list_materials&      list_materials,
    const mesh::vector3&                         starting_position,
    std::size_t                                  number_electrons_start,
    std::size_t                                  number_holes_start,
    int                                          seed_random_generator)
    : self_consistent_device_amc_simulation_base(simulation_device,
                                                 simulation_options,
                                                 self_consistent_options.m_common,
                                                 starting_position,
                                                 number_electrons_start,
                                                 number_holes_start,
                                                 seed_random_generator),
      m_self_consistent_options(self_consistent_options),
      m_poisson_solver(m_device.get_p_mesh(), m_device.get_p_mesh()->get_nb_vertices(), list_materials) {
    validate_self_consistent_options();
    initialize_contact_elements();
    initialize_poisson_solver();
    place_initial_charges_according_to_doping();
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

void self_consistent_device_amc_simulation_3d::recompute_vertex_space_charge_from_element_charges(
    std::size_t accumulation_steps) {
    if (accumulation_steps == 0) {
        throw std::invalid_argument("accumulation_steps must be positive.");
    }
    const double factor = 1.0 / static_cast<double>(accumulation_steps);
    m_device.get_p_mesh()->convert_charge_on_element_into_charge_at_vtx(factor);
}

void self_consistent_device_amc_simulation_3d::update_self_consistent_potential() {
    m_poisson_solver.update_second_member();
    m_poisson_solver.apply_dirichlet_condition_second_member("anode", anode_voltage_for_poisson());
    m_poisson_solver.apply_dirichlet_condition_second_member("cathode", cathode_voltage_for_poisson());
    m_poisson_solver.solve_system();
    constexpr bool add_gradient = true;
    m_poisson_solver.add_solution_to_mesh_functions("PoissonSolution", add_gradient);
}

void self_consistent_device_amc_simulation_3d::run_self_consistent_transport_simulation() {
    const std::size_t total_iterations =
        static_cast<std::size_t>(std::ceil(m_simulation_options.m_t_max / m_simulation_options.m_time_step));

    fmt::print("START 3D SELF-CONSISTENT AMC SIMULATION\n");
    fmt::print("Total iterations: {}\n", total_iterations);
    fmt::print("Poisson frequency: {}\n", common_options().m_poisson_frequency);

    double accumulator_ramo_current_electron = 0.0;
    double accumulator_ramo_current_hole     = 0.0;
    double ramo_current_electron             = 0.0;
    double ramo_current_hole                 = 0.0;
    double ramo_current                      = 0.0;

    while (m_state.m_time_s <= m_simulation_options.m_t_max && !m_list_particles.empty()) {
        if (m_simulation_options.m_stop_simu_when_no_electron_remaining && get_number_electrons() == 0) {
            fmt::print("Stop: no electrons remaining in device.\n");
            break;
        }

        if (has_reached_particle_limit()) {
            fmt::print("Stop: hard particle limit reached.\n");
            break;
        }

        const std::size_t impact_events_before_step =
            m_simulation_history.m_impact_ionization_positions.size();
        transport_particles_one_time_step();
        add_particle_charges_to_elements();
        const auto [electron_current, hole_current] = compute_ramo_current();
        accumulator_ramo_current_electron += electron_current;
        accumulator_ramo_current_hole += hole_current;

        const bool should_update_poisson =
            (m_state.m_iteration % common_options().m_poisson_frequency == 0) && (m_state.m_iteration != 0);

        if (should_update_poisson) {
            const std::size_t poisson_frequency = common_options().m_poisson_frequency;

            ramo_current_electron =
                accumulator_ramo_current_electron / static_cast<double>(poisson_frequency);
            ramo_current_hole = accumulator_ramo_current_hole / static_cast<double>(poisson_frequency);
            ramo_current      = ramo_current_electron + ramo_current_hole;
            ramo_current -= common_options().m_background_ramo_current_A;
            accumulator_ramo_current_electron = 0.0;
            accumulator_ramo_current_hole     = 0.0;

            if (m_simulation_options.m_scheduled_particle_injection.m_done) {
                const double circuit_dt_s =
                    m_simulation_options.m_time_step * static_cast<double>(poisson_frequency);
                const double sample_time_s = m_state.m_time_s + m_simulation_options.m_time_step;
                advance_quench_circuit(ramo_current, circuit_dt_s, sample_time_s);
            }

            add_charges_at_contacts(poisson_frequency);
            add_particle_charges_to_elements();
            add_missing_contact_charge_to_poisson_reservoir(poisson_frequency + 1);
            recompute_vertex_space_charge_from_element_charges(poisson_frequency + 1);
            update_self_consistent_potential();
            reset_element_charges();
        }

        m_state.m_time_s += m_simulation_options.m_time_step;
        ++m_state.m_iteration;
        update_successful_quench_detection(m_state.m_time_s, impact_events_before_step);

        const double max_electric_field_V_per_cm = max_particle_electric_field_V_per_cm();
        m_simulation_history.add_data_to_history(m_state.m_time_s,
                                                 get_number_electrons(),
                                                 get_number_holes(),
                                                 m_simulation_history.m_impact_ionization_positions.size(),
                                                 ramo_current_electron,
                                                 ramo_current_hole,
                                                 ramo_current,
                                                 max_electric_field_V_per_cm,
                                                 anode_voltage_for_poisson(),
                                                 cathode_voltage_for_poisson(),
                                                 quench_supply_voltage_for_history(),
                                                 quench_device_current_for_history(),
                                                 quench_resistor_current_for_history(),
                                                 quench_voltage_drop_for_history());

        if (m_state.m_iteration == 1 ||
            (m_simulation_options.m_export_time_step &&
             m_state.m_iteration % static_cast<std::size_t>(m_simulation_options.m_frequency_export_trajectory) == 0)) {
            export_current_state();
            fmt::print("\rExported iteration at time {:<10.3e}ps - {:>9d} / {} ({:.1f}%) -- number of particles: {}",
                       m_state.m_time_s * 1e12,
                       m_state.m_iteration,
                       total_iterations,
                       static_cast<double>(m_state.m_iteration) / static_cast<double>(total_iterations) * 100.0,
                       m_list_particles.size());
            std::fflush(stdout);
        }
    }

    fmt::print("END 3D SELF-CONSISTENT AMC SIMULATION\n");
}

void self_consistent_device_amc_simulation_3d::export_current_state() { export_current_snapshot(); }

}  // namespace uepm::amc
