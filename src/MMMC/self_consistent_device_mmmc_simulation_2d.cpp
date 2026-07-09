/**
 * @file self_consistent_device_mmmc_simulation_2d.cpp
 * @brief Self-consistent 2D Mixed Method Monte Carlo device simulation.
 */

#include "self_consistent_device_mmmc_simulation_2d.hpp"

#include <fmt/core.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

#include "boundary_reflection.hpp"
#include "mmmc_particle_transfer.hpp"
#include "physical_constants.hpp"
#include "unit_conversion.hpp"

namespace uepm::MMMC {
namespace {

double finite_or_zero(double value) {
    return std::isfinite(value) ? value : 0.0;
}

double local_net_doping_cm_3(const mesh::element& element, const mesh::vector3& position_um) {
    return finite_or_zero(element.interpolate_doping_at_location(position_um));
}

double admc_signed_charge_C(const ADMC::device_admc_particle& particle) {
    return particle.weight * ADMC::carrier_charge_sign(particle.particle.type()) * constants::q_e;
}

}  // namespace

void options_device_MMMC::synchronize_from_pbmc() {
    m_admc.m_simulation_name                  = m_pbmc.m_simulation_name;
    m_admc.m_output_directory                 = m_pbmc.m_output_directory;
    m_admc.m_lattice_temperature_K            = m_pbmc.m_lattice_temperature;
    m_admc.m_time_step_s                      = m_pbmc.m_time_step;
    m_admc.m_final_time_s                     = m_pbmc.m_t_max;
    m_admc.m_max_number_particles             = m_pbmc.m_max_number_particle;
    m_admc.m_stop_when_no_electrons           = m_pbmc.m_stop_simu_when_no_electron_remaining;
    m_admc.m_export_time_step                 = m_pbmc.m_export_time_step;
    m_admc.m_frequency_export                 = m_pbmc.m_frequency_export_trajectory;
    m_admc.m_prefix_export_filename           = m_pbmc.m_prefix_export_filename;
    m_admc.m_boundary_reflection_model        = m_pbmc.m_boundary_reflection_model;
    m_admc.m_current_probe.m_enabled          = m_pbmc.m_current_probe.m_enabled;
    m_admc.m_current_probe.m_box_um           = m_pbmc.m_current_probe.m_box_um;
}

void options_device_MMMC::validate() const {
    m_pbmc.validate();
    m_admc.validate();
    if (m_pbmc.m_time_step != m_admc.m_time_step_s) {
        throw std::invalid_argument("MMMC requires identical PBMC and ADMC time steps.");
    }
    if (m_pbmc.m_t_max != m_admc.m_final_time_s) {
        throw std::invalid_argument("MMMC requires identical PBMC and ADMC final times.");
    }
}

void options_self_consistent_device_MMMC_2d::validate() const {
    m_common.validate();
    m_policy.validate();
    if (!std::isfinite(m_effective_depth_um) || m_effective_depth_um <= 0.0) {
        throw std::invalid_argument("MMMC effective depth must be positive and finite.");
    }
}

self_consistent_device_mmmc_simulation_2d::self_consistent_device_mmmc_simulation_2d(
    const device::device&                         simulation_device,
    const options_device_MMMC&                    simulation_options,
    const options_self_consistent_device_MMMC_2d& self_consistent_options,
    const physics::material_database&             material_database,
    int                                           seed_random_generator)
    : PBMC::self_consistent_device_pbmc_simulation_base(simulation_device,
                                                        simulation_options.m_pbmc,
                                                        self_consistent_options.m_common,
                                                        seed_random_generator),
      m_mmmc_options(simulation_options),
      m_self_consistent_options(self_consistent_options),
      m_poisson_solver(m_device.get_p_mesh(), m_device.get_p_mesh()->get_nb_vertices(), material_database),
      m_contact_rng(seed_random_generator + 1),
      m_admc_rng(static_cast<std::uint64_t>(seed_random_generator) + 7919u) {
    m_mmmc_options.synchronize_from_pbmc();
    validate_self_consistent_options();
    apply_scheduled_contact_voltage_events(0.0);
    initialize_contact_elements();
    initialize_poisson_solver();
    initialize_particles_for_self_consistent_run();
}

self_consistent_device_mmmc_simulation_2d::self_consistent_device_mmmc_simulation_2d(
    const device::device&                         simulation_device,
    const options_device_MMMC&                    simulation_options,
    const options_self_consistent_device_MMMC_2d& self_consistent_options,
    const physics::material_database&             material_database,
    const mesh::vector3&                          starting_position,
    std::size_t                                   number_electrons_start,
    std::size_t                                   number_holes_start,
    int                                           seed_random_generator)
    : PBMC::self_consistent_device_pbmc_simulation_base(simulation_device,
                                                        simulation_options.m_pbmc,
                                                        self_consistent_options.m_common,
                                                        starting_position,
                                                        number_electrons_start,
                                                        number_holes_start,
                                                        seed_random_generator),
      m_mmmc_options(simulation_options),
      m_self_consistent_options(self_consistent_options),
      m_poisson_solver(m_device.get_p_mesh(), m_device.get_p_mesh()->get_nb_vertices(), material_database),
      m_contact_rng(seed_random_generator + 1),
      m_admc_rng(static_cast<std::uint64_t>(seed_random_generator) + 7919u) {
    m_mmmc_options.synchronize_from_pbmc();
    validate_self_consistent_options();
    apply_scheduled_contact_voltage_events(0.0);
    initialize_contact_elements();
    initialize_poisson_solver();
    initialize_particles_for_self_consistent_run();
}

void self_consistent_device_mmmc_simulation_2d::validate_self_consistent_options() const {
    if (m_dimension != 2) {
        throw std::invalid_argument("self_consistent_device_mmmc_simulation_2d requires a 2D device.");
    }
    m_mmmc_options.validate();
    m_self_consistent_options.validate();
}

double self_consistent_device_mmmc_simulation_2d::scale_integrated_2d_doping_to_carriers(
    double integrated_doping) const {
    return integrated_doping * m_self_consistent_options.m_effective_depth_um * units::micron_to_cm;
}

double self_consistent_device_mmmc_simulation_2d::charge_deposition_factor(std::size_t accumulation_steps) const {
    if (accumulation_steps == 0) {
        throw std::invalid_argument("MMMC accumulation steps must be positive.");
    }
    return 1.0 / (static_cast<double>(accumulation_steps) * m_self_consistent_options.m_effective_depth_um);
}

double self_consistent_device_mmmc_simulation_2d::ramo_current_scale_factor() const {
    return 1.0 / m_self_consistent_options.m_effective_depth_um;
}

double self_consistent_device_mmmc_simulation_2d::current_density_cell_volume_m3(
    const mesh::element& element) const {
    const double area_um2             = std::abs(element.get_measure());
    const double effective_volume_um3 = area_um2 * m_self_consistent_options.m_effective_depth_um;
    return effective_volume_um3 * std::pow(units::micron_to_meter, 3);
}

void self_consistent_device_mmmc_simulation_2d::initialize_contact_elements() {
    const auto list_bulk_elements = m_device.get_p_mesh()->get_list_bulk_element();

    m_list_element_contact.clear();
    m_list_element_contact_ptr.clear();
    m_list_element_contact_equilibrium_charge.clear();

    for (const auto& device_contact : m_device.get_list_contacts()) {
        const std::string contact_name = device_contact.get_contact_name();
        if (!contact_voltages_V().contains(contact_name)) {
            throw std::runtime_error("MMMC particle-collection contact '" + contact_name +
                                     "' has no configured Poisson voltage.");
        }
        const auto element_indices =
            m_device.get_p_mesh()->get_idx_bulk_elements_adjacent_to_contact_region(contact_name);
        std::vector<std::shared_ptr<mesh::element>> contact_elements;
        contact_elements.reserve(element_indices.size());

        for (const auto element_index : element_indices) {
            if (element_index >= list_bulk_elements.size()) {
                throw std::runtime_error("MMMC contact-adjacent bulk element index is out of range.");
            }
            auto element = list_bulk_elements[element_index];
            if (!is_transport_material_element(*element)) {
                continue;
            }
            contact_elements.push_back(element);
            m_list_element_contact.push_back(element_index);
            m_list_element_contact_ptr.push_back(element);
            m_list_element_contact_equilibrium_charge.push_back(
                scale_integrated_2d_doping_to_carriers(element->integrate_scalar("DopingConcentration")));
        }

        if (!contact_elements.empty()) {
            update_built_in_contact_voltage_offset(contact_name, contact_elements);
        }
    }
}

void self_consistent_device_mmmc_simulation_2d::compute_unitary_potential() {
    m_poisson_solver.compute_second_member(0.0);
    for (const auto& [contact_name, unused_voltage] : contact_voltages_V()) {
        static_cast<void>(unused_voltage);
        m_poisson_solver.apply_dirichlet_condition(contact_name, contact_name == ramo_electrode() ? 1.0 : 0.0);
    }
    m_poisson_solver.decompose_matrix();
    m_poisson_solver.solve_system();
    m_poisson_solver.add_solution_to_mesh_functions("RamoUnitaryPotential", true);

    const double unitary_potential_gradient_max =
        m_device.get_p_mesh()->get_argmax_max_of_function("RamoUnitaryPotential_gradient_norm").second;
    const double unitary_potential_gradient_min =
        m_device.get_p_mesh()->get_argmin_min_of_function("RamoUnitaryPotential_gradient_norm").second;
    if (std::abs(unitary_potential_gradient_max - unitary_potential_gradient_min) < 1e-6) {
        m_state.m_use_constant_RamoUnitaryElectricField = true;
        m_state.m_RamoUnitaryElectricField_Vm_per_cm =
            m_device.get_p_mesh()->interpolate_vector_at_location("RamoUnitaryPotential_gradient", {1e-3, 1e-3, 0.0});
    }
}

void self_consistent_device_mmmc_simulation_2d::initialize_poisson_solver() {
    m_poisson_solver.compute_stiffness_matrix();
    compute_unitary_potential();

    m_poisson_solver.compute_stiffness_matrix();
    m_poisson_solver.compute_second_member(0.0);
    for (const auto& [contact_name, unused_voltage] : contact_voltages_V()) {
        static_cast<void>(unused_voltage);
        m_poisson_solver.apply_dirichlet_condition(contact_name, contact_voltage_for_poisson(contact_name));
    }
    m_poisson_solver.decompose_matrix();
}

void self_consistent_device_mmmc_simulation_2d::initialize_particles_for_self_consistent_run() {
    if (!common_options().m_initial_particle_state_file.empty()) {
        load_particles_from_state_csv(common_options().m_initial_particle_state_file);
        flatten_particle_positions_for_2d();
    } else if (common_options().m_initialize_particles_from_doping) {
        place_initial_charges_according_to_doping(common_options().m_initial_particle_weight);
    }

    apply_transport_policy();
    reset_element_charges();
    add_particle_charges_to_elements();
    recompute_vertex_space_charge_from_element_charges(1);
    update_self_consistent_potential();
    reset_element_charges();
}

void self_consistent_device_mmmc_simulation_2d::place_initial_charges_according_to_doping(double particle_weight) {
    if (particle_weight <= 0.0) {
        throw std::invalid_argument("MMMC particle weight must be positive.");
    }

    const std::string donor_field_name    = "DonorConcentration";
    const std::string acceptor_field_name = "AcceptorConcentration";
    auto*             mesh                = m_device.get_p_mesh();
    const auto integrate_carriers_over_2d_mesh = [&](const std::string& field_name) {
        double total_charge = 0.0;
        for (const auto& element : mesh->get_list_bulk_element()) {
            total_charge += scale_integrated_2d_doping_to_carriers(element->integrate_scalar(field_name));
        }
        return total_charge;
    };

    const double      total_donor_charge    = integrate_carriers_over_2d_mesh(donor_field_name);
    const double      total_acceptor_charge = integrate_carriers_over_2d_mesh(acceptor_field_name);
    const std::size_t number_electrons = static_cast<std::size_t>(std::floor(total_donor_charge / particle_weight));
    const std::size_t number_holes     = static_cast<std::size_t>(std::floor(total_acceptor_charge / particle_weight));

    constexpr double min_probability = 50e-2;
    const double max_donor_concentration = mesh->get_argmax_max_of_function(donor_field_name).second;
    const double max_acceptor_concentration = mesh->get_argmax_max_of_function(acceptor_field_name).second;
    const mesh::bbox active_region_bbox = mesh->get_p_region("Silicon_1")->compute_bounding_box();
    std::uniform_real_distribution<double> uniform01(0.0, 1.0);

    std::vector<mesh::vector3> electron_positions;
    std::vector<mesh::vector3> hole_positions;
    electron_positions.reserve(number_electrons);
    hole_positions.reserve(number_holes);

    const std::size_t max_attempts_electrons = 1000 * std::max<std::size_t>(number_electrons, 1);
    for (std::size_t attempts = 0; electron_positions.size() < number_electrons && attempts < max_attempts_electrons;
         ++attempts) {
        const mesh::vector3 position = active_region_bbox.draw_uniform_random_point_inside_box(m_contact_rng);
        const double donor_density = mesh->interpolate_scalar_at_location(donor_field_name, position);
        const double acceptor_density = mesh->interpolate_scalar_at_location(acceptor_field_name, position);
        if (acceptor_density > donor_density || max_donor_concentration <= 0.0) {
            continue;
        }
        const double probability = donor_density / max_donor_concentration;
        if (probability > min_probability && uniform01(m_contact_rng) < probability) {
            electron_positions.push_back(position);
        }
    }

    const std::size_t max_attempts_holes = 1000 * std::max<std::size_t>(number_holes, 1);
    for (std::size_t attempts = 0; hole_positions.size() < number_holes && attempts < max_attempts_holes;
         ++attempts) {
        const mesh::vector3 position = active_region_bbox.draw_uniform_random_point_inside_box(m_contact_rng);
        const double acceptor_density = mesh->interpolate_scalar_at_location(acceptor_field_name, position);
        const double donor_density = mesh->interpolate_scalar_at_location(donor_field_name, position);
        if (donor_density > acceptor_density || max_acceptor_concentration <= 0.0) {
            continue;
        }
        const double probability = acceptor_density / max_acceptor_concentration;
        if (probability > min_probability && uniform01(m_contact_rng) < probability) {
            hole_positions.push_back(position);
        }
    }

    add_particles_at_positions(electron_positions, PBMC::particle_type::electron, particle_weight);
    add_particles_at_positions(hole_positions, PBMC::particle_type::hole, particle_weight);
    flatten_particle_positions_for_2d();

    fmt::print("MMMC initial particles placed: PBMC pool before policy split = {}\n", m_list_particles.size());
}

void self_consistent_device_mmmc_simulation_2d::reset_element_charges() {
    m_device.get_p_mesh()->for_each_bulk_element([](mesh::element& element) { element.reset_charge(); });
}

void self_consistent_device_mmmc_simulation_2d::add_particle_charges_to_elements() {
    for (const auto& p_particle : m_list_particles) {
        auto* element = p_particle->get_containing_element();
        if (element == nullptr) {
            continue;
        }
        if (p_particle->type() == PBMC::particle_type::electron) {
            element->add_n_charge(p_particle->weight());
        } else {
            element->add_p_charge(p_particle->weight());
        }
    }

    for (const auto& particle : m_admc_particles) {
        auto* element = particle.containing_element;
        if (element == nullptr) {
            continue;
        }
        if (particle.particle.type() == ADMC::carrier_type::electron) {
            element->add_n_charge(particle.weight);
        } else {
            element->add_p_charge(particle.weight);
        }
    }
}

void self_consistent_device_mmmc_simulation_2d::recompute_vertex_space_charge_from_element_charges(
    std::size_t accumulation_steps) {
    m_device.get_p_mesh()->convert_charge_on_element_into_charge_at_vtx(charge_deposition_factor(accumulation_steps));
}

void self_consistent_device_mmmc_simulation_2d::update_self_consistent_potential(bool publish_mesh_functions) {
    m_poisson_solver.update_second_member();
    for (const auto& [contact_name, unused_voltage] : contact_voltages_V()) {
        static_cast<void>(unused_voltage);
        m_poisson_solver.apply_dirichlet_condition_second_member(contact_name,
                                                                 contact_voltage_for_poisson(contact_name));
    }
    m_poisson_solver.solve_system();
    if (common_options().m_enable_poisson_mixing) {
        if (m_previous_poisson_solution.size() == m_poisson_solver.get_solution().size() &&
            m_previous_poisson_solution.allFinite()) {
            m_poisson_solver.mix_solution_with(m_previous_poisson_solution,
                                               common_options().m_poisson_mixing_old_solution_fraction);
        }
        m_previous_poisson_solution = m_poisson_solver.get_solution();
    }
    if (publish_mesh_functions) {
        m_poisson_solver.add_solution_to_mesh_functions("PoissonSolution", true);
    } else {
        m_poisson_solver.update_mesh_electric_field_from_solution();
    }
}

void self_consistent_device_mmmc_simulation_2d::add_charges_at_contacts(std::size_t poisson_frequency_value) {
    if (poisson_frequency_value == 0) {
        throw std::invalid_argument("MMMC Poisson frequency must be positive.");
    }
    const double particle_weight = common_options().m_contact_injection_particle_weight;
    if (particle_weight <= 0.0) {
        throw std::invalid_argument("MMMC contact injection particle weight must be positive.");
    }
    if (m_list_element_contact_ptr.empty()) {
        return;
    }

    std::vector<double> electron_charge_to_add(m_list_element_contact_ptr.size(), 0.0);
    std::vector<double> hole_charge_to_add(m_list_element_contact_ptr.size(), 0.0);
    double total_electron_charge_to_add = 0.0;
    double total_hole_charge_to_add = 0.0;

    for (std::size_t i = 0; i < m_list_element_contact_ptr.size(); ++i) {
        auto& element = m_list_element_contact_ptr[i];
        const double element_charge = element->get_n_charge() - element->get_p_charge();
        const double averaged_element_charge = element_charge / static_cast<double>(poisson_frequency_value);
        const double equilibrium_charge = m_list_element_contact_equilibrium_charge[i];
        const double charge_to_add = equilibrium_charge - averaged_element_charge;
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
    std::uniform_int_distribution<std::size_t> contact_index_distribution(0, m_list_element_contact_ptr.size() - 1);

    for (std::size_t i = 0; i < number_electrons_to_place && !has_reached_particle_limit(); ++i) {
        const std::size_t contact_index = contact_index_distribution(m_contact_rng);
        if (electron_charge_to_add[contact_index] <= 0.0) {
            continue;
        }
        add_particle_at_position(m_list_element_contact_ptr[contact_index]->draw_uniform_random_point_inside_element(
                                     m_contact_rng),
                                 PBMC::particle_type::electron,
                                 particle_weight);
        electron_charge_to_add[contact_index] -= particle_weight;
    }
    for (std::size_t i = 0; i < number_holes_to_place && !has_reached_particle_limit(); ++i) {
        const std::size_t contact_index = contact_index_distribution(m_contact_rng);
        if (hole_charge_to_add[contact_index] <= 0.0) {
            continue;
        }
        add_particle_at_position(m_list_element_contact_ptr[contact_index]->draw_uniform_random_point_inside_element(
                                     m_contact_rng),
                                 PBMC::particle_type::hole,
                                 particle_weight);
        hole_charge_to_add[contact_index] -= particle_weight;
    }

    apply_transport_policy();
}

void self_consistent_device_mmmc_simulation_2d::add_missing_contact_charge_to_poisson_reservoir(
    std::size_t accumulation_steps) {
    const double accumulation_factor = static_cast<double>(accumulation_steps);
    for (std::size_t i = 0; i < m_list_element_contact_ptr.size(); ++i) {
        auto& element = m_list_element_contact_ptr[i];
        const double equilibrium_charge = m_list_element_contact_equilibrium_charge[i];
        const double accumulated_mobile_charge = element->get_n_charge() - element->get_p_charge();
        const double target_accumulated_charge = equilibrium_charge * accumulation_factor;
        const double correction = target_accumulated_charge - accumulated_mobile_charge;
        if (correction > 0.0) {
            element->add_n_charge(correction);
        } else if (correction < 0.0) {
            element->add_p_charge(-correction);
        }
    }
}

ADMC::vector3 self_consistent_device_mmmc_simulation_2d::draw_standard_normal() {
    return {m_standard_normal(m_admc_rng), m_standard_normal(m_admc_rng), m_standard_normal(m_admc_rng)};
}

mesh::vector3 self_consistent_device_mmmc_simulation_2d::to_mesh_position_um(const ADMC::vector3& position_m) const {
    mesh::vector3 position_um = position_m * units::meter_to_micron;
    position_um.to_2d_inplace();
    return position_um;
}

ADMC::vector3 self_consistent_device_mmmc_simulation_2d::to_admc_position_m(const mesh::vector3& position_um) const {
    mesh::vector3 normalized = position_um;
    normalized.to_2d_inplace();
    return normalized * units::micron_to_meter;
}

ADMC::carrier_type self_consistent_device_mmmc_simulation_2d::to_admc_type(PBMC::particle_type type) const {
    return to_admc_carrier_type(type);
}

PBMC::particle_type self_consistent_device_mmmc_simulation_2d::to_pbmc_type(ADMC::carrier_type type) const {
    return to_pbmc_particle_type(type);
}

PBMC::pbmc_transport_kernel& self_consistent_device_mmmc_simulation_2d::pbmc_transport_for(
    ADMC::carrier_type type) {
    return transport_for(to_pbmc_type(type));
}

ADMC::admc_local_environment self_consistent_device_mmmc_simulation_2d::local_admc_environment(
    const ADMC::device_admc_particle& particle) const {
    if (particle.containing_element == nullptr) {
        throw std::runtime_error("MMMC ADMC particle has no containing element.");
    }
    const mesh::vector3 position_um = to_mesh_position_um(particle.particle.state().position_m);
    const mesh::vector3 electric_field_V_per_cm =
        particle.containing_element->interpolate_electric_field_at_location(position_um);
    const double net_doping_cm_3 = local_net_doping_cm_3(*particle.containing_element, position_um);
    return {
        .electric_field_V_per_m    = electric_field_V_per_cm * units::electric_field_V_per_cm_to_V_per_m,
        .doping_concentration_cm_3 = std::abs(net_doping_cm_3),
        .lattice_temperature_K     = m_mmmc_options.m_admc.m_lattice_temperature_K,
    };
}

void self_consistent_device_mmmc_simulation_2d::advance_admc_particles_one_time_step(double dt_s) {
    for (auto& particle : m_admc_particles) {
        const auto environment = local_admc_environment(particle);
        m_admc_transport.step(particle.particle, environment, dt_s, draw_standard_normal());
        particle.particle.state().position_m.set_z(0.0);
        update_admc_element_and_check_boundary(particle);
    }
    remove_collected_admc_particles();
}

void self_consistent_device_mmmc_simulation_2d::update_admc_element_and_check_boundary(
    ADMC::device_admc_particle& particle) {
    if (particle.containing_element == nullptr) {
        particle.crossed_contact = true;
        return;
    }

    const mesh::vector3 current_position_um = to_mesh_position_um(particle.particle.state().position_m);
    const mesh::vector3 previous_position_um = to_mesh_position_um(particle.particle.state().previous_position_m);
    if (m_device.check_enters_contact(current_position_um) ||
        m_device.check_crossing_contact(previous_position_um, current_position_um)) {
        particle.crossed_contact = true;
        return;
    }
    if (particle.containing_element->is_location_inside_element(current_position_um)) {
        return;
    }

    auto* new_element = m_device.find_element_at_location(current_position_um);
    if (new_element == nullptr || !is_transport_material_element(*new_element)) {
        auto& state = particle.particle.state();
        const auto hit = mesh::find_boundary_exit_hit(*particle.containing_element,
                                                      previous_position_um,
                                                      current_position_um,
                                                      m_dimension);
        if (m_mmmc_options.m_admc.m_boundary_reflection_model == mesh::boundary_reflection_model::reverse ||
            !hit.has_value()) {
            state.position_m = state.previous_position_m;
            state.total_velocity_m_per_s *= -1.0;
            state.drift_velocity_m_per_s *= -1.0;
            return;
        }

        const mesh::vector3 remaining_displacement_um = current_position_um - hit->position;
        mesh::vector3 outgoing_displacement_um;
        if (m_mmmc_options.m_admc.m_boundary_reflection_model == mesh::boundary_reflection_model::specular) {
            state.total_velocity_m_per_s =
                mesh::reflect_vector_specular(state.total_velocity_m_per_s, hit->inward_normal);
            state.drift_velocity_m_per_s =
                mesh::reflect_vector_specular(state.drift_velocity_m_per_s, hit->inward_normal);
            outgoing_displacement_um = mesh::reflect_vector_specular(remaining_displacement_um, hit->inward_normal);
        } else {
            state.total_velocity_m_per_s = mesh::draw_diffuse_reflection_vector(state.total_velocity_m_per_s,
                                                                                hit->inward_normal,
                                                                                m_dimension,
                                                                                m_admc_rng);
            state.drift_velocity_m_per_s = mesh::draw_diffuse_reflection_vector(state.drift_velocity_m_per_s,
                                                                                hit->inward_normal,
                                                                                m_dimension,
                                                                                m_admc_rng);
            outgoing_displacement_um =
                mesh::align_displacement_with_direction(remaining_displacement_um, state.total_velocity_m_per_s);
        }
        state.position_m = to_admc_position_m(mesh::place_reflected_position_inside(*particle.containing_element,
                                                                                    previous_position_um,
                                                                                    current_position_um,
                                                                                    *hit,
                                                                                    outgoing_displacement_um,
                                                                                    m_dimension));
        return;
    }
    particle.containing_element = new_element;
}

void self_consistent_device_mmmc_simulation_2d::remove_collected_admc_particles() {
    std::erase_if(m_admc_particles, [](const auto& particle) { return particle.crossed_contact; });
}

void self_consistent_device_mmmc_simulation_2d::apply_transport_policy() {
    m_last_transfer_counters = {};

    for (std::size_t i = 0; i < m_list_particles.size();) {
        const auto& particle = *m_list_particles[i];
        if (m_self_consistent_options.m_policy.method_for_position(particle.state().position) ==
            transport_method::admc) {
            auto converted = convert_pbmc_to_admc(particle, particle.index());
            if (converted.containing_element != nullptr) {
                const auto environment = local_admc_environment(converted);
                converted.particle.state().electric_field_V_per_m = environment.electric_field_V_per_m;
                converted.particle.state().doping_concentration_cm_3 = environment.doping_concentration_cm_3;
                converted.particle.state().lattice_temperature_K = environment.lattice_temperature_K;
            }
            m_admc_particles.push_back(std::move(converted));
            m_list_particles[i] = std::move(m_list_particles.back());
            m_list_particles.pop_back();
            ++m_last_transfer_counters.pbmc_to_admc;
            continue;
        }
        ++i;
    }

    for (std::size_t i = 0; i < m_admc_particles.size();) {
        const auto& particle = m_admc_particles[i];
        if (m_self_consistent_options.m_policy.method_for_position(to_mesh_position_um(
                particle.particle.state().position_m)) == transport_method::pbmc) {
            auto& transport = pbmc_transport_for(particle.particle.type());
            auto converted = convert_admc_to_pbmc(particle, particle.particle.index(), transport);
            m_list_particles.push_back(std::make_unique<PBMC::pbmc_particle>(std::move(converted)));
            m_admc_particles[i] = std::move(m_admc_particles.back());
            m_admc_particles.pop_back();
            ++m_last_transfer_counters.admc_to_pbmc;
            continue;
        }
        ++i;
    }

    m_total_transfer_counters.pbmc_to_admc += m_last_transfer_counters.pbmc_to_admc;
    m_total_transfer_counters.admc_to_pbmc += m_last_transfer_counters.admc_to_pbmc;
}

void self_consistent_device_mmmc_simulation_2d::advance_mmmc_particles_one_time_step() {
    const double dt_s = m_simulation_options.m_time_step;
    transport_particles_one_time_step();
    advance_admc_particles_one_time_step(dt_s);
    apply_transport_policy();
}

std::pair<double, double> self_consistent_device_mmmc_simulation_2d::compute_admc_ramo_current() const {
    double electron_current = 0.0;
    double hole_current = 0.0;
    for (const auto& particle : m_admc_particles) {
        const double current = compute_admc_ramo_current_for_particle(particle);
        if (particle.particle.type() == ADMC::carrier_type::electron) {
            electron_current += current;
        } else {
            hole_current += current;
        }
    }
    return {electron_current, hole_current};
}

double self_consistent_device_mmmc_simulation_2d::compute_admc_ramo_current_for_particle(
    const ADMC::device_admc_particle& particle) const {
    const mesh::vector3 position_um = to_mesh_position_um(particle.particle.state().position_m);
    const mesh::vector3 ramo_field = get_RamoUnitaryElectricField_at_position(position_um, particle.containing_element);
    const double signed_charge_C = particle.weight * ADMC::carrier_charge_sign(particle.particle.type()) *
                                   constants::q_e;
    return -signed_charge_C * particle.particle.state().total_velocity_m_per_s.dot(ramo_field) *
           ramo_current_scale_factor();
}

double self_consistent_device_mmmc_simulation_2d::max_admc_particle_electric_field_V_per_cm() const {
    double max_field = 0.0;
    for (const auto& particle : m_admc_particles) {
        max_field = std::max(max_field,
                             particle.particle.state().electric_field_V_per_m.norm() * 0.01);
    }
    return max_field;
}

std::size_t self_consistent_device_mmmc_simulation_2d::get_total_number_electrons() const {
    return get_number_electrons() +
           static_cast<std::size_t>(std::count_if(m_admc_particles.begin(), m_admc_particles.end(), [](const auto& p) {
               return p.particle.type() == ADMC::carrier_type::electron;
           }));
}

std::size_t self_consistent_device_mmmc_simulation_2d::get_total_number_holes() const {
    return get_number_holes() +
           static_cast<std::size_t>(std::count_if(m_admc_particles.begin(), m_admc_particles.end(), [](const auto& p) {
               return p.particle.type() == ADMC::carrier_type::hole;
           }));
}

std::size_t self_consistent_device_mmmc_simulation_2d::get_number_pbmc_particles() const noexcept {
    return m_list_particles.size();
}

std::size_t self_consistent_device_mmmc_simulation_2d::get_number_admc_particles() const noexcept {
    return m_admc_particles.size();
}

std::size_t self_consistent_device_mmmc_simulation_2d::total_pbmc_to_admc_transfers() const noexcept {
    return m_total_transfer_counters.pbmc_to_admc;
}

std::size_t self_consistent_device_mmmc_simulation_2d::total_admc_to_pbmc_transfers() const noexcept {
    return m_total_transfer_counters.admc_to_pbmc;
}

void self_consistent_device_mmmc_simulation_2d::export_current_mmmc_snapshot() const {
    PBMC::device_pbmc_simulation::export_current_snapshot();
    const std::filesystem::path base_directory(m_simulation_options.m_prefix_export_filename);
    export_current_mmmc_particles_as_vtp((base_directory / "particles").string());
}

void self_consistent_device_mmmc_simulation_2d::export_mmmc_particle_state_csv(const std::string& filename) const {
    std::filesystem::path output_path(filename);
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    std::ofstream stream(filename);
    if (!stream.is_open()) {
        throw std::runtime_error(fmt::format("Could not open MMMC particle CSV file '{}'", filename));
    }
    stream << std::setprecision(std::numeric_limits<double>::max_digits10);
    stream << "particle_index,carrier_type_code,carrier_type_name,transport_method_code,transport_method_name,"
              "has_pbmc_state,has_admc_state,time_s,x_um,y_um,z_um,vx_m_per_s,vy_m_per_s,vz_m_per_s,"
              "drift_vx_m_per_s,drift_vy_m_per_s,drift_vz_m_per_s,local_kx_1_per_m,local_ky_1_per_m,"
              "local_kz_1_per_m,energy_eV,gamma_eV,valley_index,mobility_m2_per_V_s,diffusion_m2_per_s,"
              "doping_concentration_cm_3,lattice_temperature_K,electric_field_x_V_per_m,"
              "electric_field_y_V_per_m,electric_field_z_V_per_m,electric_field_norm_V_per_m,"
              "electric_field_x_V_per_cm,electric_field_y_V_per_cm,electric_field_z_V_per_cm,"
              "electric_field_norm_V_per_cm,weight,signed_charge_C\n";

    for (const auto& p_particle : m_list_particles) {
        const auto& particle = *p_particle;
        const auto& state = particle.state();
        const auto electric_field_V_per_m = state.electric_field * units::electric_field_V_per_cm_to_V_per_m;
        stream << particle.index() << ',' << static_cast<int>(particle.type()) << ','
               << PBMC::carrier_type_to_string(particle.type()) << ",0,PBMC,1,0," << state.time << ','
               << state.position.x() << ',' << state.position.y() << ',' << state.position.z() << ','
               << state.velocity.x() << ',' << state.velocity.y() << ',' << state.velocity.z()
               << ",0,0,0," << state.local_k.x() << ',' << state.local_k.y() << ',' << state.local_k.z() << ','
               << state.kinetic_energy << ',' << state.gamma << ',' << static_cast<long long>(state.valley_index)
               << ",0,0," << state.doping_concentration_cm_3 << ',' << state.lattice_temperature_K << ','
               << electric_field_V_per_m.x() << ',' << electric_field_V_per_m.y() << ','
               << electric_field_V_per_m.z() << ',' << electric_field_V_per_m.norm() << ','
               << state.electric_field.x() << ',' << state.electric_field.y() << ',' << state.electric_field.z()
               << ',' << state.electric_field.norm() << ',' << particle.weight() << ','
               << particle.get_signed_charge() << '\n';
    }

    for (const auto& particle : m_admc_particles) {
        const auto& state = particle.particle.state();
        const auto position_um = to_mesh_position_um(state.position_m);
        const auto electric_field_V_per_cm = state.electric_field_V_per_m * 0.01;
        stream << particle.particle.index() << ',' << static_cast<int>(to_pbmc_type(particle.particle.type()))
               << ',' << ADMC::carrier_type_name(particle.particle.type()) << ",1,ADMC,0,1," << state.time_s
               << ',' << position_um.x() << ',' << position_um.y() << ',' << position_um.z() << ','
               << state.total_velocity_m_per_s.x() << ',' << state.total_velocity_m_per_s.y() << ','
               << state.total_velocity_m_per_s.z() << ',' << state.drift_velocity_m_per_s.x() << ','
               << state.drift_velocity_m_per_s.y() << ',' << state.drift_velocity_m_per_s.z()
               << ",0,0,0,0,0,-1," << state.mobility_m2_per_V_s << ',' << state.diffusion_m2_per_s << ','
               << state.doping_concentration_cm_3 << ',' << state.lattice_temperature_K << ','
               << state.electric_field_V_per_m.x() << ',' << state.electric_field_V_per_m.y() << ','
               << state.electric_field_V_per_m.z() << ',' << state.electric_field_V_per_m.norm() << ','
               << electric_field_V_per_cm.x() << ',' << electric_field_V_per_cm.y() << ','
               << electric_field_V_per_cm.z() << ',' << electric_field_V_per_cm.norm() << ',' << particle.weight
               << ',' << admc_signed_charge_C(particle) << '\n';
    }
}

void self_consistent_device_mmmc_simulation_2d::export_current_mmmc_particles_as_vtp(
    const std::string& directory) const {
    const std::filesystem::path output_directory(directory);
    std::filesystem::create_directories(output_directory);

    const std::string filename = fmt::format("particles_{:012d}.vtp", m_state.m_iteration);
    const std::filesystem::path vtp_path = output_directory / filename;
    const std::filesystem::path pvd_path = output_directory / "particles.pvd";
    std::ofstream stream(vtp_path);
    if (!stream.is_open()) {
        throw std::runtime_error(fmt::format("Could not open MMMC particle VTP file '{}'", vtp_path.string()));
    }

    const std::size_t number_particles = m_list_particles.size() + m_admc_particles.size();
    stream << std::setprecision(std::numeric_limits<double>::max_digits10);
    stream << "<?xml version=\"1.0\"?>\n";
    stream << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    stream << "  <PolyData>\n";
    stream << "    <Piece NumberOfPoints=\"" << number_particles << "\" NumberOfVerts=\"" << number_particles
           << "\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
    stream << "      <PointData Scalars=\"energy_eV\" Vectors=\"velocity_m_per_s\">\n";

    stream << "        <DataArray type=\"Int64\" Name=\"particle_index\" format=\"ascii\">\n          ";
    for (const auto& p_particle : m_list_particles) {
        stream << static_cast<long long>(p_particle->index()) << ' ';
    }
    for (const auto& particle : m_admc_particles) {
        stream << static_cast<long long>(particle.particle.index()) << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Int32\" Name=\"particle_type\" format=\"ascii\">\n          ";
    for (const auto& p_particle : m_list_particles) {
        stream << static_cast<int>(p_particle->type()) << ' ';
    }
    for (const auto& particle : m_admc_particles) {
        stream << static_cast<int>(to_pbmc_type(particle.particle.type())) << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Int32\" Name=\"transport_method\" format=\"ascii\">\n          ";
    for (std::size_t i = 0; i < m_list_particles.size(); ++i) {
        stream << 0 << ' ';
    }
    for (std::size_t i = 0; i < m_admc_particles.size(); ++i) {
        stream << 1 << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Int32\" Name=\"has_pbmc_state\" format=\"ascii\">\n          ";
    for (std::size_t i = 0; i < m_list_particles.size(); ++i) {
        stream << 1 << ' ';
    }
    for (std::size_t i = 0; i < m_admc_particles.size(); ++i) {
        stream << 0 << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Int32\" Name=\"has_admc_state\" format=\"ascii\">\n          ";
    for (std::size_t i = 0; i < m_list_particles.size(); ++i) {
        stream << 0 << ' ';
    }
    for (std::size_t i = 0; i < m_admc_particles.size(); ++i) {
        stream << 1 << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"time_s\" format=\"ascii\">\n          ";
    for (const auto& p_particle : m_list_particles) {
        stream << p_particle->state().time << ' ';
    }
    for (const auto& particle : m_admc_particles) {
        stream << particle.particle.state().time_s << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"energy_eV\" format=\"ascii\">\n          ";
    for (const auto& p_particle : m_list_particles) {
        stream << p_particle->state().kinetic_energy << ' ';
    }
    for (std::size_t i = 0; i < m_admc_particles.size(); ++i) {
        stream << 0.0 << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"lattice_temperature_K\" format=\"ascii\">\n          ";
    for (const auto& p_particle : m_list_particles) {
        stream << p_particle->state().lattice_temperature_K << ' ';
    }
    for (const auto& particle : m_admc_particles) {
        stream << particle.particle.state().lattice_temperature_K << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"weight\" format=\"ascii\">\n          ";
    for (const auto& p_particle : m_list_particles) {
        stream << p_particle->weight() << ' ';
    }
    for (const auto& particle : m_admc_particles) {
        stream << particle.weight << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"mobility_m2_per_V_s\" format=\"ascii\">\n          ";
    for (std::size_t i = 0; i < m_list_particles.size(); ++i) {
        stream << 0.0 << ' ';
    }
    for (const auto& particle : m_admc_particles) {
        stream << particle.particle.state().mobility_m2_per_V_s << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"diffusion_m2_per_s\" format=\"ascii\">\n          ";
    for (std::size_t i = 0; i < m_list_particles.size(); ++i) {
        stream << 0.0 << ' ';
    }
    for (const auto& particle : m_admc_particles) {
        stream << particle.particle.state().diffusion_m2_per_s << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"doping_concentration_cm_3\" format=\"ascii\">\n          ";
    for (const auto& p_particle : m_list_particles) {
        stream << p_particle->state().doping_concentration_cm_3 << ' ';
    }
    for (const auto& particle : m_admc_particles) {
        stream << particle.particle.state().doping_concentration_cm_3 << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Int32\" Name=\"valley_index\" format=\"ascii\">\n          ";
    for (const auto& p_particle : m_list_particles) {
        stream << static_cast<int>(p_particle->state().valley_index) << ' ';
    }
    for (std::size_t i = 0; i < m_admc_particles.size(); ++i) {
        stream << -1 << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"signed_charge_C\" format=\"ascii\">\n          ";
    for (const auto& p_particle : m_list_particles) {
        stream << p_particle->get_signed_charge() << ' ';
    }
    for (const auto& particle : m_admc_particles) {
        stream << admc_signed_charge_C(particle) << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"electric_field_norm_V_per_cm\" format=\"ascii\">\n          ";
    for (const auto& p_particle : m_list_particles) {
        stream << p_particle->state().electric_field.norm() << ' ';
    }
    for (const auto& particle : m_admc_particles) {
        stream << particle.particle.state().electric_field_V_per_m.norm() * 0.01 << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"velocity_m_per_s\" NumberOfComponents=\"3\" "
              "format=\"ascii\">\n          ";
    for (const auto& p_particle : m_list_particles) {
        const auto& velocity = p_particle->state().velocity;
        stream << velocity.x() << ' ' << velocity.y() << ' ' << velocity.z() << ' ';
    }
    for (const auto& particle : m_admc_particles) {
        const auto& velocity = particle.particle.state().total_velocity_m_per_s;
        stream << velocity.x() << ' ' << velocity.y() << ' ' << velocity.z() << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"drift_velocity_m_per_s\" NumberOfComponents=\"3\" "
              "format=\"ascii\">\n          ";
    for (std::size_t i = 0; i < m_list_particles.size(); ++i) {
        stream << "0 0 0 ";
    }
    for (const auto& particle : m_admc_particles) {
        const auto& velocity = particle.particle.state().drift_velocity_m_per_s;
        stream << velocity.x() << ' ' << velocity.y() << ' ' << velocity.z() << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"local_k_1_per_m\" NumberOfComponents=\"3\" "
              "format=\"ascii\">\n          ";
    for (const auto& p_particle : m_list_particles) {
        const auto& local_k = p_particle->state().local_k;
        stream << local_k.x() << ' ' << local_k.y() << ' ' << local_k.z() << ' ';
    }
    for (std::size_t i = 0; i < m_admc_particles.size(); ++i) {
        stream << "0 0 0 ";
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"electric_field_V_per_m\" NumberOfComponents=\"3\" "
              "format=\"ascii\">\n          ";
    for (const auto& p_particle : m_list_particles) {
        const auto electric_field = p_particle->state().electric_field *
                                    units::electric_field_V_per_cm_to_V_per_m;
        stream << electric_field.x() << ' ' << electric_field.y() << ' ' << electric_field.z() << ' ';
    }
    for (const auto& particle : m_admc_particles) {
        const auto& electric_field = particle.particle.state().electric_field_V_per_m;
        stream << electric_field.x() << ' ' << electric_field.y() << ' ' << electric_field.z() << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"electric_field_V_per_cm\" NumberOfComponents=\"3\" "
              "format=\"ascii\">\n          ";
    for (const auto& p_particle : m_list_particles) {
        const auto& electric_field = p_particle->state().electric_field;
        stream << electric_field.x() << ' ' << electric_field.y() << ' ' << electric_field.z() << ' ';
    }
    for (const auto& particle : m_admc_particles) {
        const auto electric_field = particle.particle.state().electric_field_V_per_m * 0.01;
        stream << electric_field.x() << ' ' << electric_field.y() << ' ' << electric_field.z() << ' ';
    }
    stream << "\n        </DataArray>\n";

    stream << "      </PointData>\n";
    stream << "      <Points>\n";
    stream << "        <DataArray type=\"Float64\" Name=\"position_um\" NumberOfComponents=\"3\" format=\"ascii\">\n"
              "          ";
    for (const auto& p_particle : m_list_particles) {
        const auto& position = p_particle->state().position;
        stream << position.x() << ' ' << position.y() << ' ' << position.z() << ' ';
    }
    for (const auto& particle : m_admc_particles) {
        const auto position = to_mesh_position_um(particle.particle.state().position_m);
        stream << position.x() << ' ' << position.y() << ' ' << position.z() << ' ';
    }
    stream << "\n        </DataArray>\n";
    stream << "      </Points>\n";
    stream << "      <Verts>\n";
    stream << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n          ";
    for (std::size_t i = 0; i < number_particles; ++i) {
        stream << static_cast<long long>(i) << ' ';
    }
    stream << "\n        </DataArray>\n";
    stream << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n          ";
    for (std::size_t i = 0; i < number_particles; ++i) {
        stream << static_cast<long long>(i + 1) << ' ';
    }
    stream << "\n        </DataArray>\n";
    stream << "      </Verts>\n";
    stream << "    </Piece>\n";
    stream << "  </PolyData>\n";
    stream << "</VTKFile>\n";

    const auto already_recorded =
        std::find_if(m_particle_vtp_export_records.begin(),
                     m_particle_vtp_export_records.end(),
                     [&](const PBMC::vtk_time_series_record& record) { return record.m_filename == filename; });
    if (already_recorded == m_particle_vtp_export_records.end()) {
        m_particle_vtp_export_records.push_back(
            PBMC::vtk_time_series_record{.m_time_s = m_state.m_time_s, .m_filename = filename});
    }
    write_particle_vtp_time_collection(pvd_path.string());
}

void self_consistent_device_mmmc_simulation_2d::run_self_consistent_transport_simulation() {
    const std::size_t total_iterations =
        static_cast<std::size_t>(std::ceil(m_simulation_options.m_t_max / m_simulation_options.m_time_step));
    const std::size_t poisson_frequency_value = poisson_frequency();
    const double      sim_poisson_frequency = static_cast<double>(poisson_frequency_value);

    const std::string history_filename = initialize_simulation_history_file();
    std::fstream history_stream(history_filename, std::ios::app);
    history_stream << std::setprecision(std::numeric_limits<double>::max_digits10);
    if (!history_stream.is_open()) {
        throw std::runtime_error("Could not open MMMC simulation history CSV file '" + history_filename + "'.");
    }

    fmt::print("START 2D SELF-CONSISTENT MMMC SIMULATION\n");
    fmt::print("Total iterations: {}\n", total_iterations);
    fmt::print("Poisson frequency: {}\n", poisson_frequency_value);
    fmt::print("PBMC bbox: x=[{:.6e},{:.6e}] y=[{:.6e},{:.6e}] z=[{:.6e},{:.6e}] um\n",
               m_self_consistent_options.m_policy.m_pbmc_region_um.get_x_min(),
               m_self_consistent_options.m_policy.m_pbmc_region_um.get_x_max(),
               m_self_consistent_options.m_policy.m_pbmc_region_um.get_y_min(),
               m_self_consistent_options.m_policy.m_pbmc_region_um.get_y_max(),
               m_self_consistent_options.m_policy.m_pbmc_region_um.get_z_min(),
               m_self_consistent_options.m_policy.m_pbmc_region_um.get_z_max());

    double accumulator_ramo_current_electron = 0.0;
    double accumulator_ramo_current_hole = 0.0;
    double ramo_current_electron = 0.0;
    double ramo_current_hole = 0.0;
    double ramo_current = 0.0;

    reset_element_charges();
    while (m_state.m_time_s <= m_simulation_options.m_t_max) {
        if (m_state.m_iteration > 10 && m_simulation_options.m_stop_simu_when_no_electron_remaining &&
            get_total_number_electrons() == 0 && !has_pending_scheduled_particle_injection()) {
            fmt::print("Stop: no electrons remaining in device.\n");
            break;
        }
        if (m_list_particles.size() + m_admc_particles.size() >= m_simulation_options.m_max_number_particle) {
            fmt::print("Stop: hard particle limit reached.\n");
            break;
        }

        const std::size_t impact_events_before_step = m_simulation_history.m_impact_ionization_positions.size();
        advance_mmmc_particles_one_time_step();
        add_particle_charges_to_elements();

        const auto pbmc_currents = compute_ramo_current();
        const auto admc_currents = compute_admc_ramo_current();
        accumulator_ramo_current_electron += pbmc_currents.first + admc_currents.first;
        accumulator_ramo_current_hole += pbmc_currents.second + admc_currents.second;

        const bool should_update_poisson =
            (m_state.m_iteration % poisson_frequency_value == 0) && (m_state.m_iteration != 0);
        if (should_update_poisson) {
            const double poisson_sample_time_s = m_state.m_time_s + m_simulation_options.m_time_step;
            apply_scheduled_contact_voltage_events(poisson_sample_time_s);

            ramo_current_electron = accumulator_ramo_current_electron / sim_poisson_frequency;
            ramo_current_hole = accumulator_ramo_current_hole / sim_poisson_frequency;
            ramo_current = ramo_current_electron + ramo_current_hole - common_options().m_background_ramo_current_A;
            accumulator_ramo_current_electron = 0.0;
            accumulator_ramo_current_hole = 0.0;

            add_charges_at_contacts(poisson_frequency_value);
            add_particle_charges_to_elements();
            add_missing_contact_charge_to_poisson_reservoir(poisson_frequency_value + 1);
            recompute_vertex_space_charge_from_element_charges(poisson_frequency_value + 1);
            update_self_consistent_potential(false);
            reset_element_charges();
        }

        m_state.m_time_s += m_simulation_options.m_time_step;
        ++m_state.m_iteration;
        update_successful_quench_detection(m_state.m_time_s, impact_events_before_step);

        const double max_electric_field_V_per_cm =
            std::max(max_particle_electric_field_V_per_cm(), max_admc_particle_electric_field_V_per_cm());
        m_simulation_history.add_data_to_history(m_state.m_time_s,
                                                 get_total_number_electrons(),
                                                 get_total_number_holes(),
                                                 m_simulation_history.m_impact_ionization_positions.size(),
                                                 ramo_current_electron,
                                                 ramo_current_hole,
                                                 ramo_current,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 max_electric_field_V_per_cm,
                                                 ramo_electrode_voltage_for_history(),
                                                 reference_electrode_voltage_for_history(),
                                                 quench_supply_voltage_for_history(),
                                                 quench_device_current_for_history(),
                                                 quench_resistor_current_for_history(),
                                                 quench_voltage_drop_for_history(),
                                                 m_simulation_history.contact_voltage_values_from_map(
                                                     contact_voltages_V()));

        if (m_state.m_iteration == 1 || m_state.m_iteration % 10 == 0) {
            m_simulation_history.append_last_iter_to_csv(history_stream);
        }

        if (m_state.m_iteration == 1 ||
            m_state.m_iteration % static_cast<std::size_t>(m_simulation_options.m_frequency_export_trajectory) == 0) {
            fmt::print("\rMMMC t={:<10.3e} ps - {:>9d} / {} ({:.1f}%) -- PBMC: {} ADMC: {} -- "
                       "transfers P->A/A->P: {}/{} -- Ramo current: {:.6e} A",
                       m_state.m_time_s * 1e12,
                       m_state.m_iteration,
                       total_iterations,
                       total_iterations == 0
                           ? 100.0
                           : static_cast<double>(m_state.m_iteration) / static_cast<double>(total_iterations) * 100.0,
                       m_list_particles.size(),
                       m_admc_particles.size(),
                       m_total_transfer_counters.pbmc_to_admc,
                       m_total_transfer_counters.admc_to_pbmc,
                       ramo_current);
            std::fflush(stdout);
            history_stream.flush();
            if (m_simulation_options.m_export_time_step) {
                m_poisson_solver.add_solution_to_mesh_functions("PoissonSolution", true);
                export_current_mmmc_snapshot();
            }
        }
    }

    history_stream.close();
    m_poisson_solver.add_solution_to_mesh_functions("PoissonSolution", true);
    export_current_mmmc_snapshot();
    fmt::print("\nEND 2D SELF-CONSISTENT MMMC SIMULATION\n");
}

}  // namespace uepm::MMMC
