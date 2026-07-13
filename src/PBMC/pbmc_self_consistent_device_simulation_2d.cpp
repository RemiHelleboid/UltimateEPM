/**
 * @file pbmc_self_consistent_device_simulation_2d.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-27
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "pbmc_self_consistent_device_simulation_2d.hpp"

#include <fmt/core.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "unit_conversion.hpp"
#include "physical_constants.hpp"

namespace uepm::PBMC {

void options_self_consistent_device_pbmc_2d::validate() const {
    if (m_effective_depth_um <= 0.0) {
        throw std::invalid_argument("--effective-depth must be positive.");
    }
}

double self_consistent_device_pbmc_simulation_2d::scale_integrated_2d_doping_to_carriers(
    double integrated_doping) const {
    constexpr double micron_to_cm = 1.0e-4;
    return integrated_doping * m_self_consistent_options.m_effective_depth_um * micron_to_cm;
}

double self_consistent_device_pbmc_simulation_2d::charge_deposition_factor(std::size_t accumulation_steps) const {
    if (accumulation_steps == 0) {
        throw std::invalid_argument("accumulation_steps must be positive.");
    }
    return 1.0 / (static_cast<double>(accumulation_steps) * m_self_consistent_options.m_effective_depth_um);
}

void self_consistent_device_pbmc_simulation_2d::validate_self_consistent_options() const {
    if (m_dimension != 2) {
        throw std::invalid_argument("self_consistent_device_pbmc_simulation_2d requires a 2D device.");
    }
    if (poisson_frequency() == 0) {
        throw std::invalid_argument("Poisson frequency must be positive.");
    }
    m_self_consistent_options.validate();
}

double self_consistent_device_pbmc_simulation_2d::ramo_current_scale_factor() const {
    return 1.0 / m_self_consistent_options.m_effective_depth_um;
}

double self_consistent_device_pbmc_simulation_2d::current_density_cell_volume_m3(const mesh::element& element) const {
    const double area_um2             = std::abs(element.get_measure());
    const double effective_volume_um3 = area_um2 * m_self_consistent_options.m_effective_depth_um;
    return effective_volume_um3 * std::pow(uepm::units::micron_to_meter, 3);
}

/**
 * @brief We need to identify the elements adjacent to the contacts, and compute their equilibrium charge (integral of
 * doping concentration) to be able to properly update their charge during the simulation (so that they are not fixed at
 * the equilibrium charge, but can exchange charge with the device while still having a reference equilibrium charge for
 * the poisson solver).
 *
 */
void self_consistent_device_pbmc_simulation_2d::initialize_contact_elements() {
    const auto list_bulk_elements = m_device.get_p_mesh()->get_list_bulk_element();

    m_list_element_contact.clear();
    m_list_element_contact_ptr.clear();
    m_list_element_contact_equilibrium_charge.clear();

    for (const auto& device_contact : m_device.get_list_contacts()) {
        const std::string contact_name = device_contact.get_contact_name();
        if (!contact_voltages_V().contains(contact_name)) {
            throw std::runtime_error("Particle-collection contact '" + contact_name +
                                     "' has no configured Poisson voltage.");
        }
        const auto element_indices =
            m_device.get_p_mesh()->get_idx_bulk_elements_adjacent_to_contact_region(contact_name);
        std::vector<std::shared_ptr<mesh::element>> contact_elements;
        contact_elements.reserve(element_indices.size());

        for (const auto element_index : element_indices) {
            if (element_index >= list_bulk_elements.size()) {
                throw std::runtime_error("Contact-adjacent bulk element index is out of range.");
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

void self_consistent_device_pbmc_simulation_2d::place_initial_charges_according_to_doping(double particle_weight) {
    if (particle_weight <= 0.0) {
        throw std::invalid_argument("particle_weight must be positive.");
    }

    const std::string donor_field_name                = "DonorConcentration";
    const std::string acceptor_field_name             = "AcceptorConcentration";
    auto*             mesh                            = m_device.get_p_mesh();
    const auto        integrate_carriers_over_2d_mesh = [&](const std::string& field_name) {
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

    // Pour éviter des trucs relou bref
    constexpr double min_probability = 50e-2;

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

    const double max_donor_concentration    = mesh->get_argmax_max_of_function(donor_field_name).second;
    const double max_acceptor_concentration = mesh->get_argmax_max_of_function(acceptor_field_name).second;

    if (number_electrons > 0 && max_donor_concentration <= 0.0) {
        throw std::runtime_error("Donor concentration maximum is non-positive.");
    }

    if (number_holes > 0 && max_acceptor_concentration <= 0.0) {
        throw std::runtime_error("Acceptor concentration maximum is non-positive.");
    }
    const mesh::bbox  active_region_bbox = mesh->get_p_region("Silicon_1")->compute_bounding_box();
    const std::size_t max_trial          = 100 * number_electrons;
    std::size_t       count_trial        = 0;

    while (electron_positions.size() < number_electrons) {
        const mesh::vector3 position = active_region_bbox.draw_uniform_random_point_inside_box(m_contact_rng);
        count_trial++;
        const double donor_density    = mesh->interpolate_scalar_at_location(donor_field_name, position);
        const double acceptor_density = mesh->interpolate_scalar_at_location(acceptor_field_name, position);
        if (acceptor_density > donor_density) {
            continue;
        }

        const double probability = donor_density / max_donor_concentration;

        if (probability > min_probability && uniform01(m_contact_rng) < probability) {
            electron_positions.push_back(position);
        }
        if (count_trial > max_trial) {
            fmt::print("Error\n");
        }
    }
    fmt::print("Initial electrons placed according to doping.\n");

    const std::size_t max_trial_holes = 100 * number_holes;
    count_trial                       = 0;

    while (hole_positions.size() < number_holes) {
        const mesh::vector3 position         = active_region_bbox.draw_uniform_random_point_inside_box(m_contact_rng);
        const double        acceptor_density = mesh->interpolate_scalar_at_location(acceptor_field_name, position);
        const double        donor_density    = mesh->interpolate_scalar_at_location(donor_field_name, position);
        count_trial++;
        if (donor_density > acceptor_density) {
            continue;
        }

        const double probability = acceptor_density / max_acceptor_concentration;
        if (probability > min_probability && uniform01(m_contact_rng) < probability) {
            hole_positions.push_back(position);
        }
        if (count_trial > max_trial_holes) {
            fmt::print("Could not place all initial holes according to doping. Placed {} out of {}.\n",
                       hole_positions.size(),
                       number_holes);
            break;
        }
    }
    fmt::print("Initial holes placed according to doping.\n");

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

void self_consistent_device_pbmc_simulation_2d::add_charges_at_contacts(std::size_t poisson_frequency) {
    if (poisson_frequency == 0) {
        throw std::invalid_argument("Poisson frequency must be positive.");
    }
    const double particle_weight = common_options().m_contact_injection_particle_weight;
    if (particle_weight <= 0.0) {
        throw std::invalid_argument("Contact injection particle weight must be positive.");
    }
    const std::size_t number_contact_elements = m_list_element_contact_ptr.size();
    if (number_contact_elements == 0) {
        return;
    }

    std::vector<double> electron_charge_to_add(number_contact_elements, 0.0);
    std::vector<double> hole_charge_to_add(number_contact_elements, 0.0);

    for (std::size_t i = 0; i < number_contact_elements; ++i) {
        auto& element = m_list_element_contact_ptr[i];

        const double element_charge          = element->get_n_charge() - element->get_p_charge();
        const double averaged_element_charge = element_charge / static_cast<double>(poisson_frequency);
        const double equilibrium_charge      = m_list_element_contact_equilibrium_charge[i];
        const double charge_to_add           = equilibrium_charge - averaged_element_charge;

        if (equilibrium_charge > 0.0 && charge_to_add > 0.0) {
            electron_charge_to_add[i] = charge_to_add;
        } else if (equilibrium_charge < 0.0 && charge_to_add < 0.0) {
            hole_charge_to_add[i] = -charge_to_add;
        }
    }

    std::vector<mesh::vector3> electron_positions;
    std::vector<mesh::vector3> hole_positions;
    std::uniform_real_distribution<double> uniform01(0.0, 1.0);

    const auto remaining_capacity = [&]() {
        const std::size_t reserved = m_list_particles.size() + electron_positions.size() + hole_positions.size();
        return reserved < m_simulation_options.m_max_number_particle
                   ? m_simulation_options.m_max_number_particle - reserved
                   : std::size_t{0};
    };

    const auto allocate_positions_by_local_deficit = [&](const std::vector<double>& local_deficits,
                                                          std::vector<mesh::vector3>& positions) {
        for (std::size_t i = 0; i < local_deficits.size(); ++i) {
            const double exact_particle_count = local_deficits[i] / particle_weight;
            if (!(exact_particle_count > 0.0) || !std::isfinite(exact_particle_count)) {
                continue;
            }

            const double      integral_count_as_double = std::floor(exact_particle_count);
            const std::size_t capacity                 = remaining_capacity();
            const std::size_t integral_count = integral_count_as_double >= static_cast<double>(capacity)
                                                   ? capacity
                                                   : static_cast<std::size_t>(integral_count_as_double);

            positions.reserve(positions.size() + integral_count + 1);
            for (std::size_t particle_index = 0; particle_index < integral_count; ++particle_index) {
                positions.push_back(
                    m_list_element_contact_ptr[i]->draw_uniform_random_point_inside_element(m_contact_rng));
            }

            if (remaining_capacity() == 0 || integral_count_as_double > static_cast<double>(integral_count)) {
                continue;
            }
            const double fractional_count = exact_particle_count - integral_count_as_double;
            if (fractional_count > 0.0 && uniform01(m_contact_rng) < fractional_count) {
                positions.push_back(
                    m_list_element_contact_ptr[i]->draw_uniform_random_point_inside_element(m_contact_rng));
            }
        }
    };

    allocate_positions_by_local_deficit(electron_charge_to_add, electron_positions);
    allocate_positions_by_local_deficit(hole_charge_to_add, hole_positions);

    add_particles_at_positions(electron_positions, particle_type::electron, particle_weight);
    add_particles_at_positions(hole_positions, particle_type::hole, particle_weight);
}

void self_consistent_device_pbmc_simulation_2d::add_missing_contact_charge_to_poisson_reservoir(
    std::size_t accumulation_steps) {
    if (accumulation_steps == 0) {
        throw std::invalid_argument("accumulation_steps must be positive.");
    }
    const double accumulation_factor = static_cast<double>(accumulation_steps);
    for (std::size_t i = 0; i < m_list_element_contact_ptr.size(); ++i) {
        auto&        element                   = m_list_element_contact_ptr[i];
        const double equilibrium_charge        = m_list_element_contact_equilibrium_charge[i];
        const double accumulated_mobile_charge = element->get_n_charge() - element->get_p_charge();
        const double target_accumulated_charge = equilibrium_charge * accumulation_factor;
        const double correction                = target_accumulated_charge - accumulated_mobile_charge;
        if (correction > 0.0) {
            element->add_n_charge(correction);
        } else if (correction < 0.0) {
            element->add_p_charge(-correction);
        }
    }
}

void self_consistent_device_pbmc_simulation_2d::compute_unitary_potential() {
    m_poisson_solver.compute_second_member(0.0);
    for (const auto& [contact_name, unused_voltage] : contact_voltages_V()) {
        static_cast<void>(unused_voltage);
        m_poisson_solver.apply_dirichlet_condition(contact_name, contact_name == ramo_electrode() ? 1.0 : 0.0);
    }
    m_poisson_solver.decompose_matrix();
    m_poisson_solver.solve_system();

    constexpr bool add_gradient = true;
    m_poisson_solver.add_solution_to_mesh_functions("RamoUnitaryPotential", add_gradient);

    // Check if the Unitary is constant everywhere.
    const double unitary_potential_max =
        m_device.get_p_mesh()->get_argmax_max_of_function("RamoUnitaryPotential_gradient_norm").second;
    const double unitary_potential_min =
        m_device.get_p_mesh()->get_argmin_min_of_function("RamoUnitaryPotential_gradient_norm").second;
    if (std::abs(unitary_potential_max - unitary_potential_min) < 1e-6) {
        fmt::print("INFO : Unitary potential is constant across the device (max = {:.6e}, min = {:.6e}). The Constant "
                   "Unitary Electric Field will be used.\n",
                   unitary_potential_max,
                   unitary_potential_min);
        m_state.m_use_constant_RamoUnitaryElectricField = true;
        m_state.m_RamoUnitaryElectricField_Vm_per_cm =
            m_device.get_p_mesh()->interpolate_vector_at_location("RamoUnitaryPotential_gradient", {1e-3, 1e-3, 0.0});
    } else {
        fmt::print("Unitary potential computed. Max value: {:.6e}, Min value: {:.6e}\n",
                   unitary_potential_max,
                   unitary_potential_min);
    }
}

void self_consistent_device_pbmc_simulation_2d::initialize_poisson_solver() {
    m_poisson_solver.compute_stiffness_matrix();
    compute_unitary_potential();
    const double unitary_potential_max =
        m_device.get_p_mesh()->get_argmax_max_of_function("RamoUnitaryPotential").second;
    if (unitary_potential_max <= 0.0) {
        throw std::runtime_error("Unitary potential maximum is non-positive.");
    }
    fmt::print("Unitary potential computed. Max value: {:.6e}\n", unitary_potential_max);

    // The weighting-potential solve modifies the matrix. Reassemble it for the
    // physical Poisson problem and constrain every configured contact.
    m_poisson_solver.compute_stiffness_matrix();
    m_unconstrained_poisson_stiffness = m_poisson_solver.get_lhs_matrix();
    m_poisson_solver.compute_second_member(0.0);
    for (const auto& [contact_name, unused_voltage] : contact_voltages_V()) {
        static_cast<void>(unused_voltage);
        m_poisson_solver.apply_dirichlet_condition(contact_name, contact_voltage_for_poisson(contact_name));
    }
    m_poisson_solver.decompose_matrix();
}

self_consistent_device_pbmc_simulation_2d::self_consistent_device_pbmc_simulation_2d(
    const device::device&                         simulation_device,
    const options_device_PBMC&                    simulation_options,
    const options_self_consistent_device_pbmc_2d& self_consistent_options,
    const physics::material_database&             material_database,
    int                                           seed_random_generator)
    : self_consistent_device_pbmc_simulation_base(simulation_device,
                                                  simulation_options,
                                                  self_consistent_options.m_common,
                                                  seed_random_generator),
      m_self_consistent_options(self_consistent_options),
      m_poisson_solver(m_device.get_p_mesh(), m_device.get_p_mesh()->get_nb_vertices(), material_database),
      m_contact_rng(seed_random_generator + 1) {
    validate_self_consistent_options();
    apply_scheduled_contact_voltage_events(0.0);
    initialize_contact_elements();
    initialize_poisson_solver();
    initialize_particles_for_self_consistent_run();
}

self_consistent_device_pbmc_simulation_2d::self_consistent_device_pbmc_simulation_2d(
    const device::device&                         simulation_device,
    const options_device_PBMC&                    simulation_options,
    const options_self_consistent_device_pbmc_2d& self_consistent_options,
    const physics::material_database&             material_database,
    const mesh::vector3&                          starting_position,
    std::size_t                                   number_electrons_start,
    std::size_t                                   number_holes_start,
    int                                           seed_random_generator)
    : self_consistent_device_pbmc_simulation_base(simulation_device,
                                                  simulation_options,
                                                  self_consistent_options.m_common,
                                                  starting_position,
                                                  number_electrons_start,
                                                  number_holes_start,
                                                  seed_random_generator),
      m_self_consistent_options(self_consistent_options),
      m_poisson_solver(m_device.get_p_mesh(), m_device.get_p_mesh()->get_nb_vertices(), material_database),
      m_contact_rng(seed_random_generator + 1) {
    validate_self_consistent_options();
    apply_scheduled_contact_voltage_events(0.0);
    initialize_contact_elements();
    initialize_poisson_solver();
    initialize_particles_for_self_consistent_run();
}

void self_consistent_device_pbmc_simulation_2d::initialize_particles_for_self_consistent_run() {
    if (!common_options().m_initial_particle_state_file.empty()) {
        load_particles_from_state_csv(common_options().m_initial_particle_state_file);
        flatten_particle_positions_for_2d();
        reset_element_charges();
        add_particle_charges_to_elements();
        recompute_vertex_space_charge_from_element_charges(1);
        update_self_consistent_potential();
        reset_element_charges();
        return;
    }

    if (common_options().m_initialize_particles_from_doping) {
        place_initial_charges_according_to_doping(common_options().m_initial_particle_weight);
    }
    flatten_particle_positions_for_2d();
}

void self_consistent_device_pbmc_simulation_2d::add_particle_charges_to_elements() {
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

void self_consistent_device_pbmc_simulation_2d::reset_element_charges() {
    m_device.get_p_mesh()->for_each_bulk_element([](mesh::element& element) { element.reset_charge(); });
}

void self_consistent_device_pbmc_simulation_2d::recompute_vertex_space_charge_from_element_charges(
    std::size_t accumulation_steps) {
    if (accumulation_steps == 0) {
        throw std::invalid_argument("accumulation_steps must be positive.");
    }

    const double factor = charge_deposition_factor(accumulation_steps);
    m_device.get_p_mesh()->convert_charge_on_element_into_charge_at_vtx(factor);
}

void self_consistent_device_pbmc_simulation_2d::update_self_consistent_potential(bool publish_mesh_functions) {
    const bool nonlinear_ready = common_options().m_nonlinear_steady_state_poisson &&
                                 m_previous_poisson_solution.size() ==
                                     static_cast<Eigen::Index>(m_device.get_p_mesh()->get_nb_vertices());
    if (nonlinear_ready) {
        update_nonlinear_self_consistent_potential();
    } else {
        m_poisson_solver.update_second_member();
        for (const auto& [contact_name, unused_voltage] : contact_voltages_V()) {
            static_cast<void>(unused_voltage);
            m_poisson_solver.apply_dirichlet_condition_second_member(contact_name,
                                                                     contact_voltage_for_poisson(contact_name));
        }
        m_poisson_solver.solve_system();
    }
    if (common_options().m_enable_poisson_mixing) {
        if (m_previous_poisson_solution.size() == m_poisson_solver.get_solution().size()) {
            if (m_previous_poisson_solution.allFinite()) {
                m_poisson_solver.mix_solution_with(m_previous_poisson_solution,
                                                   common_options().m_poisson_mixing_old_solution_fraction);
            } else {
                fmt::print("WARNING: previous PBMC Poisson solution is non-finite; skipping Poisson mixing for this "
                           "update.\n");
            }
        }
    }
    m_previous_poisson_solution = m_poisson_solver.get_solution();
    if (publish_mesh_functions) {
        constexpr bool add_gradient = true;
        m_poisson_solver.add_solution_to_mesh_functions("PoissonSolution", add_gradient);
    } else {
        m_poisson_solver.update_mesh_electric_field_from_solution();
    }
}

void self_consistent_device_pbmc_simulation_2d::update_nonlinear_self_consistent_potential() {
    auto* mesh = m_device.get_p_mesh();
    const Eigen::Index size = static_cast<Eigen::Index>(mesh->get_nb_vertices());
    Eigen::VectorXd fixed_source = Eigen::VectorXd::Zero(size);
    Eigen::VectorXd lumped_volume = Eigen::VectorXd::Zero(size);
    Eigen::VectorXd electron_density(size);
    Eigen::VectorXd hole_density(size);
    Eigen::VectorXd thermal_voltage(size);
    constexpr double density_conversion = 1.0e-6;

    for (const auto& element : mesh->get_list_bulk_element()) {
        const double nodal_measure = std::abs(element->get_measure()) / 3.0;
        for (const auto* vertex : element->get_vertices()) {
            lumped_volume(vertex->get_index()) += nodal_measure * density_conversion;
        }
    }
    for (Eigen::Index i = 0; i < size; ++i) {
        const auto* vertex = mesh->get_p_vertex(static_cast<std::size_t>(i));
        electron_density(i) = vertex->get_electron_density();
        hole_density(i) = vertex->get_hole_density();
        fixed_source(i) = lumped_volume(i) * vertex->get_doping_concentration();
        thermal_voltage(i) = uepm::constants::k_B * vertex->get_temperature() / uepm::constants::q_e;
    }

    std::vector<uepm::fem::nonlinear_poisson_dirichlet> boundary_conditions;
    std::vector<bool> constrained(static_cast<std::size_t>(size), false);
    for (const auto& [contact_name, unused_voltage] : contact_voltages_V()) {
        static_cast<void>(unused_voltage);
        const auto* region = mesh->get_p_region(contact_name);
        if (region == nullptr) {
            throw std::runtime_error("Poisson contact region '" + contact_name + "' does not exist.");
        }
        const double value = contact_voltage_for_poisson(contact_name);
        for (const auto index : region->get_unique_vertices()) {
            if (!constrained[index]) {
                boundary_conditions.push_back({static_cast<Eigen::Index>(index), value});
                constrained[index] = true;
            }
        }
    }

    const auto result = uepm::fem::solve_nonlinear_poisson_lumped(
        m_unconstrained_poisson_stiffness,
        fixed_source,
        lumped_volume,
        electron_density,
        hole_density,
        m_previous_poisson_solution,
        m_previous_poisson_solution,
        thermal_voltage,
        boundary_conditions,
        common_options().m_nonlinear_poisson_options);
    if (!result.converged) {
        throw std::runtime_error("Nonlinear steady-state Poisson did not converge after " +
                                 std::to_string(result.iterations) + " iterations.");
    }
    fmt::print("Nonlinear Poisson converged in {} iterations, residual={:.6e}, max update={:.6e} V\n",
               result.iterations,
               result.final_residual_norm,
               result.maximum_correction_V);
    m_poisson_solver.set_solution(result.potential_V);
}

void self_consistent_device_pbmc_simulation_2d::run_self_consistent_transport_simulation() {
    const std::size_t total_iterations =
        static_cast<std::size_t>(std::ceil(m_simulation_options.m_t_max / m_simulation_options.m_time_step));

    std::string  history_filename = initialize_simulation_history_file();
    std::fstream stream(history_filename, std::ios::app);

    fmt::print("START 2D SELF-CONSISTENT PBMC SIMULATION\n");
    fmt::print("Total iterations: {}\n", total_iterations);
    fmt::print("Poisson frequency: {}\n", poisson_frequency());

    double accumulator_ramo_current_electron       = 0.0;
    double accumulator_ramo_current_hole           = 0.0;
    double ramo_current_electron                   = 0.0;
    double ramo_current_hole                       = 0.0;
    double ramo_current                            = 0.0;
    double probe_ramo_current_electron             = 0.0;
    double probe_ramo_current_hole                 = 0.0;
    double probe_ramo_current                      = 0.0;

    std::size_t poisson_sample_count = 0;
    std::size_t nonlinear_warmup_steps_remaining =
        common_options().m_nonlinear_steady_state_poisson ? common_options().m_nonlinear_poisson_warmup_steps : 0;
    bool terminated_early = false;

    const auto finalize_poisson_batch = [&](double sample_time_s) {
        if (poisson_sample_count == 0) {
            // A voltage event immediately following a completed batch still
            // needs to update the field before the next transport step.
            update_self_consistent_potential(false);
            return;
        }

        const double sample_count = static_cast<double>(poisson_sample_count);
        const double averaged_ramo_current_electron = accumulator_ramo_current_electron / sample_count;
        const double averaged_ramo_current_hole     = accumulator_ramo_current_hole / sample_count;
        const double averaged_ramo_current = averaged_ramo_current_electron + averaged_ramo_current_hole -
                                             common_options().m_background_ramo_current_A;

        if (m_simulation_options.m_scheduled_particle_injection.m_done) {
            const double circuit_dt_s = m_simulation_options.m_time_step * sample_count;
            advance_quench_circuit(averaged_ramo_current, circuit_dt_s, sample_time_s);
        }

        add_charges_at_contacts(poisson_sample_count);
        add_missing_contact_charge_to_poisson_reservoir(poisson_sample_count);
        recompute_vertex_space_charge_from_element_charges(poisson_sample_count);
        update_self_consistent_potential(false);
        reset_element_charges();

        poisson_sample_count = 0;
        accumulator_ramo_current_electron = 0.0;
        accumulator_ramo_current_hole     = 0.0;
        nonlinear_warmup_steps_remaining = common_options().m_nonlinear_steady_state_poisson
                                               ? common_options().m_nonlinear_poisson_warmup_steps
                                               : 0;

    };

    while (m_state.m_iteration < total_iterations) {
        if (m_state.m_iteration > 10 && m_simulation_options.m_stop_simu_when_no_electron_remaining &&
            get_number_electrons() == 0 && !has_pending_scheduled_particle_injection()) {
            fmt::print("Stop: no electrons remaining in device.\n");
            terminated_early = true;
            break;
        }

        if (has_reached_particle_limit()) {
            fmt::print("Stop: hard particle limit reached.\n");
            terminated_early = true;
            break;
        }

        // Apply scheduled voltages at transport-step resolution. If an event
        // cuts a Poisson window short, finalize that partial window with its
        // actual sample count before propagating under the new field.
        if (apply_scheduled_contact_voltage_events(m_state.m_time_s)) {
            finalize_poisson_batch(m_state.m_time_s);
        }

        if (m_common_options.m_auto_background_ramo_current && has_pending_scheduled_particle_injection() &&
            m_state.m_time_s >= m_simulation_options.m_scheduled_particle_injection.m_time_s &&
            m_common_options.m_background_ramo_current_A == 0.0) {
            constexpr double time_window_s = 1e-12;  // 1 ps
            const double background_current_A = m_simulation_history.extract_final_current(time_window_s);
            fmt::print("Extracted pre-injection background Ramo current: {:.3e} A\n", background_current_A);
            m_common_options.m_background_ramo_current_A = background_current_A;
        }

        const std::size_t impact_events_before_step = m_simulation_history.m_impact_ionization_positions.size();
        transport_particles_one_time_step();

        const auto currents = compute_ramo_currents(true, true);
        const bool is_sampling_step = nonlinear_warmup_steps_remaining == 0;
        if (is_sampling_step) {
            add_particle_charges_to_elements();
            ++poisson_sample_count;
            accumulator_ramo_current_electron += currents.electron;
            accumulator_ramo_current_hole += currents.hole;
        } else {
            --nonlinear_warmup_steps_remaining;
        }

        // History stores the current of this transport step. Batch-averaged
        // currents are computed separately below for the external circuit.
        ramo_current_electron       = currents.electron;
        ramo_current_hole           = currents.hole;
        ramo_current                = ramo_current_electron + ramo_current_hole - common_options().m_background_ramo_current_A;
        probe_ramo_current_electron = currents.probe_electron;
        probe_ramo_current_hole     = currents.probe_hole;
        probe_ramo_current          = probe_ramo_current_electron + probe_ramo_current_hole;

        const bool should_update_poisson = poisson_sample_count == poisson_frequency();

        if (should_update_poisson) {
            const double poisson_sample_time_s = m_state.m_time_s + m_simulation_options.m_time_step;
            finalize_poisson_batch(poisson_sample_time_s);
        }

        m_state.m_time_s += m_simulation_options.m_time_step;
        ++m_state.m_iteration;
        update_successful_quench_detection(m_state.m_time_s, impact_events_before_step);

        const double max_electric_field_V_per_cm = max_particle_electric_field_V_per_cm();
        const auto   mean_energies_eV            = get_mean_kinetic_energies_eV();
        m_simulation_history.add_data_to_history(
            m_state.m_time_s,
            get_number_electrons(),
            get_number_holes(),
            mean_energies_eV[0],
            mean_energies_eV[1],
            mean_energies_eV[2],
            m_simulation_history.m_impact_ionization_positions.size(),
            ramo_current_electron,
            ramo_current_hole,
            ramo_current,
            probe_ramo_current_electron,
            probe_ramo_current_hole,
            probe_ramo_current,
            max_electric_field_V_per_cm,
            ramo_electrode_voltage_for_history(),
            reference_electrode_voltage_for_history(),
            quench_supply_voltage_for_history(),
            quench_device_current_for_history(),
            quench_resistor_current_for_history(),
            quench_voltage_drop_for_history(),
            m_simulation_history.contact_voltage_values_from_map(contact_voltages_V()));

        // Periodic checkpoint. The complete in-memory history overwrites this
        // file after the loop finishes successfully or stops early.
        if (m_state.m_iteration % 10 == 0) {
            m_simulation_history.append_last_iter_to_csv(stream);
        }
        if (m_state.m_iteration == 1 ||
            m_state.m_iteration % static_cast<std::size_t>(m_simulation_options.m_frequency_export_trajectory) == 0) {
            fmt::print("\rExported iteration at time {:<10.3e}ps - {:>9d} / {} ({:.1f}%) -- number of particles: {}",
                       m_state.m_time_s * 1e12,
                       m_state.m_iteration,
                       total_iterations,
                       static_cast<double>(m_state.m_iteration) / static_cast<double>(total_iterations) * 100.0,
                       m_list_particles.size());
            std::fflush(stdout);
            stream.flush();
            if (m_simulation_options.m_export_time_step) {
                constexpr bool add_gradient = true;
                m_poisson_solver.add_solution_to_mesh_functions("PoissonSolution", add_gradient);
                export_current_state();
            }
        }
    }

    if (!terminated_early && poisson_sample_count > 0) {
        finalize_poisson_batch(m_state.m_time_s);
    }
    stream.close();
    m_simulation_history.export_to_csv(history_filename);
    constexpr bool add_gradient = true;
    m_poisson_solver.add_solution_to_mesh_functions("PoissonSolution", add_gradient);
    export_current_state();
    fmt::print("END 2D SELF-CONSISTENT PBMC SIMULATION\n");
}

void self_consistent_device_pbmc_simulation_2d::export_current_state() { export_current_snapshot(); }

}  // namespace uepm::PBMC
