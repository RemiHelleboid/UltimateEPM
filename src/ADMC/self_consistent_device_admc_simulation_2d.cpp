/**
 * @file self_consistent_device_admc_simulation_2d.cpp
 * @brief Self-consistent 2D ADMC device simulation.
 */

#include "self_consistent_device_admc_simulation_2d.hpp"

#include <fmt/core.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <stdexcept>

#include "physical_constants.hpp"
#include "unit_conversion.hpp"

namespace uepm::ADMC {
namespace {

double silicon_band_gap_varshni_eV(double temperature_K) {
    constexpr double eg_0_eV  = 1.17;
    constexpr double alpha_eV = 4.73e-4;
    constexpr double beta_K   = 636.0;
    return eg_0_eV - alpha_eV * temperature_K * temperature_K / (temperature_K + beta_K);
}

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

double contact_equilibrium_voltage_offset_V(const contact_doping_summary& summary, double temperature_K) {
    const double donor_excess    = summary.donor_cm_3 - summary.acceptor_cm_3;
    const double acceptor_excess = summary.acceptor_cm_3 - summary.donor_cm_3;
    const double thermal_voltage = uepm::constants::k_B * temperature_K / uepm::constants::q_e;
    const double intrinsic_concentration_cm_3 = silicon_intrinsic_concentration_admc_cm_3(temperature_K);

    if (donor_excess > 0.0) {
        return thermal_voltage * std::log(donor_excess / intrinsic_concentration_cm_3);
    }
    if (acceptor_excess > 0.0) {
        return -thermal_voltage * std::log(acceptor_excess / intrinsic_concentration_cm_3);
    }
    return 0.0;
}

}  // namespace

double silicon_intrinsic_concentration_admc_cm_3(double temperature_K) {
    if (!std::isfinite(temperature_K) || temperature_K <= 0.0) {
        throw std::invalid_argument("Intrinsic concentration temperature must be positive and finite.");
    }

    constexpr double reference_temperature_K      = 300.0;
    constexpr double reference_concentration_cm_3 = 1.0e10;
    const double     reference_band_gap_eV        = silicon_band_gap_varshni_eV(reference_temperature_K);
    const double     band_gap_eV                  = silicon_band_gap_varshni_eV(temperature_K);
    const double     effective_density_ratio      = std::pow(temperature_K / reference_temperature_K, 1.5);
    const double     band_gap_boltzmann_correction =
        std::exp(reference_band_gap_eV / (2.0 * uepm::constants::k_b_eV * reference_temperature_K) -
                 band_gap_eV / (2.0 * uepm::constants::k_b_eV * temperature_K));

    return reference_concentration_cm_3 * effective_density_ratio * band_gap_boltzmann_correction;
}

void options_self_consistent_device_ADMC_common::validate() const {
    if (m_poisson_frequency == 0) {
        throw std::invalid_argument("ADMC Poisson frequency must be positive.");
    }
    if (m_contact_voltages_V.empty()) {
        throw std::invalid_argument("ADMC contacts.voltages_V must define at least one contact.");
    }
    for (const auto& [contact_name, voltage_V] : m_contact_voltages_V) {
        if (contact_name.empty()) {
            throw std::invalid_argument("ADMC contact names must not be empty.");
        }
        if (!std::isfinite(voltage_V)) {
            throw std::invalid_argument("ADMC contact voltages must be finite.");
        }
    }
    if (m_ramo_electrode.empty() || !m_contact_voltages_V.contains(m_ramo_electrode)) {
        throw std::invalid_argument("ADMC Ramo electrode must name a configured contact.");
    }
    if (!std::isfinite(m_built_in_contact_voltage_scale)) {
        throw std::invalid_argument("ADMC built-in contact voltage scale must be finite.");
    }
    if (!std::isfinite(m_poisson_mixing_old_solution_fraction) ||
        m_poisson_mixing_old_solution_fraction < 0.0 || m_poisson_mixing_old_solution_fraction > 1.0) {
        throw std::invalid_argument("ADMC Poisson mixing old solution fraction must be in [0, 1].");
    }
    if (!std::isfinite(m_initial_particle_weight) || m_initial_particle_weight <= 0.0) {
        throw std::invalid_argument("ADMC initial particle weight must be positive and finite.");
    }
    if (!std::isfinite(m_contact_injection_particle_weight) || m_contact_injection_particle_weight <= 0.0) {
        throw std::invalid_argument("ADMC contact injection particle weight must be positive and finite.");
    }
}

void options_self_consistent_device_ADMC_2d::validate() const {
    m_common.validate();
    if (!std::isfinite(m_effective_depth_um) || m_effective_depth_um <= 0.0) {
        throw std::invalid_argument("ADMC effective depth must be positive and finite.");
    }
}

self_consistent_device_admc_simulation_2d::self_consistent_device_admc_simulation_2d(
    const device::device& simulation_device,
    const options_device_ADMC& simulation_options,
    const options_self_consistent_device_ADMC_2d& self_consistent_options,
    const physics::material_database& material_database,
    const mesh::vector3& starting_position_um,
    std::size_t number_electrons_start,
    std::size_t number_holes_start,
    std::uint64_t random_seed)
    : device_admc_simulation(simulation_device,
                             simulation_options,
                             starting_position_um,
                             number_electrons_start,
                             number_holes_start,
                             random_seed),
      m_self_consistent_options(self_consistent_options),
      m_poisson_solver(m_device.get_p_mesh(), m_device.get_p_mesh()->get_nb_vertices(), material_database),
      m_contact_rng(static_cast<unsigned long>(random_seed + 1)) {
    validate_self_consistent_options();
    initialize_contact_elements();
    initialize_poisson_solver();
    initialize_particles_for_self_consistent_run();
}

void self_consistent_device_admc_simulation_2d::validate_self_consistent_options() const {
    if (m_dimension != 2) {
        throw std::invalid_argument("self_consistent_device_admc_simulation_2d requires a 2D device.");
    }
    m_self_consistent_options.validate();
}

double self_consistent_device_admc_simulation_2d::scale_integrated_2d_doping_to_carriers(
    double integrated_doping) const {
    return integrated_doping * m_self_consistent_options.m_effective_depth_um * uepm::units::micron_to_cm;
}

double self_consistent_device_admc_simulation_2d::contact_voltage_for_poisson(const std::string& contact_name) const {
    const auto contact_it = m_self_consistent_options.m_common.m_contact_voltages_V.find(contact_name);
    if (contact_it == m_self_consistent_options.m_common.m_contact_voltages_V.end()) {
        throw std::invalid_argument("Unknown ADMC configured contact '" + contact_name + "'.");
    }
    const auto offset_it = m_built_in_contact_voltage_offsets_V.find(contact_name);
    return contact_it->second + (offset_it == m_built_in_contact_voltage_offsets_V.end() ? 0.0 : offset_it->second);
}

void self_consistent_device_admc_simulation_2d::update_built_in_contact_voltage_offset(
    const std::string& contact_name,
    const std::vector<std::shared_ptr<mesh::element>>& contact_elements) {
    m_built_in_contact_voltage_offsets_V[contact_name] = 0.0;
    if (!m_self_consistent_options.m_common.m_enable_built_in_potential) {
        return;
    }
    const auto   doping   = summarize_contact_doping(contact_elements);
    const double offset_V = m_self_consistent_options.m_common.m_built_in_contact_voltage_scale *
                            contact_equilibrium_voltage_offset_V(doping, m_options.m_lattice_temperature_K);
    m_built_in_contact_voltage_offsets_V[contact_name] = offset_V;
    fmt::print("ADMC built-in contact '{}': Nd={:.6e} cm^-3, Na={:.6e} cm^-3, offset={:.6e} V\n",
               contact_name,
               doping.donor_cm_3,
               doping.acceptor_cm_3,
               offset_V);
}

void self_consistent_device_admc_simulation_2d::initialize_contact_elements() {
    const auto list_bulk_elements = m_device.get_p_mesh()->get_list_bulk_element();

    m_list_element_contact.clear();
    m_list_element_contact_ptr.clear();
    m_list_element_contact_equilibrium_charge.clear();

    for (const auto& device_contact : m_device.get_list_contacts()) {
        const std::string contact_name = device_contact.get_contact_name();
        if (!m_self_consistent_options.m_common.m_contact_voltages_V.contains(contact_name)) {
            throw std::runtime_error("ADMC collecting contact '" + contact_name +
                                     "' has no configured Poisson voltage.");
        }
        const auto element_indices =
            m_device.get_p_mesh()->get_idx_bulk_elements_adjacent_to_contact_region(contact_name);
        std::vector<std::shared_ptr<mesh::element>> contact_elements;
        contact_elements.reserve(element_indices.size());

        for (const auto element_index : element_indices) {
            if (element_index >= list_bulk_elements.size()) {
                throw std::runtime_error("ADMC contact-adjacent bulk element index is out of range.");
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

void self_consistent_device_admc_simulation_2d::compute_unitary_potential() {
    m_poisson_solver.compute_second_member(0.0);
    for (const auto& [contact_name, unused_voltage] : m_self_consistent_options.m_common.m_contact_voltages_V) {
        static_cast<void>(unused_voltage);
        m_poisson_solver.apply_dirichlet_condition(
            contact_name, contact_name == m_self_consistent_options.m_common.m_ramo_electrode ? 1.0 : 0.0);
    }
    m_poisson_solver.decompose_matrix();
    m_poisson_solver.solve_system();
    m_poisson_solver.add_solution_to_mesh_functions("RamoUnitaryPotential", true);

    const double unitary_potential_gradient_max =
        m_device.get_p_mesh()->get_argmax_max_of_function("RamoUnitaryPotential_gradient_norm").second;
    const double unitary_potential_gradient_min =
        m_device.get_p_mesh()->get_argmin_min_of_function("RamoUnitaryPotential_gradient_norm").second;
    if (std::abs(unitary_potential_gradient_max - unitary_potential_gradient_min) < 1e-6) {
        fmt::print("INFO : ADMC unitary potential gradient is constant across the device "
                   "(max = {:.6e}, min = {:.6e}). The constant Ramo field will be used.\n",
                   unitary_potential_gradient_max,
                   unitary_potential_gradient_min);
        m_state.m_use_constant_RamoUnitaryElectricField = true;
        m_state.m_RamoUnitaryElectricField_Vm_per_cm =
            m_device.get_p_mesh()->interpolate_vector_at_location("RamoUnitaryPotential_gradient", {1e-3, 1e-3, 0.0});
    } else {
        m_state.m_use_constant_RamoUnitaryElectricField = false;
        fmt::print("ADMC unitary potential gradient computed. Max value: {:.6e}, Min value: {:.6e}\n",
                   unitary_potential_gradient_max,
                   unitary_potential_gradient_min);
    }
}

void self_consistent_device_admc_simulation_2d::initialize_poisson_solver() {
    m_poisson_solver.compute_stiffness_matrix();
    compute_unitary_potential();

    m_poisson_solver.compute_stiffness_matrix();
    m_poisson_solver.compute_second_member(0.0);
    for (const auto& [contact_name, unused_voltage] : m_self_consistent_options.m_common.m_contact_voltages_V) {
        static_cast<void>(unused_voltage);
        m_poisson_solver.apply_dirichlet_condition(contact_name, contact_voltage_for_poisson(contact_name));
    }
    m_poisson_solver.decompose_matrix();
}

void self_consistent_device_admc_simulation_2d::update_self_consistent_potential(bool publish_mesh_functions) {
    m_poisson_solver.update_second_member();
    for (const auto& [contact_name, unused_voltage] : m_self_consistent_options.m_common.m_contact_voltages_V) {
        static_cast<void>(unused_voltage);
        m_poisson_solver.apply_dirichlet_condition_second_member(contact_name,
                                                                 contact_voltage_for_poisson(contact_name));
    }
    m_poisson_solver.solve_system();
    if (m_self_consistent_options.m_common.m_enable_poisson_mixing) {
        if (m_previous_poisson_solution.size() == m_poisson_solver.get_solution().size()) {
            if (m_previous_poisson_solution.allFinite()) {
                m_poisson_solver.mix_solution_with(
                    m_previous_poisson_solution,
                    m_self_consistent_options.m_common.m_poisson_mixing_old_solution_fraction);
            } else {
                fmt::print("WARNING: previous ADMC Poisson solution is non-finite; skipping Poisson mixing for this "
                           "update.\n");
            }
        }
        m_previous_poisson_solution = m_poisson_solver.get_solution();
    }
    if (publish_mesh_functions) {
        m_poisson_solver.add_solution_to_mesh_functions("PoissonSolution", true);
    } else {
        m_poisson_solver.update_mesh_electric_field_from_solution();
    }
}

void self_consistent_device_admc_simulation_2d::initialize_particles_for_self_consistent_run() {
    if (m_self_consistent_options.m_common.m_initialize_particles_from_doping) {
        place_initial_charges_according_to_doping(m_self_consistent_options.m_common.m_initial_particle_weight);
    }
    reset_element_charges();
    add_particle_charges_to_elements();
    recompute_vertex_space_charge_from_element_charges(1);
    update_self_consistent_potential();
    reset_element_charges();
}

void self_consistent_device_admc_simulation_2d::place_initial_charges_according_to_doping(double particle_weight) {
    if (particle_weight <= 0.0) {
        throw std::invalid_argument("ADMC particle weight must be positive.");
    }

    auto* mesh = m_device.get_p_mesh();
    const auto integrate_carriers_over_2d_mesh = [&](const std::string& field_name) {
        double total_charge = 0.0;
        for (const auto& element : mesh->get_list_bulk_element()) {
            total_charge += scale_integrated_2d_doping_to_carriers(element->integrate_scalar(field_name));
        }
        return total_charge;
    };

    const std::string donor_field_name    = "DonorConcentration";
    const std::string acceptor_field_name = "AcceptorConcentration";
    const double      total_donor_charge = integrate_carriers_over_2d_mesh(donor_field_name);
    const double      total_acceptor_charge = integrate_carriers_over_2d_mesh(acceptor_field_name);
    const std::size_t number_electrons = static_cast<std::size_t>(std::floor(total_donor_charge / particle_weight));
    const std::size_t number_holes     = static_cast<std::size_t>(std::floor(total_acceptor_charge / particle_weight));

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
    const double max_donor_concentration = mesh->get_argmax_max_of_function(donor_field_name).second;
    const double max_acceptor_concentration = mesh->get_argmax_max_of_function(acceptor_field_name).second;
    const mesh::bbox active_region_bbox = mesh->get_p_region("Silicon_1")->compute_bounding_box();

    if (number_electrons > 0 && max_donor_concentration <= 0.0) {
        throw std::runtime_error("Donor concentration maximum is non-positive.");
    }

    if (number_holes > 0 && max_acceptor_concentration <= 0.0) {
        throw std::runtime_error("Acceptor concentration maximum is non-positive.");
    }

    while (electron_positions.size() < number_electrons) {
        const mesh::vector3 position = active_region_bbox.draw_uniform_random_point_inside_box(m_contact_rng);
        const double donor_density = mesh->interpolate_scalar_at_location(donor_field_name, position);
        const double acceptor_density = mesh->interpolate_scalar_at_location(acceptor_field_name, position);
        if (acceptor_density > donor_density) {
            continue;
        }

        const double probability = donor_density / max_donor_concentration;
        if (probability > min_probability && uniform01(m_contact_rng) < probability) {
            electron_positions.push_back(position);
        }
    }
    while (hole_positions.size() < number_holes) {
        const mesh::vector3 position = active_region_bbox.draw_uniform_random_point_inside_box(m_contact_rng);
        const double acceptor_density = mesh->interpolate_scalar_at_location(acceptor_field_name, position);
        const double donor_density = mesh->interpolate_scalar_at_location(donor_field_name, position);
        if (donor_density > acceptor_density) {
            continue;
        }

        const double probability = acceptor_density / max_acceptor_concentration;
        if (probability > min_probability && uniform01(m_contact_rng) < probability) {
            hole_positions.push_back(position);
        }
    }

    m_particles.clear();
    m_particles.reserve(electron_positions.size() + hole_positions.size());
    for (const auto& position : electron_positions) {
        add_particle_at_position(position, carrier_type::electron, particle_weight);
    }
    for (const auto& position : hole_positions) {
        add_particle_at_position(position, carrier_type::hole, particle_weight);
    }

    fmt::print("Initial particles placed according to doping.\n");
}

void self_consistent_device_admc_simulation_2d::add_charges_at_contacts(std::size_t poisson_frequency) {
    const double particle_weight = m_self_consistent_options.m_common.m_contact_injection_particle_weight;
    if (m_list_element_contact_ptr.empty()) {
        return;
    }

    std::vector<double> electron_charge_to_add(m_list_element_contact_ptr.size(), 0.0);
    std::vector<double> hole_charge_to_add(m_list_element_contact_ptr.size(), 0.0);
    double total_electron_charge_to_add = 0.0;
    double total_hole_charge_to_add     = 0.0;

    for (std::size_t i = 0; i < m_list_element_contact_ptr.size(); ++i) {
        auto& element = m_list_element_contact_ptr[i];
        const double element_charge = element->get_n_charge() - element->get_p_charge();
        const double averaged_element_charge = element_charge / static_cast<double>(poisson_frequency);
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
    for (std::size_t count = 0; count < number_electrons_to_place && !has_reached_particle_limit(); ++count) {
        const std::size_t i = contact_index_distribution(m_contact_rng);
        if (electron_charge_to_add[i] <= 0.0) {
            continue;
        }
        add_particle_at_position(m_list_element_contact_ptr[i]->draw_uniform_random_point_inside_element(m_contact_rng),
                                 carrier_type::electron,
                                 particle_weight);
        electron_charge_to_add[i] -= particle_weight;
    }
    for (std::size_t count = 0; count < number_holes_to_place && !has_reached_particle_limit(); ++count) {
        const std::size_t i = contact_index_distribution(m_contact_rng);
        if (hole_charge_to_add[i] <= 0.0) {
            continue;
        }
        add_particle_at_position(m_list_element_contact_ptr[i]->draw_uniform_random_point_inside_element(m_contact_rng),
                                 carrier_type::hole,
                                 particle_weight);
        hole_charge_to_add[i] -= particle_weight;
    }
}

void self_consistent_device_admc_simulation_2d::add_missing_contact_charge_to_poisson_reservoir(
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

void self_consistent_device_admc_simulation_2d::reset_element_charges() {
    m_device.get_p_mesh()->for_each_bulk_element([](mesh::element& element) { element.reset_charge(); });
}

void self_consistent_device_admc_simulation_2d::add_particle_charges_to_elements() {
    for (const auto& device_particle : m_particles) {
        auto* element = device_particle.containing_element;
        if (element == nullptr) {
            continue;
        }
        if (device_particle.particle.type() == carrier_type::electron) {
            element->add_n_charge(device_particle.weight);
        } else {
            element->add_p_charge(device_particle.weight);
        }
    }
}

void self_consistent_device_admc_simulation_2d::recompute_vertex_space_charge_from_element_charges(
    std::size_t accumulation_steps) {
    if (accumulation_steps == 0) {
        throw std::invalid_argument("ADMC accumulation steps must be positive.");
    }
    const double factor =
        1.0 / (static_cast<double>(accumulation_steps) * m_self_consistent_options.m_effective_depth_um);
    m_device.get_p_mesh()->convert_charge_on_element_into_charge_at_vtx(factor);
}

void self_consistent_device_admc_simulation_2d::run_self_consistent_transport_simulation() {
    const std::size_t total_iterations =
        static_cast<std::size_t>(std::ceil(m_options.m_final_time_s / m_options.m_time_step_s));
    const std::size_t poisson_frequency = m_self_consistent_options.m_common.m_poisson_frequency;
    double accumulator_ramo_current_electron = 0.0;
    double accumulator_ramo_current_hole     = 0.0;
    double ramo_current_electron             = 0.0;
    double ramo_current_hole                 = 0.0;
    double ramo_current                      = 0.0;
    const std::string history_filename       = initialize_simulation_history_file();
    std::fstream      history_stream(history_filename, std::ios::app);
    history_stream << std::setprecision(std::numeric_limits<double>::max_digits10);
    if (!history_stream.is_open()) {
        throw std::runtime_error("Could not open ADMC simulation history CSV file '" + history_filename + "'.");
    }
    std::size_t last_history_append_iteration = 0;

    fmt::print("START 2D SELF-CONSISTENT ADMC SIMULATION\n");
    fmt::print("Total iterations: {}\n", total_iterations);
    fmt::print("Poisson frequency: {}\n", poisson_frequency);

    reset_element_charges();
    while (m_state.m_time_s < m_options.m_final_time_s) {
        if (m_particles.empty() && !has_pending_scheduled_particle_injection()) {
            fmt::print("Stop: no particles remaining in device.\n");
            break;
        }
        if (m_options.m_stop_when_no_electrons && get_number_electrons() == 0 &&
            !has_pending_scheduled_particle_injection()) {
            fmt::print("Stop: no electrons remaining in device.\n");
            break;
        }
        if (has_reached_particle_limit()) {
            fmt::print("Stop: hard particle limit reached.\n");
            break;
        }

        advance_particles_one_time_step();
        add_particle_charges_to_elements();

        const auto [electron_current, hole_current] = last_ramo_current();
        accumulator_ramo_current_electron += electron_current;
        accumulator_ramo_current_hole += hole_current;

        if (m_state.m_iteration == 1 || m_state.m_iteration % 10 == 0) {
            m_history.append_last_iter_to_csv(history_stream);
            last_history_append_iteration = m_state.m_iteration;
        }

        const bool should_update_poisson = (m_state.m_iteration % poisson_frequency == 0) && m_state.m_iteration != 0;
        if (should_update_poisson) {
            ramo_current_electron = accumulator_ramo_current_electron / static_cast<double>(poisson_frequency);
            ramo_current_hole     = accumulator_ramo_current_hole / static_cast<double>(poisson_frequency);
            ramo_current          = ramo_current_electron + ramo_current_hole;
            accumulator_ramo_current_electron = 0.0;
            accumulator_ramo_current_hole     = 0.0;

            add_charges_at_contacts(poisson_frequency);
            add_particle_charges_to_elements();
            add_missing_contact_charge_to_poisson_reservoir(poisson_frequency + 1);
            recompute_vertex_space_charge_from_element_charges(poisson_frequency + 1);
            update_self_consistent_potential(false);
            reset_element_charges();

        }

        if (m_state.m_iteration == 1 ||
            m_state.m_iteration % static_cast<std::size_t>(m_options.m_frequency_export) == 0) {
            fmt::print("\rExported iteration at time {:<10.3e}ps - {:>9d} / {} ({:.1f}%) -- "
                       "number of particles: {} -- Ramo current: {:.6e} A",
                       m_state.m_time_s * 1e12,
                       m_state.m_iteration,
                       total_iterations,
                       total_iterations == 0
                           ? 100.0
                           : static_cast<double>(m_state.m_iteration) / static_cast<double>(total_iterations) * 100.0,
                       m_particles.size(),
                       ramo_current);
            std::fflush(stdout);
            history_stream.flush();
            if (m_options.m_export_time_step) {
                m_poisson_solver.add_solution_to_mesh_functions("PoissonSolution", true);
                export_current_snapshot();
            }
        }
    }
    if (m_state.m_iteration > 0 && last_history_append_iteration != m_state.m_iteration) {
        m_history.append_last_iter_to_csv(history_stream);
    }
    history_stream.close();
    m_poisson_solver.add_solution_to_mesh_functions("PoissonSolution", true);
    export_current_snapshot();
    fmt::print("\n");
    fmt::print("END 2D SELF-CONSISTENT ADMC SIMULATION\n");
}

}  // namespace uepm::ADMC
