/**
 * @file device_amc_simulation.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-04
 *
 *
 */

#include "device_amc_simulation.hpp"

#include <fmt/chrono.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <fmt/xchar.h>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "physical_constants.hpp"

namespace uepm::amc {

amc_transport_config device_amc_simulation::make_transport_config(const options_device_amc &options, particle_type carrier_type) {
    amc_transport_config cfg;
    cfg.m_carrier_type                  = carrier_type;
    cfg.m_lattice_temperature           = options.m_lattice_temperature;
    cfg.m_max_energy_eV                 = options.m_max_energy_eV;
    cfg.m_self_scattering_safety_factor = options.m_self_scattering_safety_factor;
    cfg.m_gamma_max_energy_samples      = options.m_gamma_max_energy_samples;
    return cfg;
}

amc_transport_kernel &device_amc_simulation::transport_for(particle_type type) {
    if (type == particle_type::electron) {
        return m_electron_transport;
    }

    return m_hole_transport;
}

const amc_transport_kernel &device_amc_simulation::transport_for(particle_type type) const {
    if (type == particle_type::electron) {
        return m_electron_transport;
    }

    return m_hole_transport;
}

void device_amc_simulation::initialize_particle_transport_state(particle_amc &particle) {
    auto &transport = transport_for(particle.type());

    if (transport.valleys().empty()) {
        throw std::runtime_error("transport kernel has no valleys/bands");
    }

    particle.state().valley_index = particle.index() % transport.valleys().size();

    transport.initialize_particle_state(particle);

    if (m_simulation_options.m_keep_particles_history) {
        particle.record_state();
    }
}

device_amc_simulation::device_amc_simulation(const device::device     &simulation_device,
                                             const options_device_amc &simulation_option,
                                             const std::string        &simulation_name,
                                             int                       seed_random_generator)
    : m_device(simulation_device),
      m_electron_transport(make_transport_config(simulation_option, particle_type::electron), seed_random_generator),
      m_hole_transport(make_transport_config(simulation_option, particle_type::hole), seed_random_generator + 1),
      m_dimension(m_device.get_dimension()),
      m_simulation_options(simulation_option),
      m_simulation_name(simulation_name),
      m_iteration{0}{
    m_simulation_history.m_initial_seed_rng = seed_random_generator;
    m_electron_transport.initialize();
    m_hole_transport.initialize();
}

device_amc_simulation::device_amc_simulation(const device::device     &device_simulation,
                                             const options_device_amc &simulation_option,
                                             const std::string        &simulation_name,
                                             const mesh::vector3      &starting_position,
                                             std::size_t               number_electrons_start,
                                             std::size_t               number_holes_start,
                                             int                       seed_random_generator)
    : m_device(device_simulation),
      m_electron_transport(make_transport_config(simulation_option, particle_type::electron), seed_random_generator),
      m_hole_transport(make_transport_config(simulation_option, particle_type::hole), seed_random_generator + 1),
      m_dimension(m_device.get_dimension()),
      m_simulation_options(simulation_option),
      m_simulation_name(simulation_name),
      m_iteration{0} {
    m_electron_transport.initialize();
    m_hole_transport.initialize();
    mesh::element *first_element{nullptr};
    if (m_dimension == 2) {
        first_element = m_device.find_element_at_location(starting_position.to_2d());
    } else {
        first_element = m_device.find_element_at_location(starting_position);
    }
    if (first_element == nullptr) {
        std::cout << "Error : particle can't find its first element. No particle created.    " << starting_position << std::endl;
        return;
    }
    // Creation of electrons and then holes
    m_list_particles.reserve(number_electrons_start + number_holes_start);

    for (std::size_t i = 0; i < number_electrons_start; ++i) {
        const std::size_t particle_index = m_list_particles.size();
        m_list_particles.push_back(std::make_unique<particle_amc>(particle_index, particle_type::electron));
    }

    for (std::size_t i = 0; i < number_holes_start; ++i) {
        const std::size_t particle_index = m_list_particles.size();
        m_list_particles.push_back(std::make_unique<particle_amc>(particle_index, particle_type::hole));
    }

    //  Setup the initial position (unique), random number for RPLA and set containing elements for each particles.
    for (auto &p_particle : m_list_particles) {
        p_particle->set_position(starting_position);
        p_particle->set_weight(1.0);
        p_particle->set_containing_element(first_element);

        initialize_particle_transport_state(*p_particle);
    }

    m_simulation_history.m_initial_seed_rng = seed_random_generator;
}

void device_amc_simulation::add_particle_at_position(const mesh::vector3 &location, particle_type type_of_particle, double weight) {
    mesh::element *first_element{nullptr};
    if (m_dimension == 2) {
        first_element = m_device.find_element_at_location(location.to_2d());
    } else {
        first_element = m_device.find_element_at_location(location);
    }
    if (first_element == nullptr) {
        std::cout << "Error : particle can't find its first element. No particle created.    " << location << std::endl;
        return;
    }
    const std::size_t idx_particle = m_list_particles.size();
    particle_state    initial_state{};
    initial_state.position = location;
    if (type_of_particle == particle_type::electron) {
        m_list_particles.push_back(std::make_unique<particle_amc>(idx_particle, particle_type::electron, initial_state, weight));
    } else {
        m_list_particles.push_back(std::make_unique<particle_amc>(idx_particle, particle_type::hole, initial_state, weight));
    }
    auto &particle = *m_list_particles.back();

    particle.set_containing_element(first_element);
    particle.set_weight(weight);

    initialize_particle_transport_state(particle);
}

void device_amc_simulation::add_particles_at_positions(const std::vector<mesh::vector3> &positions,
                                                       particle_type                     type_of_particle,
                                                       double                            weight) {
    // Reserve memory for all particles upfront
    m_list_particles.reserve(m_list_particles.size() + positions.size());
    for (const auto &location : positions) {
        mesh::element *first_element{nullptr};
        if (m_dimension == 2) {
            first_element = m_device.find_element_at_location(location.to_2d());
        } else {
            first_element = m_device.find_element_at_location(location);
        }

        if (first_element == nullptr) {
            std::cout << "Error : particle can't find its first element. No particle created.    " << location << std::endl;
            continue;
        }

        const std::size_t idx_particle       = m_list_particles.size();
        particle_state    initial_state{};
        initial_state.position = location;
        if (type_of_particle == particle_type::electron) {
            m_list_particles.push_back(std::make_unique<particle_amc>(idx_particle, particle_type::electron, initial_state, weight));
        } else {
            m_list_particles.push_back(std::make_unique<particle_amc>(idx_particle, particle_type::hole, initial_state, weight));
        }
        auto &particle = *m_list_particles.back();

        particle.set_containing_element(first_element);
        particle.set_weight(weight);

        initialize_particle_transport_state(particle);
    }
}

std::size_t device_amc_simulation::get_number_electrons() const {
    std::size_t nb_electron =
        std::accumulate(m_list_particles.begin(), m_list_particles.end(), 0, [](const std::size_t nb_part, const auto &p_part) {
            return nb_part + static_cast<std::size_t>(p_part->type() == particle_type::electron);
        });
    return nb_electron;
}

std::size_t device_amc_simulation::get_number_holes() const {
    std::size_t nb_hole =
        std::accumulate(m_list_particles.begin(), m_list_particles.end(), 0, [](const std::size_t nb_part, const auto &p_part_2) {
            return nb_part + static_cast<std::size_t>(p_part_2->type() == particle_type::hole);
        });
    return nb_hole;
}

double device_amc_simulation::compute_ramo_current() const {
    double current = 0.0;
    for (const auto &p_particle : m_list_particles) {
        // current += p_particle->get_signed_charge() * p_particle->state().velocity.dot(p_particle->state().);
        current += 0.0;  // TODO : compute the weighting field and use it to compute the Ramo current
    }

    current /= static_cast<double>(m_list_particles.size());
    return current;
}

void device_amc_simulation::transport_particles_one_time_step() {
    const double dt = m_simulation_options.m_time_step;

    for (auto &p_particle : m_list_particles) {
        auto &particle = *p_particle;

        particle.set_data_from_device(m_dimension);

        auto &transport = transport_for(particle.type());

        transport.drift_particle(particle, particle.state().electric_field, dt);
        transport.scatter_particle(particle, dt);

        if (m_simulation_options.m_keep_particles_history) {
            particle.record_state();
        }
    }

    update_element_and_check_boundary();
    remove_collected_particles();
}

void device_amc_simulation::advance_particles_one_time_step() {
    transport_particles_one_time_step();
    m_time += m_simulation_options.m_time_step;
    ++m_iteration;
}

void device_amc_simulation::set_particles_transport_data_from_device() {
    for (auto &p_particle : m_list_particles) {
        p_particle->set_data_from_device(m_dimension);
    }
}

void device_amc_simulation::update_element_and_check_boundary() {
    const bool is_2d = m_dimension == 2;
    for (auto &p_particle : m_list_particles) {
        const mesh::element *old_element      = p_particle->get_containing_element();
        mesh::vector3        current_position = p_particle->state().position;
        if (is_2d) {
            current_position.to_2d_inplace();
        }
        if (old_element->is_location_inside_element(current_position)) {
            continue;
        }
        if (m_device.check_enters_contact(current_position)) {
            p_particle->set_crossed_contact(true);
            continue;
        }
        auto *new_element = m_device.find_element_at_location(current_position);
        if (new_element == nullptr) {
            // p_particle->reset_to_previous_position();
            // TODO IMPLEMEMENT
            continue;
        } else if (m_device.get_material_name_at_element(new_element) != "Silicon") {
            // p_particle->reset_to_previous_position();
            // TODO IMPLEMEMENT
            continue;
        } else {
            p_particle->set_containing_element(new_element);
        }
    }
}

void device_amc_simulation::remove_collected_particles() {
    m_anode_current      = 0.0;
    m_cathode_current    = 0.0;
    bool remove_particle = false;
    int  nb_part_erased  = 0;
    for (const auto &p_particle : m_list_particles) {
        if (p_particle->state().m_crossed_contact) {
            remove_particle = true;
            if (m_simulation_options.m_export_time_step) {
            }
            nb_part_erased++;
        }
    }
    if (remove_particle) {
        std::erase_if(m_list_particles, [](auto &&p_part) { return p_part->state().m_crossed_contact; });
    }
    std::vector<double> currents = m_device.get_electrode_currents();
    m_anode_current              = uepm::constants::q_e * currents[0] / m_simulation_options.m_time_step;
    m_cathode_current            = uepm::constants::q_e * currents[1] / m_simulation_options.m_time_step;
}

void device_amc_simulation::run() {
    while (m_time < m_simulation_options.m_t_max) {
        if (m_list_particles.empty()) {
            break;
        }

        if (m_simulation_options.m_stop_simu_when_no_electron_remaining && get_number_electrons() == 0) {
            break;
        }

        if (has_reached_avalanche()) {
            break;
        }

        advance_particles_one_time_step();

        const auto nb_electrons         = get_number_electrons();
        const auto nb_holes             = get_number_holes();
        const auto nb_impact_ionization = m_simulation_history.m_impact_ionization_positions.size();

        m_simulation_history
            .add_data_to_history(m_time, nb_electrons, nb_holes, nb_impact_ionization, m_anode_current, m_cathode_current, 0.0);

        if (m_simulation_options.m_export_time_step &&
            m_iteration % static_cast<std::size_t>(m_simulation_options.m_frequency_export_trajectory) == 0) {
            export_current_time_step_as_csv(m_prefix_export_filename);
        }
    }
}

std::vector<mesh::vector3> device_amc_simulation::get_all_particles_position() const {
    std::vector<mesh::vector3> all_positions(m_list_particles.size());
    std::transform(m_list_particles.begin(), m_list_particles.end(), all_positions.begin(), [](auto &&p_particle) {
        return p_particle->state().position;
    });
    return all_positions;
}

// std::vector<mesh::vector3> device_amc_simulation::get_all_global_velocities() const {
//     std::vector<mesh::vector3> all_velocities(m_list_particles.size());
//     std::transform(m_list_particles.begin(), m_list_particles.end(), all_velocities.begin(), [](auto &&p_particle) {
//         return p_particle->compute_global_velocities();
//     });
//     return all_velocities;
// }

// std::vector<std::size_t> device_amc_simulation::get_all_number_impact_ionization() const {
//     std::vector<std::size_t> all_number_impact_ionization(m_list_particles.size());
//     std::transform(m_list_particles.begin(), m_list_particles.end(), all_number_impact_ionization.begin(), [](auto &&p_particle) {
//         return p_particle->get_total_number_impact_ionization();
//     });
//     return all_number_impact_ionization;
// }

// std::vector<mesh::vector3> device_amc_simulation::get_all_positions_impact_ionization() const {
//     return m_simulation_history.m_impact_ionization_positions;
// }

// std::vector<mesh::vector3> device_amc_simulation::get_all_first_impact_ionization_position() const {
//     std::vector<mesh::vector3> all_impact_ionization_position;
//     for (const auto &p_particle : m_list_particles) {
//         auto impact_ionization_position = p_particle->get_first_impact_ionization_position();
//         if (impact_ionization_position.has_value()) {
//             mesh::vector3 first_i_position = impact_ionization_position.value();
//             all_impact_ionization_position.push_back(first_i_position);
//         }
//     }
//     return all_impact_ionization_position;
// }

std::pair<double, double> device_amc_simulation::compute_depletion_region() const {
    double x_min = std::numeric_limits<double>::max();
    double x_max = std::numeric_limits<double>::max();
    // Compute the distance of the particles to the center of the device. xmin is the maximal distance to the center towards the anode
    // and xmax is the maximal distance to the center towards the cathode.
    double center_x = m_device.get_p_mesh()->get_bounding_box().get_x_min() +
                      0.5 * (m_device.get_p_mesh()->get_bounding_box().get_x_max() - m_device.get_p_mesh()->get_bounding_box().get_x_min());
    for (const auto &p_particle : m_list_particles) {
        double distance_to_center = p_particle->state().position.x() - center_x;
        if (distance_to_center < 0) {
            x_min = std::min(std::abs(x_min), std::abs(distance_to_center));
        } else {
            x_max = std::min(x_max, distance_to_center);
        }
    }
    x_min = center_x - std::abs(x_min);
    x_max = center_x + std::abs(x_max);
    return {x_min, x_max};
}

void device_amc_simulation::export_current_time_step_as_csv(const std::string &prefix_filename) const {
    std::vector<double> x_positions;
    std::vector<double> y_positions;
    std::vector<double> z_positions;
    std::vector<double> particles_times;
    std::vector<double> particle_electric_field;
    std::vector<double> particles_charge;
    for (auto &&p_particle : m_list_particles) {
        x_positions.push_back(p_particle->state().position.x());
        y_positions.push_back(p_particle->state().position.y());
        z_positions.push_back(p_particle->state().position.z());
        particles_times.push_back(p_particle->state().time);
        particle_electric_field.push_back(p_particle->state().electric_field.norm());
        particles_charge.push_back(p_particle->get_signed_charge());
    }
    std::vector<std::string> list_column_names  = {"time", "X", "Y", "Z", "electric field", "type"};
    std::string              iteration_filename = fmt::format("{}.{:09d}.csv", prefix_filename, m_iteration);
    utils::export_multiple_vector_to_csv(
        iteration_filename,
        list_column_names,
        {particles_times, x_positions, y_positions, z_positions, particle_electric_field, particles_charge});
}

void device_amc_simulation::export_all_trajectories_as_csv(const std::string &prefix_filename) const {
    for (auto &&p_particle : m_list_particles) {
        // std::cout << "\rExporting trajectory of particle " << p_particle->get_index() << std::flush;
        std::string filename = fmt::format("{}particle_{:06}_trajectory.csv", prefix_filename, p_particle->index());
        p_particle->export_trajectory_as_csv(filename);
    }
    // std::cout << std::endl;
}

}  // namespace uepm::amc
