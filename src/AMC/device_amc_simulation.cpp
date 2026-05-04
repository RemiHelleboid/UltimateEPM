/**
 * @file device_amc_simulation.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-05-04
 * 
 * 
 */

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

#include "device_amc_simulation.hpp"

namespace uepm::amc {

std::optional<double> history_device_amc::get_time_at_number_particle(std::size_t nb_particle_reached) const {
    std::vector<std::size_t> history_total_number_particles = get_history_total_nb_particles();
    auto                     it_nb_particle_reached = std::find_if(history_total_number_particles.begin(),
                                                                   history_total_number_particles.end(),
                                                                   [&](const auto nb_part) { return (nb_part >= nb_particle_reached); });
    if (it_nb_particle_reached == history_total_number_particles.end()) {
        return std::nullopt;
    }
    return m_list_times[std::distance(history_total_number_particles.begin(), it_nb_particle_reached)];
}

device_amc_simulation::device_amc_simulation(const device::device      &simulation_device,
                                             const options_device_amc &simulation_option,
                                             const std::string         &simulation_name,
                                             int                        seed_random_generator)
    : m_device(simulation_device),
      m_dimension(m_device.get_dimension()),
      m_simulation_options(simulation_option),
      m_simulation_name(simulation_name),
      m_iteration{0},
      m_random_generator{seed_random_generator},
      m_uniform_distribution{std::uniform_real_distribution<double>{0.0, 1.0}} {
    m_simulation_history.m_initial_seed_rng = seed_random_generator;
}

device_amc_simulation::device_amc_simulation(const device::device      &device_simulation,
                                               const options_device_amc &simulation_option,
                                               const std::string         &simulation_name,
                                               const mesh::vector3       &starting_position,
                                               std::size_t                number_electrons_start,
                                               std::size_t                number_holes_start,
                                               int                        seed_random_generator)
    : m_device(device_simulation),
      m_dimension(m_device.get_dimension()),
      m_simulation_options(simulation_option),
      m_simulation_name(simulation_name),
      m_iteration{0},
      m_random_generator{seed_random_generator},
      m_uniform_distribution{std::uniform_real_distribution<double>{0.0, 1.0}} {
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
    double time_step = m_simulation_options.m_time_step;
    for (std::size_t idx_electron = 0; idx_electron < number_electrons_start; ++idx_electron) {
        m_list_particles.push_back(std::make_unique<electron>(idx_electron, starting_position, m_time, time_step));
    }
    for (std::size_t idx_hole = 0; idx_hole < number_holes_start; ++idx_hole) {
        m_list_particles.push_back(std::make_unique<hole>(idx_hole, starting_position, m_time, time_step));
    }
    //  Setup the random number for RPLA and set containing elements for each particles.
    for (auto &p_particle : m_list_particles) {
        const double new_random_number = m_uniform_distribution(m_random_generator);
        p_particle->reset_random_path_length_number(new_random_number);
        p_particle->set_containing_element(first_element);
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
    // const std::size_t idx_particle       = type_of_particle == particle_type::electron ? get_number_electrons() : get_number_holes();
    const std::size_t idx_particle       = m_list_particles.size();
    const double      new_r_rpl_particle = m_uniform_distribution(m_random_generator);
    double            time_step          = m_simulation_options.m_time_step;
    if (type_of_particle == particle_type::electron) {
        m_list_particles.push_back(std::make_unique<electron>(idx_particle, location, 0.0, time_step));
    } else {
        m_list_particles.push_back(std::make_unique<hole>(idx_particle, location, 0.0, time_step));
    }
    m_list_particles.back()->reset_random_path_length_number(new_r_rpl_particle);
    m_list_particles.back()->set_containing_element(first_element);
    m_list_particles.back()->set_weight(weight);
}

void device_amc_simulation::add_particles_at_positions(const std::vector<mesh::vector3> &positions,
                                                        particle_type                     type_of_particle,
                                                        double                            weight) {
    // Reserve memory for all particles upfront
    m_list_particles.reserve(m_list_particles.size() + positions.size());
    double time_step = m_simulation_options.m_time_step;
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
        const double      new_r_rpl_particle = m_uniform_distribution(m_random_generator);
        if (type_of_particle == particle_type::electron) {
            m_list_particles.push_back(std::make_unique<electron>(idx_particle, location, 0.0, time_step));
        } else {
            m_list_particles.push_back(std::make_unique<hole>(idx_particle, location, 0.0, time_step));
        }
        m_list_particles[m_list_particles.size() - 1]->reset_random_path_length_number(new_r_rpl_particle);
        m_list_particles[m_list_particles.size() - 1]->set_containing_element(first_element);
        m_list_particles[m_list_particles.size() - 1]->set_weight(weight);
    }
}

std::size_t device_admc_simulation::get_number_electrons() const {
    std::size_t nb_electron =
        std::accumulate(m_list_particles.begin(), m_list_particles.end(), 0, [](const std::size_t nb_part, const auto &p_part) {
            return nb_part + static_cast<std::size_t>(p_part->get_type() == particle_type::electron);
        });
    return nb_electron;
}

std::size_t device_admc_simulation::get_number_holes() const {
    std::size_t nb_electron =
        std::accumulate(m_list_particles.begin(), m_list_particles.end(), 0, [](const std::size_t nb_part, const auto &p_part_2) {
            return nb_part + static_cast<std::size_t>(p_part_2->get_type() == particle_type::hole);
        });
    return nb_electron;
}

double device_admc_simulation::compute_ramo_current() const {
    double current = 0.0;
    for (const auto &p_particle : m_list_particles) {
        current += p_particle->get_charge_sign() * p_particle->get_drift_velocity().dot(p_particle->get_electric_field());
    }

    current /= static_cast<double>(m_list_particles.size());
    return current;
}

void device_admc_simulation::perform_drift_diffusion_transport() {
    for (std::size_t idx_particle = 0; idx_particle < m_list_particles.size(); ++idx_particle) {
        m_list_particles[idx_particle]->set_data_from_device(m_dimension);

        const mesh::vector3 random_gaussian_vector(m_normal_distribution(m_random_generator),
                                                   m_normal_distribution(m_random_generator),
                                                   m_normal_distribution(m_random_generator));
        m_list_particles[idx_particle]->perform_transport_step(m_simulation_options.m_time_step,*
                                                               random_gaussian_vector,
                                                               m_device.get_lattice_temperature());
    }
}

void device_admc_simulation::set_particles_transport_data_from_device() {
    for (auto &p_particle : m_list_particles) {
        p_particle->set_data_from_device(m_dimension);
        p_particle->set_transport_data_from_device( m_device.get_lattice_temperature());
    }
}

void device_admc_simulation::perform_drift_diffusion_transport_without_data_setting() {
    for (std::size_t idx_particle = 0; idx_particle < m_list_particles.size(); ++idx_particle) {
        const mesh::vector3 random_gaussian_vector(m_normal_distribution(m_random_generator),
                                                   m_normal_distribution(m_random_generator),
                                                   m_normal_distribution(m_random_generator));
        m_list_particles[idx_particle]->perform_transport_step_with_current_data(random_gaussian_vector);
    }
}

void device_admc_simulation::update_element_and_check_boundary() {
    const bool is_2d = m_dimension == 2;
    for (auto &p_particle : m_list_particles) {
        const mesh::element *old_element = p_particle->get_containing_element();
        // const mesh::vector3  current_position = m_dimension == 2 ? p_particle->get_position().to_2d() : p_particle->get_position();
        mesh::vector3 current_position = p_particle->get_position();
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
            p_particle->reset_to_previous_position();
            continue;
        } else if (m_device.get_material_name_at_element(new_element) != "Silicon") {
            p_particle->reset_to_previous_position();
            p_particle->reset_cumulative_impact_ionization();
            continue;
        } else {
            p_particle->set_containing_element(new_element);
        }
    }
}

void device_admc_simulation::remove_collected_particles() {
    m_anode_current      = 0.0;
    m_cathode_current    = 0.0;
    bool remove_particle = false;
    int  nb_part_erased  = 0;
    for (const auto &p_particle : m_list_particles) {
        if (p_particle->get_crossed_contact()) {
            remove_particle = true;
            if (m_simulation_options.m_export_time_step) {
                // std::string filename = fmt::format("{}/particle_{:06}_trajectory.csv", m_simulation_name, p_particle->get_index());
                // p_particle->export_trajectory_as_csv(filename);
            }
            nb_part_erased++;
        }
    }
    if (remove_particle) {
        std::erase_if(m_list_particles, [](auto &&p_part) { return p_part->get_crossed_contact(); });
    }
    std::vector<double> currents = m_device.get_electrode_currents();
    m_anode_current              = physic::constant::elementary_charge * currents[0] / m_simulation_options.m_time_step;
    m_cathode_current            = physic::constant::elementary_charge * currents[1] / m_simulation_options.m_time_step;
    // std::cout << "Anode current : " << m_anode_current << "   Cathode current : " << m_cathode_current << std::endl;
    // std::cout << "Nb particles erased : " << nb_part_erased << std::endl;
}

void device_admc_simulation::perform_impact_ionization(double weight_poisson_new_particle) {
    std::size_t number_particles = m_list_particles.size();
    for (std::size_t idx_particle = 0; idx_particle < number_particles; ++idx_particle) {
        m_list_particles[idx_particle]->perform_impact_ionization_step(m_simulation_options.m_time_step,
                                                                       m_simulation_options.m_impact_ionization_model);
        if (m_list_particles[idx_particle]->has_impact_ionized()) {
            m_number_impact_ionizations++;
            m_list_particles[idx_particle]->reset_cumulative_impact_ionization();
            m_list_particles[idx_particle]->add_impact_ionization_to_history();
            m_simulation_history.m_impact_ionization_positions.push_back(m_list_particles[idx_particle]->get_position());
            m_simulation_history.m_last_impact_ionization_position = m_list_particles[idx_particle]->get_position();
            std::size_t index_new_particles                        = m_list_particles.size();
            double      new_r_parent_particle                      = m_uniform_distribution(m_random_generator);
            m_list_particles[idx_particle]->reset_random_path_length_number(new_r_parent_particle);
            if (m_simulation_options.m_particle_creation_activated) {
                double                    new_r_electron = m_uniform_distribution(m_random_generator);
                double                    new_r_hole     = m_uniform_distribution(m_random_generator);
                std::unique_ptr<particle> new_electron   = std::make_unique<electron>(index_new_particles,
                                                                                      m_list_particles[idx_particle]->get_position(),
                                                                                      0.0,
                                                                                      m_simulation_options.m_time_step);
                std::unique_ptr<particle> new_hole       = std::make_unique<hole>(index_new_particles + 1,
                                                                                  m_list_particles[idx_particle]->get_position(),
                                                                                  0.0,
                                                                                  m_simulation_options.m_time_step);
                new_electron->reset_random_path_length_number(new_r_electron);
                new_hole->reset_random_path_length_number(new_r_hole);
                new_electron->set_containing_element(m_list_particles[idx_particle]->get_containing_element());
                new_hole->set_containing_element(m_list_particles[idx_particle]->get_containing_element());
                new_electron->set_poisson_weight(weight_poisson_new_particle);
                new_hole->set_poisson_weight(weight_poisson_new_particle);
                new_electron->set_data_from_device(m_dimension);
                new_electron->set_transport_data_from_device(m_simulation_options.m_mobility_model, m_device.get_lattice_temperature());
                new_hole->set_data_from_device(m_dimension);
                new_hole->set_transport_data_from_device(m_simulation_options.m_mobility_model, m_device.get_lattice_temperature());
                m_list_particles.push_back(std::move(new_electron));
                m_list_particles.push_back(std::move(new_hole));
            }
        }
    }
}

void device_admc_simulation::run_transport_simulation() {
    while (m_time <= m_simulation_options.m_t_max && !m_list_particles.empty() &&
           m_list_particles.size() <= m_simulation_options.m_max_number_particle && get_number_electrons() >= 1) {
        perform_drift_diffusion_transport();
        update_element_and_check_boundary();
        m_time += m_simulation_options.m_time_step;
        ++m_iteration;
        if (m_simulation_options.m_export_time_step && m_iteration % m_simulation_options.m_frequency_export_trajectory == 0) {
            export_current_time_step_as_csv(m_prefix_export_filename);
        }
        remove_collected_particles();
    }
}

void device_admc_simulation::run_transport_ionization_simulation() {
    bool populate_particle_history = m_simulation_options.m_keep_particles_history;
    add_data_to_history();
    if (m_simulation_options.m_export_time_step) {
        export_current_time_step_as_csv(m_prefix_export_filename);
    }
    while (m_time <= m_simulation_options.m_t_max && !m_list_particles.empty() &&
           m_list_particles.size() <= m_simulation_options.m_max_number_particle && (get_number_electrons() > 0)) {
        perform_drift_diffusion_transport();
        update_element_and_check_boundary();
        perform_impact_ionization();
        m_time += m_simulation_options.m_time_step;
        ++m_iteration;
        // std::cout << "\r Nb particles : " << m_list_particles.size() << "    " << std::flush;
        if (m_simulation_options.m_export_time_step && m_iteration % m_simulation_options.m_frequency_export_trajectory == 0) {
            // std::cout << "Exporting iteration : " << m_iteration << std::endl;
            export_current_time_step_as_csv(m_prefix_export_filename);
        }
        remove_collected_particles();
        if (m_simulation_options.m_keep_particles_history) {
            add_data_to_history(populate_particle_history);
        }
    }
    if (m_list_particles.size() >= m_simulation_options.m_max_number_particle) {
        m_time_to_avalanche = m_time;
    }
    // std::cout << "Simulation finished at time : " << m_time << std::endl;
    // std::cout << "Nb particles : " << m_list_particles.size() << std::endl;
}

bool device_admc_simulation::run_transport_to_doping_threshold(const double doping_threshold) {
    if (m_list_particles.size() != 1) {
        throw std::invalid_argument("ERROR: This method is only suitable for a single initial particle (can be electron or hole)");
    }
    add_data_to_history();
    if (m_simulation_options.m_export_time_step) {
        export_current_time_step_as_csv(m_prefix_export_filename);
    }
    while (m_time <= m_simulation_options.m_t_max) {
        // std::cout << "INFO: Time is " << m_time << ", particle position is " << m_list_particles[0]->get_position() << std::endl;
        perform_drift_diffusion_transport();
        update_element_and_check_boundary();
        m_time += m_simulation_options.m_time_step;
        ++m_iteration;
        if (m_simulation_options.m_export_time_step && m_iteration % m_simulation_options.m_frequency_export_trajectory == 0) {
            export_current_time_step_as_csv(m_prefix_export_filename);
        }

        remove_collected_particles();
        if (m_simulation_options.m_keep_particles_history) {
            std::cout << "Adding data to history" << std::endl;
            add_data_to_history();
        }

        if (m_list_particles.empty()) {
            return true;
        }
        if (m_list_particles[0]->get_doping_concentration() > doping_threshold) {
            return true;
        }
    }
    return false;
}

std::vector<mesh::vector3> device_admc_simulation::get_all_particles_position() const {
    std::vector<mesh::vector3> all_positions(m_list_particles.size());
    std::transform(m_list_particles.begin(), m_list_particles.end(), all_positions.begin(), [](auto &&p_particle) {
        return p_particle->get_position();
    });
    return all_positions;
}

std::vector<mesh::vector3> device_admc_simulation::get_all_global_velocities() const {
    std::vector<mesh::vector3> all_velocities(m_list_particles.size());
    std::transform(m_list_particles.begin(), m_list_particles.end(), all_velocities.begin(), [](auto &&p_particle) {
        return p_particle->compute_global_velocities();
    });
    return all_velocities;
}

std::vector<std::size_t> device_admc_simulation::get_all_number_impact_ionization() const {
    std::vector<std::size_t> all_number_impact_ionization(m_list_particles.size());
    std::transform(m_list_particles.begin(), m_list_particles.end(), all_number_impact_ionization.begin(), [](auto &&p_particle) {
        return p_particle->get_total_number_impact_ionization();
    });
    return all_number_impact_ionization;
}

std::vector<mesh::vector3> device_admc_simulation::get_all_positions_impact_ionization() const {
    return m_simulation_history.m_impact_ionization_positions;
}

std::vector<double> device_admc_simulation::get_all_coefficient_impact_ionization() const {
    std::vector<double> all_number_impact_ionization(m_list_particles.size());
    std::transform(m_list_particles.begin(), m_list_particles.end(), all_number_impact_ionization.begin(), [](auto &&p_particle) {
        return p_particle->get_total_number_impact_ionization() / double(p_particle->get_displacement_x());
    });
    return all_number_impact_ionization;
}

std::vector<double> device_admc_simulation::get_all_mean_dead_spaces() const {
    std::vector<double> all_mean_dead_spaces(m_list_particles.size());
    std::transform(m_list_particles.begin(), m_list_particles.end(), all_mean_dead_spaces.begin(), [](auto &&p_particle) {
        return p_particle->get_mean_dead_spaces();
    });
    return all_mean_dead_spaces;
}

std::vector<mesh::vector3> device_admc_simulation::get_all_first_impact_ionization_position() const {
    std::vector<mesh::vector3> all_impact_ionization_position;
    for (const auto &p_particle : m_list_particles) {
        auto impact_ionization_position = p_particle->get_first_impact_ionization_position();
        if (impact_ionization_position.has_value()) {
            mesh::vector3 first_i_position = impact_ionization_position.value();
            all_impact_ionization_position.push_back(first_i_position);
        }
    }
    return all_impact_ionization_position;
}

std::pair<double, double> device_admc_simulation::compute_depletion_region() const {
    double x_min = std::numeric_limits<double>::max();
    double x_max = std::numeric_limits<double>::max();
    // Compute the distance of the particles to the center of the device. xmin is the maximal distance to the center towards the anode
    // and xmax is the maximal distance to the center towards the cathode.
    double center_x = m_device.get_p_mesh()->get_bounding_box().get_x_min() +
                      0.5 * (m_device.get_p_mesh()->get_bounding_box().get_x_max() - m_device.get_p_mesh()->get_bounding_box().get_x_min());
    for (const auto &p_particle : m_list_particles) {
        double distance_to_center = p_particle->get_position().x() - center_x;
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

void device_admc_simulation::export_current_time_step_as_csv(const std::string &prefix_filename) const {
    std::vector<double> x_positions;
    std::vector<double> y_positions;
    std::vector<double> z_positions;
    std::vector<double> particles_times;
    std::vector<double> particle_electric_field;
    std::vector<double> particles_impact_ionization_probability;
    std::vector<double> particles_charge;
    for (auto &&p_particle : m_list_particles) {
        x_positions.push_back(p_particle->get_position().x());
        y_positions.push_back(p_particle->get_position().y());
        z_positions.push_back(p_particle->get_position().z());
        particles_times.push_back(p_particle->get_time());
        particle_electric_field.push_back(p_particle->get_electric_field().norm());
        particles_impact_ionization_probability.push_back(p_particle->get_probability_impact_ionization());
        particles_charge.push_back(p_particle->get_charge_sign());
    }
    std::vector<std::string> list_column_names  = {"time", "X", "Y", "Z", "electric field", "impact ionization probability", "type"};
    std::string              iteration_filename = fmt::format("{}.{:09d}.csv", prefix_filename, m_iteration);
    utils::export_multiple_vector_to_csv(iteration_filename,
                                         list_column_names,
                                         {particles_times,
                                          x_positions,
                                          y_positions,
                                          z_positions,
                                          particle_electric_field,
                                          particles_impact_ionization_probability,
                                          particles_charge});
}

void device_admc_simulation::export_all_trajectories_as_csv(const std::string &prefix_filename) const {
    for (auto &&p_particle : m_list_particles) {
        // std::cout << "\rExporting trajectory of particle " << p_particle->get_index() << std::flush;
        std::string filename = fmt::format("{}particle_{:06}_trajectory.csv", prefix_filename, p_particle->get_index());
        p_particle->export_trajectory_as_csv(filename);
    }
    // std::cout << std::endl;
}

}  // namespace admc
