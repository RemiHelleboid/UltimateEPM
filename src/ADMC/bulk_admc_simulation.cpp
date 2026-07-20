/**
 * @file bulk_admc_simulation.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-15
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "bulk_admc_simulation.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <utility>

#include "unit_conversion.hpp"

namespace uepm::ADMC {

void bulk_admc_simulation_config::validate() const {
    environment.validate();
    if (number_holes > std::numeric_limits<std::size_t>::max() - number_electrons) {
        throw std::invalid_argument("bulk ADMC particle count overflows size_t");
    }
    if (number_electrons + number_holes == 0) {
        throw std::invalid_argument("bulk ADMC simulation requires at least one particle");
    }
    if (!std::isfinite(time_step_s) || time_step_s <= 0.0) {
        throw std::invalid_argument("bulk ADMC time step must be positive and finite");
    }
    if (!std::isfinite(final_time_s) || final_time_s < 0.0) {
        throw std::invalid_argument("bulk ADMC final time must be finite and non-negative");
    }
}

bulk_admc_simulation::bulk_admc_simulation(bulk_admc_simulation_config config)
    : bulk_admc_simulation(std::move(config), silicon_arora_canali_mobility{}) {}

bulk_admc_simulation::bulk_admc_simulation(bulk_admc_simulation_config   config,
                                           silicon_arora_canali_mobility mobility_model)
    : m_config(std::move(config)),
      m_transport(std::move(mobility_model)),
      m_random_generator(m_config.random_seed) {
    m_config.validate();
}

void bulk_admc_simulation::initialize() {
    m_particles.clear();
    m_particles.reserve(m_config.number_electrons + m_config.number_holes);
    m_random_generator.seed(m_config.random_seed);
    m_standard_normal.reset();
    m_current_time_s = 0.0;
    m_diffusion      = {};

    const auto number_of_steps = static_cast<std::size_t>(std::ceil(m_config.final_time_s / m_config.time_step_s));

    for (std::size_t i = 0; i < m_config.number_electrons; ++i) {
        m_particles.emplace_back(m_particles.size(), carrier_type::electron, m_config.initial_position_m);
        m_transport.initialize_particle_state(m_particles.back(), m_config.environment);
        if (m_config.record_history) {
            m_particles.back().history().reserve(number_of_steps + 1);
            m_particles.back().record_state();
        }
    }
    for (std::size_t i = 0; i < m_config.number_holes; ++i) {
        m_particles.emplace_back(m_particles.size(), carrier_type::hole, m_config.initial_position_m);
        m_transport.initialize_particle_state(m_particles.back(), m_config.environment);
        if (m_config.record_history) {
            m_particles.back().history().reserve(number_of_steps + 1);
            m_particles.back().record_state();
        }
    }
}

vector3 bulk_admc_simulation::draw_standard_normal() {
    return {m_standard_normal(m_random_generator),
            m_standard_normal(m_random_generator),
            m_standard_normal(m_random_generator)};
}

void bulk_admc_simulation::run() {
    initialize();
    while (m_current_time_s < m_config.final_time_s) {
        const double time_step_s = std::min(m_config.time_step_s, m_config.final_time_s - m_current_time_s);
        for (auto& particle : m_particles) {
            m_transport.step(particle, m_config.environment, time_step_s, draw_standard_normal());
            if (m_config.record_history) {
                particle.record_state();
            }
        }
        m_current_time_s += time_step_s;
    }

    m_diffusion =
        statistics::estimate_directional_diffusion(m_particles.size(), m_current_time_s, [&](std::size_t index) {
            const auto displacement = m_particles[index].state().position_m - m_config.initial_position_m;
            return std::array<double, 3>{displacement.x(), displacement.y(), displacement.z()};
        });
}

void bulk_admc_simulation::export_particles_history_to_csv(const std::string& prefix_name) const {
    for (const auto& particle : m_particles) {
        const std::string filename = prefix_name + "_particle_" + std::to_string(particle.index()) + ".csv";
        std::ofstream     file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("cannot open particle history output file: " + filename);
        }
        file << std::setprecision(std::numeric_limits<double>::max_digits10);

        file << "time,"
                "position_x,position_y,position_z,"
                "local_k_x,local_k_y,local_k_z,"
                "velocity_x,velocity_y,velocity_z,"
                "lattice_temperature_K,kinetic_energy,gamma,valley_index\n";

        for (const auto& snapshot : particle.history().snapshots()) {
            const vector3 position_um = snapshot.position_m * uepm::units::meter_to_micron;
            file << snapshot.time_s << ',' << position_um.x() << ',' << position_um.y() << ',' << position_um.z()
                 << ",,,," << snapshot.total_velocity_m_per_s.x() << ',' << snapshot.total_velocity_m_per_s.y() << ','
                 << snapshot.total_velocity_m_per_s.z() << ',' << snapshot.lattice_temperature_K << ",,,\n";
        }
    }
}

}  // namespace uepm::ADMC
