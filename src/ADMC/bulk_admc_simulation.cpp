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
#include <limits>
#include <stdexcept>
#include <utility>

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

    for (std::size_t i = 0; i < m_config.number_electrons; ++i) {
        m_particles.emplace_back(m_particles.size(), carrier_type::electron, m_config.initial_position_m);
    }
    for (std::size_t i = 0; i < m_config.number_holes; ++i) {
        m_particles.emplace_back(m_particles.size(), carrier_type::hole, m_config.initial_position_m);
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
        }
        m_current_time_s += time_step_s;
    }
}

}  // namespace uepm::ADMC
