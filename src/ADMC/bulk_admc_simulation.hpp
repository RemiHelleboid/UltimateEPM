/**
 * @file bulk_admc_simulation.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-06-15
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "admc_transport.hpp"

namespace uepm::ADMC {

struct bulk_admc_simulation_config {
    admc_local_environment environment{};
    vector3                initial_position_m{};
    std::size_t            number_electrons = 1;
    std::size_t            number_holes     = 0;
    double                 time_step_s       = 1.0e-15;
    double                 final_time_s      = 1.0e-12;
    std::uint64_t          random_seed       = 5489u;

    void validate() const;
};

class bulk_admc_simulation {
 public:
    explicit bulk_admc_simulation(bulk_admc_simulation_config config);
    bulk_admc_simulation(bulk_admc_simulation_config config, silicon_arora_canali_mobility mobility_model);

    void initialize();
    void run();

    const bulk_admc_simulation_config& config() const noexcept { return m_config; }
    const std::vector<admc_particle>&  particles() const noexcept { return m_particles; }
    double                             current_time_s() const noexcept { return m_current_time_s; }

 private:
    vector3 draw_standard_normal();

    bulk_admc_simulation_config m_config;
    admc_transport_kernel       m_transport;
    std::vector<admc_particle>  m_particles;
    std::mt19937_64             m_random_generator;
    std::normal_distribution<double> m_standard_normal{0.0, 1.0};
    double                      m_current_time_s = 0.0;
};

}  // namespace uepm::ADMC
