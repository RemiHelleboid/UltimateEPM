/**
 * @file single_part_fbmc.hpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2025-09-17
 *
 * @copyright Copyright (c) 2025
 *
 */

#pragma once

#include <array>
#include <cmath>
#include <cstddef>

#include "electron_phonon.hpp"
#include "particle.hpp"
#include "vector.hpp"

namespace uepm::fbmc {

struct Bulk_environment {
    double  m_temperature;
    vector3 m_electric_field;
    double  m_doping_concentration;
};

struct Simulation_parameters {
    double      m_simulation_time;
    std::size_t m_export_frequency  = 10;
    std::size_t m_nb_openmp_threads = 1;
};

class Single_particle_simulation {
 private:
    uepm::mesh_bz::ElectronPhonon* m_ptr_mesh_bz;
    Bulk_environment               m_bulk_env;
    Simulation_parameters          m_sim_params;
    std::vector<particle>          m_list_particle;
    std::size_t                    m_nb_particles;

    double      m_time = 0.0;
    std::size_t m_step = 0;

 public:
    double a = 0;
    Single_particle_simulation(uepm::mesh_bz::ElectronPhonon* ptr_mesh_bz,
                               const Bulk_environment&        bulk_env,
                               const Simulation_parameters&   sim_params,
                               std::size_t                    nb_particles = 1);

    std::size_t get_nb_particles() const { return m_nb_particles; }
    void        set_nb_particles(std::size_t nb_particles) { m_nb_particles = nb_particles; }

    void run_simulation();

    void export_history(const std::string& filename);

    void extract_stats_and_export(const std::string& filename);
};
}  // namespace uepm::fbmc