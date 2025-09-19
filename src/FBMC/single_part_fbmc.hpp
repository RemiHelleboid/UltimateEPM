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

#include "bz_mesh.hpp"
#include "bz_states.hpp"
#include "electron_phonon.hpp"
#include "particle.hpp"
#include "vector.hpp"

namespace fbmc {

struct Bulk_environment {
    double  m_temperature;
    vector3 m_electric_field;
    double  m_doping_concentration;
};

struct Simulation_parameters {
    double      m_simulation_time;
    std::size_t m_export_frequency;
};

class Single_particle_simulation {
 private:
    bz_mesh::MeshBZ*      m_ptr_mesh_bz;
    Bulk_environment      m_bulk_env;
    Simulation_parameters m_sim_params;
    particle              m_particle;

    double m_time = 0.0;

 public:
    double a = 0;
    Single_particle_simulation(bz_mesh::MeshBZ* ptr_mesh_bz, const Bulk_environment& bulk_env, const Simulation_parameters& sim_params);

    void run_simulation();

    void export_history(const std::string& filename);
};
}  // namespace fbmc