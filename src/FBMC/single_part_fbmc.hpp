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

namespace fbmc {

struct bulk_environment {
    double   m_temperature;
    vector3 m_electric_field;
    double   m_doping_concentration;
};

struct simulation_parameters {
    double      m_simulation_time;
    std::size_t m_export_frequency;
};

class single_particle_simulation {
 private:
    bz_mesh*              m_ptr_mesh_bz;
    bulk_environment      m_bulk_env;
    simulation_parameters m_sim_params;
    // particle              m_particle;

    double m_time = 0.0;

 public:
    double a = 0;
    single_particle_simulation(mesh_bz* ptr_mesh_bz, const bulk_environment& bulk_env, const simulation_parameters& sim_params);

    void run_simulation();

    void export_history(const std::string& filename);
};
}