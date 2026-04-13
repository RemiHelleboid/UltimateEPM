/**
 * @file bulk_amc_simulation.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-04-10
 *
 *
 */

#pragma once

#include <random>
#include <vector>

#include "intervalley_phonon.hpp"
#include "particle_amc.hpp"
#include "valley_model.hpp"
#include "vector.hpp"

namespace uepm::amc {

struct bulk_amc_simulation_config {
    std::size_t            m_number_of_particles      = 10000;
    uepm::mesh_bz::vector3 m_electric_field           = {0.0, 0.0, 0.0};  // V/m
    double                 m_lattice_temperature      = 300.0;            // K
    double                 m_final_time               = 5.0e-12;          // s
    double                 m_doping_concentration     = 1.0e16;           // m^-3
    double                 m_max_self_scattering_rate = 1.0e14;           // s^-1
    bool                   m_record_history           = true;
};

struct bulk_observables {
    double      electric_field_V_per_m  = 0.0;
    double      mean_velocity_x_m_per_s = 0.0;
    double      mean_kinetic_energy_eV  = 0.0;
    std::size_t sample_count            = 0;
};

class bulk_amc_simulation {
 private:
    bulk_amc_simulation_config m_cfg;
    std::vector<valley_model>  m_valleys;
    std::vector<particle_amc>  m_particles;
    std::mt19937_64            m_rng;

    std::vector<intervalley_phonon_branch> m_intervalley_branches;

    bulk_observables m_observables;

 public:
    bulk_amc_simulation() : m_rng(std::random_device{}()) {}

    explicit bulk_amc_simulation(const bulk_amc_simulation_config& cfg) : m_cfg(cfg), m_rng(std::random_device{}()) {}

    void initialize();
    void drift_particle(particle_amc& p, double dt);
    void scatter_particle(particle_amc& p, double dt);
    void run();
    void export_particles_history_to_csv(const std::string& prefix_name) const;
    void accumulate_observables();
    void export_observables_to_csv(const std::string& filename) const;
};

}  // namespace uepm::amc