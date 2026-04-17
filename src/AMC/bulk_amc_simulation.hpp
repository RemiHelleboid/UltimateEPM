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
#include "scattering_channels.hpp"
#include "valley_model.hpp"
#include "vector.hpp"

namespace uepm::amc {

struct hole_optical_transition {
    const char* name                           = "";
    std::size_t initial_band                   = 0;
    std::size_t final_band                     = 0;
    double      phonon_energy_eV               = 0.0;
    double      deformation_potential_eV_per_m = 0.0;
    double      overlap_factor                 = 1.0;
};

struct bulk_amc_simulation_config {
    particle_type          m_carrier_type             = particle_type::electron;
    std::size_t            m_number_of_particles      = 10000;
    uepm::mesh_bz::vector3 m_electric_field           = {0.0, 0.0, 0.0};  // V/m
    double                 m_lattice_temperature      = 300.0;            // K
    double                 m_final_time               = 5.0e-12;          // s
    double                 m_doping_concentration     = 1.0e16;           // m^-3
    double                 m_max_self_scattering_rate = 1.0e15;           // s^-1
    bool                   m_record_history           = true;
    double                 m_time_step                = 5.0e-15;
    double                 m_warmup_fraction          = 0.2;

    double      m_max_energy_eV                 = 2.0;
    double      m_self_scattering_safety_factor = 1.2;
    std::size_t m_gamma_max_energy_samples      = 1000;
};

struct bulk_observables {
    double electric_field_V_per_m        = 0.0;
    double weighted_velocity_x_m2_per_s2 = 0.0;
    double weighted_kinetic_energy_eV_s  = 0.0;
    double accumulated_time_s            = 0.0;
};

class bulk_amc_simulation {
 private:
    bulk_amc_simulation_config m_cfg;
    std::vector<valley_model>  m_valleys;
    std::vector<particle_amc>  m_particles;
    std::mt19937_64            m_rng;

    std::vector<intervalley_phonon_branch> m_intervalley_branches;
    std::vector<hole_optical_transition>   m_hole_optical_transitions;
    bulk_observables                       m_observables;

    double m_gamma_max_s_1 = 0.0;

 public:
    bulk_amc_simulation() : m_rng(std::random_device{}()) {}
    explicit bulk_amc_simulation(const bulk_amc_simulation_config& cfg) : m_cfg(cfg), m_rng(std::random_device{}()) {}

    double total_scattering_rate_for_energy(std::size_t valley_index, double energy_eV) const;
    double compute_max_self_scattering_rate(double max_energy_eV, std::size_t n_samples) const;

    void initialize();
    void drift_particle(particle_amc& p, double dt);
    void scatter_particle(particle_amc& p, double dt);
    void run();

    std::vector<scattering_channel> build_scattering_channels(const particle_amc& p) const;
    double                          total_scattering_rate(const particle_amc& p) const;
    void                            apply_scattering_channel(particle_amc& p, const scattering_channel& channel);

    double             sample_free_flight_time();
    scattering_channel select_scattering_channel(const particle_amc& p);
    void               run_self_scattering_emc();

    void export_particles_history_to_csv(const std::string& prefix_name) const;
    void accumulate_observables(double dt);
    void accumulate_particle_observables(const particle_amc& p, double dt);
    void export_observables_to_csv(const std::string& filename) const;
};

}  // namespace uepm::amc