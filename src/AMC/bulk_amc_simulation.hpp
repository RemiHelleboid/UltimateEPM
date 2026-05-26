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

#include "amc_material_model.hpp"
#include "amc_scattering_model.hpp"
#include "amc_transport_kernel.hpp"
#include "intervalley_phonon.hpp"
#include "particle_amc.hpp"
#include "scattering_channels.hpp"
#include "valley_model.hpp"
#include "vector.hpp"

namespace uepm::amc {

struct bulk_amc_simulation_config {
    particle_type       m_carrier_type             = particle_type::electron;
    std::size_t         m_number_of_particles      = 10000;
    uepm::mesh::vector3 m_electric_field           = {0.0, 0.0, 0.0};  // V/m
    double              m_lattice_temperature      = 300.0;            // K
    double              m_final_time               = 5.0e-12;          // s
    double              m_doping_concentration     = 1.0e16;           // m^-3
    double              m_max_self_scattering_rate = 1.0e15;           // s^-1
    bool                m_record_history           = true;
    double              m_time_step                = 5.0e-15;
    double              m_warmup_fraction          = 0.2;

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
    amc_transport_kernel       m_transport;
    std::vector<particle_amc>  m_particles;
    bulk_observables           m_observables;

    static amc_transport_config make_transport_config(const bulk_amc_simulation_config& cfg);

 public:
    bulk_amc_simulation() : m_transport(make_transport_config(m_cfg)) {}
    explicit bulk_amc_simulation(const bulk_amc_simulation_config& cfg) : m_cfg(cfg), m_transport(make_transport_config(m_cfg)) {}

    void initialize();
    void run();
    void run_self_scattering_emc();

    void export_particles_history_to_csv(const std::string& prefix_name) const;
    void accumulate_observables(double dt);
    void accumulate_particle_observables(const particle_amc& p, double dt);
    void export_observables_to_csv(const std::string& filename) const;
};

}  // namespace uepm::amc