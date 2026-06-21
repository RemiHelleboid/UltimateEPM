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
#include <cstdint>
#include <string>
#include <vector>

#include "electron_phonon.hpp"
#include "keldysh_impactio.hpp"
#include "particle.hpp"
#include "vector_bz.hpp"

namespace uepm::fbmc {

struct Bulk_environment {
    double  m_temperature;           // K
    vector3 m_electric_field;        // V/m
    double  m_doping_concentration;  // cm^-3 (metadata only; impurity scattering is not implemented)
};

struct Simulation_parameters {
    double                  m_simulation_time;
    double                  m_warmup_fraction               = 0.2;
    std::size_t             m_nb_openmp_threads             = 1;
    double                  m_max_energy_eV                 = 10.0;
    double                  m_self_scattering_safety_factor = 1.2;
    bool                    m_enable_impact_ionization      = true;
    KeldyshImpactIonization m_impact_ionization_model;
    bool                    m_record_history = false;
    std::uint64_t           m_random_seed    = 1234;
    std::string             m_history_export_prefix;
    particle_type           m_particle_type = particle_type::electron;
};

struct bulk_fbmc_simulation_config {
    std::size_t m_number_of_particles = 10000;
    vector3     m_electric_field{0.0, 0.0, 0.0};  // V/m
    double      m_lattice_temperature       = 300.0;
    double      m_final_time                = 5.0e-12;
    double      m_doping_concentration_cm_3 = 0.0;

    double      m_warmup_fraction               = 0.2;
    double      m_max_energy_eV                 = 10.0;
    double      m_self_scattering_safety_factor = 1.2;
    std::size_t m_nb_threads                    = 1;

    bool                    m_enable_impact_ionization = false;
    KeldyshImpactIonization m_impact_ionization_model;
    bool                    m_record_history = false;
    std::uint64_t           m_random_seed    = 1234;
    std::string             m_history_export_prefix;
    particle_type           m_particle_type = particle_type::electron;
};

struct impact_ionization_coefficient_statistics {
    std::size_t m_events                         = 0;
    double      m_carrier_time_s                 = 0.0;
    double      m_drift_velocity_time_integral_m = 0.0;
    double      m_endpoint_displacement_m        = 0.0;

    double event_rate_per_carrier_s_1() const;
    double average_drift_velocity_m_per_s() const;
    double ionization_coefficient_cm_1() const;
    double endpoint_displacement_coefficient_cm_1() const;
};

struct bulk_observables {
    double m_electric_field_V_per_m       = 0.0;
    double m_weighted_velocity_x_m        = 0.0;
    double m_weighted_velocity_y_m        = 0.0;
    double m_weighted_velocity_z_m        = 0.0;
    double m_weighted_kinetic_energy_eV_s = 0.0;
    double m_accumulated_time_s           = 0.0;
};

class Single_particle_simulation {
 private:
    uepm::mesh_bz::ElectronPhonon* m_ptr_mesh_bz;
    Bulk_environment               m_bulk_env;
    Simulation_parameters          m_sim_params;
    std::vector<particle>          m_list_particle;
    std::size_t                    m_nb_particles;

    bulk_observables                         m_observables;
    impact_ionization_coefficient_statistics m_impact_ionization_statistics;
    double                                   m_gamma_max_s_1 = 0.0;
    std::size_t                              m_discarded_carriers_over_max_energy = 0;

 public:
    Single_particle_simulation(uepm::mesh_bz::ElectronPhonon* ptr_mesh_bz,
                               const Bulk_environment&        bulk_env,
                               const Simulation_parameters&   sim_params,
                               std::size_t                    nb_particles = 1);
    Single_particle_simulation(uepm::mesh_bz::ElectronPhonon* ptr_mesh_bz, const bulk_fbmc_simulation_config& config);

    std::size_t get_nb_particles() const { return m_nb_particles; }

    void run_simulation();
    void run() { run_simulation(); }

    void export_history(const std::string& filename);
    void extract_stats_and_export(const std::string& filename);
    void export_observables_to_csv(const std::string& filename) const;

    const bulk_observables&                         observables() const noexcept { return m_observables; }
    const impact_ionization_coefficient_statistics& impact_ionization_statistics() const noexcept {
        return m_impact_ionization_statistics;
    }
    std::size_t discarded_carriers_over_max_energy() const noexcept {
        return m_discarded_carriers_over_max_energy;
    }
    double gamma_max() const noexcept { return m_gamma_max_s_1; }
};

using bulk_fbmc_simulation = Single_particle_simulation;

}  // namespace uepm::fbmc
