/**
 * @file bulk_pbmc_simulation.hpp
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
#include "pbmc_material_model.hpp"
#include "pbmc_particle.hpp"
#include "pbmc_scattering_model.hpp"
#include "pbmc_transport_kernel.hpp"
#include "directional_diffusion.hpp"
#include "scattering_channels.hpp"
#include "valley_model.hpp"
#include "vector.hpp"

namespace uepm::PBMC {

struct bulk_pbmc_simulation_config {
    pbmc_material_model m_material_model      = make_silicon_pbmc_material_model();
    particle_type       m_carrier_type        = particle_type::electron;
    std::size_t         m_number_of_particles = 10000;
    uepm::mesh::vector3 m_electric_field      = {0.0, 0.0, 0.0};  // V/m
    double              m_lattice_temperature = 300.0;            // K
    double              m_final_time          = 5.0e-12;          // s

    double                    m_doping_concentration       = 1.0e16;  // m^-3
    bool                      m_enable_impurity_scattering = false;
    impurity_scattering_model m_impurity_scattering_model  = impurity_scattering_model::screened_coulomb;
    impurity_screening_model  m_impurity_screening_model   = impurity_screening_model::debye_analytic;
    double                    m_impurity_density_cm_3      = 0.0;  // cm^-3, positive scattering center density

    double m_max_self_scattering_rate = 1.0e15;  // s^-1
    bool   m_record_history           = true;
    double m_time_step                = 5.0e-15;
    double m_warmup_fraction          = 0.2;

    bool m_enable_impact_ionization = true;

    double      m_max_energy_eV                 = 2.0;
    double      m_self_scattering_safety_factor = 1.2;
    std::size_t m_gamma_max_energy_samples      = 1000;

    std::size_t m_nb_threads = 1;
};

struct impact_ionization_coefficient_statistics {
    std::size_t m_events = 0;

    double m_carrier_time_s                               = 0.0;
    double m_drift_velocity_time_integral_m_per_s_times_s = 0.0;
    double m_sampling_time_s                              = 0.0;

    double m_raw_ii_coefficient_cm_1 = 0.0;

    double event_rate_per_carrier_s_1() const {
        if (m_carrier_time_s <= 0.0) {
            return 0.0;
        }
        return static_cast<double>(m_events) / m_carrier_time_s;
    }

    double average_drift_velocity_m_per_s() const {
        if (m_sampling_time_s <= 0.0) {
            return 0.0;
        }
        return m_drift_velocity_time_integral_m_per_s_times_s / m_sampling_time_s;
    }

    double ionization_coefficient_cm_1() const {
        constexpr double meter_per_second_to_centimeter_per_second = 1.0e2;

        const double drift_velocity_cm_per_s =
            average_drift_velocity_m_per_s() * meter_per_second_to_centimeter_per_second;

        if (drift_velocity_cm_per_s <= 0.0) {
            return 0.0;
        }

        return event_rate_per_carrier_s_1() / drift_velocity_cm_per_s;
    }
};

struct bulk_observables {
    double electric_field_V_per_m        = 0.0;
    double weighted_velocity_x_m2_per_s2 = 0.0;
    double weighted_kinetic_energy_eV_s  = 0.0;
    double accumulated_time_s            = 0.0;
    statistics::directional_diffusion diffusion{};
};

class bulk_pbmc_simulation {
 private:
    bulk_pbmc_simulation_config m_cfg;
    pbmc_transport_kernel       m_transport;
    std::vector<pbmc_particle>  m_particles;
    bulk_observables            m_observables;

    impact_ionization_coefficient_statistics m_impact_ionization_coefficient_statistics;

    static pbmc_transport_config make_transport_config(const bulk_pbmc_simulation_config& cfg);
    void update_directional_diffusion();

 public:
    bulk_pbmc_simulation() : m_transport(make_transport_config(m_cfg), m_cfg.m_material_model) {}
    explicit bulk_pbmc_simulation(const bulk_pbmc_simulation_config& cfg)
        : m_cfg(cfg),
          m_transport(make_transport_config(m_cfg), m_cfg.m_material_model) {}

    void initialize();
    void run();
    void run_self_scattering_emc();

    std::size_t count_scattering_events(scattering_event event) const;
    double      average_drift_velocity_along_field_m_per_s() const;
    void        export_particles_history_to_csv(const std::string& prefix_name) const;
    void        accumulate_observables(double dt);
    void        accumulate_particle_observables(const pbmc_particle& p, double dt);
    void        export_observables_to_csv(const std::string& filename) const;

    const bulk_observables&                         observables() const noexcept { return m_observables; }
    const std::vector<pbmc_particle>&               particles() const noexcept { return m_particles; }
    const impact_ionization_coefficient_statistics& impact_ionization_statistics() const noexcept {
        return m_impact_ionization_coefficient_statistics;
    }
    double gamma_max() const noexcept { return m_transport.gamma_max(); }
};

}  // namespace uepm::PBMC
