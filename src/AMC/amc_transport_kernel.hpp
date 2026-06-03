/**
 * @file amc_transport_kernel.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-26
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <cstddef>
#include <optional>
#include <random>
#include <vector>

#include "amc_material_model.hpp"
#include "intervalley_phonon.hpp"
#include "particle_amc.hpp"
#include "scattering_channels.hpp"
#include "valley_model.hpp"
#include "vector.hpp"

namespace uepm::amc {

struct amc_transport_config {
    particle_type m_carrier_type                  = particle_type::electron;
    double        m_lattice_temperature           = 300.0;
    double        m_max_energy_eV                 = 2.0;
    double        m_self_scattering_safety_factor = 1.2;
    std::size_t   m_gamma_max_energy_samples      = 1000;

    double m_background_impurity_density_cm_3 = 0.0;

    bool m_enable_impact_ionization   = false;
    bool m_enable_impurity_scattering = false;
};

class amc_transport_kernel {
 public:
    amc_transport_kernel();
    explicit amc_transport_kernel(const amc_transport_config& cfg);
    amc_transport_kernel(const amc_transport_config& cfg, std::uint64_t seed);

    void                                ensure_gamma_max_covers(double total_rate);
    void                                initialize();
    const impact_ionization_parameters& impact_ionization_parameters_for_carrier() const;
    double                              impact_ionization_rate(double energy_eV) const;
    const std::vector<valley_model>&    valleys() const noexcept { return m_valleys; }
    double                              gamma_max() const noexcept { return m_gamma_max_s_1; }
    void                                initialize_particle_state(particle_amc& p);
    std::optional<scattering_event>     scatter_particle(particle_amc& p, double dt);
    void                            drift_particle(particle_amc& p, const mesh::vector3& electric_field_Vm, double dt);
    std::vector<scattering_channel> build_scattering_channels(const particle_amc& p) const;
    double                          total_scattering_rate(const particle_amc& p) const;
    double             total_scattering_rate_for_energy(std::size_t band_or_valley_index, double energy_eV) const;
    double             compute_max_self_scattering_rate(double max_energy_eV, std::size_t n_samples) const;
    double             sample_free_flight_time();
    scattering_channel select_scattering_channel(const particle_amc& p);
    scattering_event   apply_scattering_channel(particle_amc& p, const scattering_channel& channel);
    double             uniform01();

 private:
    amc_transport_config m_cfg;

    std::vector<valley_model>              m_valleys;
    std::vector<intervalley_phonon_branch> m_intervalley_branches;
    std::vector<hole_optical_transition>   m_hole_optical_transitions;

    carrier_impact_ionization_parameters m_impact_ionization_parameters;

    carrier_impurity_mobility_parameters m_impurity_mobility_parameters;

    std::mt19937_64 m_rng;
    double          m_gamma_max_s_1 = 0.0;
};

}  // namespace uepm::amc