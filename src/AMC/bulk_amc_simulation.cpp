/**
 * @file bulk_amc_simulation.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-04-10
 *
 *
 */

#include "bulk_amc_simulation.hpp"

#include <array>
#include <iostream>
#include <random>
#include <stdexcept>

namespace uepm::amc {

constexpr double m0 = uepm::constants::m_e;

// Local z -> global z
valley_model::mat3 rotation_local_z_to_global_z() { return valley_model::identity_matrix(); }

// Local z -> global x
valley_model::mat3 rotation_local_z_to_global_x() { return {{{{0.0, 0.0, 1.0}}, {{0.0, 1.0, 0.0}}, {{-1.0, 0.0, 0.0}}}}; }

// Local z -> global y
valley_model::mat3 rotation_local_z_to_global_y() { return {{{{1.0, 0.0, 0.0}}, {{0.0, 0.0, 1.0}}, {{0.0, -1.0, 0.0}}}}; }

std::vector<valley_model> make_silicon_delta_valleys() {
    std::vector<valley_model> valleys;
    valleys.reserve(6);

    valley_model::parameters p;
    p.transverse_effective_mass   = 0.19 * m0;
    p.longitudinal_effective_mass = 0.916 * m0;

    // Keep this simple for now.
    // Use parabolic first; switch to Kane once the loop works.
    p.dispersion = valley_model::band_type::parabolic;

    // If you want Kane immediately, set:
    // p.dispersion = valley_model::band_type::kane;
    // p.non_parabolicity = ...; // in 1/eV or 1/J depending on your convention

    p.non_parabolicity        = 0.0;
    p.energy_offset           = 0.0;
    p.phonon_reference_energy = 0.0;
    p.degeneracy              = 1;

    p.name     = "Delta_x_plus";
    p.rotation = rotation_local_z_to_global_x();
    valleys.emplace_back(p);

    p.name     = "Delta_x_minus";
    p.rotation = rotation_local_z_to_global_x();
    valleys.emplace_back(p);

    p.name     = "Delta_y_plus";
    p.rotation = rotation_local_z_to_global_y();
    valleys.emplace_back(p);

    p.name     = "Delta_y_minus";
    p.rotation = rotation_local_z_to_global_y();
    valleys.emplace_back(p);

    p.name     = "Delta_z_plus";
    p.rotation = rotation_local_z_to_global_z();
    valleys.emplace_back(p);

    p.name     = "Delta_z_minus";
    p.rotation = rotation_local_z_to_global_z();
    valleys.emplace_back(p);

    return valleys;
}

void bulk_amc_simulation::initialize() {
    m_valleys = make_silicon_delta_valleys();

    std::cout << "Initialized " << m_valleys.size() << " valleys\n";
    for (std::size_t i = 0; i < m_valleys.size(); ++i) {
        std::cout << "  valley[" << i << "] = " << m_valleys[i].name() << '\n';
    }

    for (std::size_t i = 0; i < m_cfg.m_number_of_particles; ++i) {
        particle_amc p{i, particle_type::electron, 1.0};

        p.state().time         = 0.0;
        p.state().position     = vector3{0.0, 0.0, 0.0};
        p.state().local_k      = vector3{0.0, 0.0, 0.0};
        p.state().valley_index = i % m_valleys.size();

        const auto& valley       = m_valleys[p.state().valley_index];
        p.state().gamma          = valley.gamma_from_k_valley(p.state().local_k);
        p.state().kinetic_energy = valley.kinetic_energy_from_gamma(p.state().gamma);
        p.state().velocity       = valley.velocity_from_k_valley(p.state().local_k);

        if (m_cfg.m_record_history) {
            p.record_state();
        }

        m_particles.push_back(std::move(p));
    }
    std::cout << "Initialized " << m_particles.size() << " particles\n";
}

}  // namespace uepm::amc