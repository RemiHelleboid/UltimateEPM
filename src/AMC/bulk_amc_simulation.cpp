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

#include <fmt/format.h>
#include <fmt/core.h>

#include <array>
#include <iostream>
#include <fstream>
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

double sample_thermal_energy_eV(double temperature_K, std::mt19937_64& rng) {
    const double                    kT_eV = uepm::constants::k_b_eV * temperature_K;
    std::gamma_distribution<double> dist(1.5, kT_eV);
    return dist(rng);
}

void bulk_amc_simulation::initialize() {
    m_valleys = make_silicon_delta_valleys();

    std::cout << "Initialized " << m_valleys.size() << " valleys\n";
    for (std::size_t i = 0; i < m_valleys.size(); ++i) {
        std::cout << "  valley[" << i << "] = " << m_valleys[i].name() << '\n';
    }

    m_particles.clear();
    m_particles.reserve(m_cfg.m_number_of_particles);

    for (std::size_t i = 0; i < m_cfg.m_number_of_particles; ++i) {
        particle_amc p{i, particle_type::electron, 1.0};

        p.state().time         = 0.0;
        p.state().position     = vector3{0.0, 0.0, 0.0};
        p.state().valley_index = i % m_valleys.size();

        const auto& valley = m_valleys[p.state().valley_index];

        const double energy_eV   = sample_thermal_energy_eV(m_cfg.m_lattice_temperature, m_rng);
        p.state().local_k        = valley.draw_random_k_valley_at_energy(energy_eV, m_rng);
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

void bulk_amc_simulation::drift_particle(particle_amc& p, double dt) {
    if (dt < 0.0) {
        throw std::invalid_argument("drift time step must be non-negative");
    }

    const auto valley_index = p.state().valley_index;
    if (valley_index >= m_valleys.size()) {
        throw std::out_of_range("invalid valley index in drift_particle");
    }

    const auto&   valley                = m_valleys[valley_index];
    const vector3 electric_field_valley = valley.to_valley_frame(m_cfg.m_electric_field);

    const double  prefactor    = p.signed_charge() / uepm::constants::h_bar;
    const vector3 old_velocity = p.state().velocity;

    p.state().local_k += electric_field_valley * (prefactor * dt);

    p.state().gamma          = valley.gamma_from_k_valley(p.state().local_k);
    p.state().kinetic_energy = valley.kinetic_energy_from_gamma(p.state().gamma);
    p.state().velocity       = valley.velocity_from_k_valley(p.state().local_k);

    const vector3 avg_velocity = 0.5 * (old_velocity + p.state().velocity);
    p.state().position += avg_velocity * dt;
    p.state().time += dt;
}
void bulk_amc_simulation::run() {
    if (m_particles.empty()) {
        throw std::runtime_error("simulation not initialized");
    }

    constexpr double  dt      = 1.0e-15;
    const std::size_t n_steps = static_cast<std::size_t>(std::ceil(m_cfg.m_final_time / dt));

    for (std::size_t step = 0; step < n_steps; ++step) {
        for (auto& p : m_particles) {
            drift_particle(p, dt);

            if (m_cfg.m_record_history) {
                p.record_state();
            }
        }
    }
    fmt::print("Completed {} steps of {} particles\n", n_steps, m_particles.size());

    double  avg_energy = 0.0;
    vector3 avg_velocity{0.0, 0.0, 0.0};

    for (const auto& p : m_particles) {
        avg_energy += p.state().kinetic_energy;
        avg_velocity += p.state().velocity;
    }

    avg_energy /= static_cast<double>(m_particles.size());
    avg_velocity /= static_cast<double>(m_particles.size());

    fmt::print("Average kinetic energy: {:.6f} eV\n", avg_energy);
    fmt::print("Average velocity: ({:.6e}, {:.6e}, {:.6e}) m/s\n", avg_velocity.x(), avg_velocity.y(), avg_velocity.z());
}

void bulk_amc_simulation::export_particles_history_to_csv(const std::string& prefix_name) const {
    // Implementation for exporting particle history to CSV
    for (const auto& p : m_particles) {
        const auto& history = p.history();
        const std::string filename = fmt::format("{}_particle_{}.csv", prefix_name, p.index());
        std::ofstream    file(filename);
        if (!file.is_open()) {
            fmt::print(stderr, "Failed to open file for writing: {}\n", filename);
            continue;
        }

        // Write CSV header
        file << "time,position_x,position_y,position_z,local_k_x,local_k_y,local_k_z,velocity_x,velocity_y,velocity_z,kinetic_energy,gamma,valley_index\n";

        // Write particle history
        for (const auto& snapshot : history.snapshots()) {
            file << fmt::format("{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
                                snapshot.time,
                                snapshot.position.x(),
                                snapshot.position.y(),
                                snapshot.position.z(),
                                snapshot.local_k.x(),
                                snapshot.local_k.y(),
                                snapshot.local_k.z(),
                                snapshot.velocity.x(),
                                snapshot.velocity.y(),
                                snapshot.velocity.z(),
                                snapshot.kinetic_energy,
                                snapshot.gamma,
                                snapshot.valley_index);
        }
    }
}

}  // namespace uepm::amc