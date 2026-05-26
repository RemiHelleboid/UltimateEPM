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

#include <fmt/core.h>
#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace uepm::amc {

amc_transport_config bulk_amc_simulation::make_transport_config(const bulk_amc_simulation_config& cfg) {
    amc_transport_config transport_cfg;

    transport_cfg.m_carrier_type                  = cfg.m_carrier_type;
    transport_cfg.m_lattice_temperature           = cfg.m_lattice_temperature;
    transport_cfg.m_max_energy_eV                 = cfg.m_max_energy_eV;
    transport_cfg.m_self_scattering_safety_factor = cfg.m_self_scattering_safety_factor;
    transport_cfg.m_gamma_max_energy_samples      = cfg.m_gamma_max_energy_samples;

    return transport_cfg;
}

void bulk_amc_simulation::initialize() {
    m_transport.initialize();

    m_particles.clear();
    m_particles.reserve(m_cfg.m_number_of_particles);

    for (std::size_t i = 0; i < m_cfg.m_number_of_particles; ++i) {
        particle_amc p{i, m_cfg.m_carrier_type, 1.0};

        p.state().time         = 0.0;
        p.state().position     = vector3{0.0, 0.0, 0.0};
        p.state().valley_index = i % m_transport.valleys().size();

        m_transport.initialize_particle_state(p);

        if (m_cfg.m_record_history) {
            p.record_state();
        }

        m_particles.push_back(std::move(p));
    }

    fmt::print("Computed max self-scattering rate: {:.6e} s^-1\n", m_transport.gamma_max());

    std::cout << "Initialized " << m_particles.size() << " particles\n";
}

void bulk_amc_simulation::accumulate_observables(double dt) {
    if (dt <= 0.0) {
        return;
    }

    for (const auto& p : m_particles) {
        m_observables.weighted_velocity_x_m2_per_s2 += p.state().velocity.x() * dt;
        m_observables.weighted_kinetic_energy_eV_s += p.state().kinetic_energy * dt;
    }

    m_observables.accumulated_time_s += static_cast<double>(m_particles.size()) * dt;
}

void bulk_amc_simulation::accumulate_particle_observables(const particle_amc& p, double dt) {
    if (dt <= 0.0) {
        return;
    }

    m_observables.weighted_velocity_x_m2_per_s2 += p.state().velocity.x() * dt;
    m_observables.weighted_kinetic_energy_eV_s += p.state().kinetic_energy * dt;
    m_observables.accumulated_time_s += dt;
}

void bulk_amc_simulation::run() {
    if (m_particles.empty()) {
        throw std::runtime_error("simulation not initialized");
    }

    const double      dt           = m_cfg.m_time_step;
    const std::size_t n_steps      = static_cast<std::size_t>(std::ceil(m_cfg.m_final_time / dt));
    const std::size_t warmup_steps = static_cast<std::size_t>(m_cfg.m_warmup_fraction * n_steps);

    m_observables = {};

    for (std::size_t step = 0; step < n_steps; ++step) {
        for (auto& p : m_particles) {
            m_transport.drift_particle(p, m_cfg.m_electric_field, dt);
            m_transport.scatter_particle(p, dt);

            if (m_cfg.m_record_history) {
                p.record_state();
            }
        }

        if (step >= warmup_steps) {
            accumulate_observables(dt);
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

    std::size_t total_acoustic_events               = 0;
    std::size_t total_intervalley_absorption_events = 0;
    std::size_t total_intervalley_emission_events   = 0;
    std::size_t total_self_scattering_events        = 0;

    for (const auto& p : m_particles) {
        const auto& events = p.history().scattering_events();

        total_acoustic_events += events[static_cast<std::size_t>(scattering_event::acoustic)];
        total_intervalley_absorption_events += events[static_cast<std::size_t>(scattering_event::intervalley_absorption)];
        total_intervalley_emission_events += events[static_cast<std::size_t>(scattering_event::intervalley_emission)];
        total_self_scattering_events += events[static_cast<std::size_t>(scattering_event::self_scattering)];
    }

    fmt::print("Total acoustic events: {}\n", total_acoustic_events);
    fmt::print("Total intervalley absorption events: {}\n", total_intervalley_absorption_events);
    fmt::print("Total intervalley emission events: {}\n", total_intervalley_emission_events);
    fmt::print("Total self-scattering events: {}\n", total_self_scattering_events);

    m_observables.electric_field_V_per_m = m_cfg.m_electric_field.norm();

    if (m_observables.accumulated_time_s == 0.0) {
        throw std::runtime_error("no steady-state samples accumulated");
    }

    const double avg_vx_m_per_s = m_observables.weighted_velocity_x_m2_per_s2 / m_observables.accumulated_time_s;

    const double avg_energy_eV = m_observables.weighted_kinetic_energy_eV_s / m_observables.accumulated_time_s;

    fmt::print("Steady-state average vx: {:.6e} m/s\n", avg_vx_m_per_s);
    fmt::print("Steady-state average energy: {:.6f} eV\n", avg_energy_eV);
}

void bulk_amc_simulation::run_self_scattering_emc() {
    if (m_particles.empty()) {
        throw std::runtime_error("simulation not initialized");
    }

    if (m_cfg.m_final_time <= 0.0) {
        throw std::invalid_argument("final time must be > 0");
    }

    if (m_cfg.m_warmup_fraction < 0.0 || m_cfg.m_warmup_fraction >= 1.0) {
        throw std::invalid_argument("warmup fraction must be in [0, 1)");
    }

    if (m_transport.gamma_max() <= 0.0) {
        throw std::invalid_argument("max self-scattering rate must be > 0");
    }

    m_observables = {};

    const double warmup_time             = m_cfg.m_warmup_fraction * m_cfg.m_final_time;
    double       max_observed_total_rate = 0.0;

    for (auto& p : m_particles) {
        while (p.state().time < m_cfg.m_final_time) {
            const double tau            = m_transport.sample_free_flight_time();
            const double remaining_time = m_cfg.m_final_time - p.state().time;
            const double drift_time     = std::min(tau, remaining_time);

            m_transport.drift_particle(p, m_cfg.m_electric_field, drift_time);

            if (m_cfg.m_record_history) {
                p.record_state();
            }

            if (p.state().time >= warmup_time) {
                accumulate_particle_observables(p, drift_time);
            }

            if (drift_time < tau) {
                break;
            }

            const double total_rate = m_transport.total_scattering_rate(p);
            max_observed_total_rate = std::max(max_observed_total_rate, total_rate);

            m_transport.ensure_gamma_max_covers(total_rate);

            const double u = m_transport.uniform01();

            if (u < total_rate / m_transport.gamma_max()) {
                const auto channel = m_transport.select_scattering_channel(p);
                m_transport.apply_scattering_channel(p, channel);
            } else {
                p.increment_scattering_event_count();
                p.add_scattering_event(scattering_event::self_scattering);
            }

            if (m_cfg.m_record_history) {
                p.record_state();
            }
        }
    }

    fmt::print("Completed self-scattering EMC run with {} particles\n", m_particles.size());
    fmt::print("Maximum observed total scattering rate: {:.6e} s^-1\n", max_observed_total_rate);

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

    std::size_t total_acoustic_events               = 0;
    std::size_t total_intervalley_absorption_events = 0;
    std::size_t total_intervalley_emission_events   = 0;
    std::size_t total_self_scattering_events        = 0;

    for (const auto& p : m_particles) {
        const auto& events = p.history().scattering_events();

        total_acoustic_events += events[static_cast<std::size_t>(scattering_event::acoustic)];
        total_intervalley_absorption_events += events[static_cast<std::size_t>(scattering_event::intervalley_absorption)];
        total_intervalley_emission_events += events[static_cast<std::size_t>(scattering_event::intervalley_emission)];
        total_self_scattering_events += events[static_cast<std::size_t>(scattering_event::self_scattering)];
    }

    fmt::print("Total acoustic events: {}\n", total_acoustic_events);
    fmt::print("Total intervalley absorption events: {}\n", total_intervalley_absorption_events);
    fmt::print("Total intervalley emission events: {}\n", total_intervalley_emission_events);
    fmt::print("Total self-scattering events: {}\n", total_self_scattering_events);

    m_observables.electric_field_V_per_m = m_cfg.m_electric_field.norm();

    if (m_observables.accumulated_time_s == 0.0) {
        throw std::runtime_error("no steady-state samples accumulated");
    }

    const double avg_vx_m_per_s = m_observables.weighted_velocity_x_m2_per_s2 / m_observables.accumulated_time_s;

    const double avg_energy_eV = m_observables.weighted_kinetic_energy_eV_s / m_observables.accumulated_time_s;

    fmt::print("Steady-state average vx: {:.6e} m/s\n", avg_vx_m_per_s);
    fmt::print("Steady-state average energy: {:.6f} eV\n", avg_energy_eV);
}

/**
 * @brief Append the current observables to a CSV file. If the file does not exist, it will be created with a header. If it already exists,
 * a new line will be appended with the current observables values.
 *
 * @param filename
 */
void bulk_amc_simulation::export_observables_to_csv(const std::string& filename) const {
    std::ofstream file(filename, std::ios::app);
    if (!file.is_open()) {
        fmt::print(stderr, "Failed to open file for writing: {}\n", filename);
        return;
    }
    // Check if the file is empty to write the header
    if (file.tellp() == 0) {
        file << "electric_field_V_per_m,mean_velocity_x_m_per_s,mean_kinetic_energy_eV,sample_count\n";
    }
    file << fmt::format("{},{},{},{}\n",
                        m_observables.electric_field_V_per_m,
                        m_observables.weighted_velocity_x_m2_per_s2 / m_observables.accumulated_time_s,
                        m_observables.weighted_kinetic_energy_eV_s / m_observables.accumulated_time_s,
                        m_observables.accumulated_time_s);
}

void bulk_amc_simulation::export_particles_history_to_csv(const std::string& prefix_name) const {
    // Implementation for exporting particle history to CSV
    for (const auto& p : m_particles) {
        const auto&       history  = p.history();
        const std::string filename = fmt::format("{}_particle_{}.csv", prefix_name, p.index());
        std::ofstream     file(filename);
        if (!file.is_open()) {
            fmt::print(stderr, "Failed to open file for writing: {}\n", filename);
            continue;
        }

        // Write CSV header
        file << "time,position_x,position_y,position_z,local_k_x,local_k_y,local_k_z,velocity_x,velocity_y,velocity_z,kinetic_energy,gamma,"
                "valley_index\n";

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