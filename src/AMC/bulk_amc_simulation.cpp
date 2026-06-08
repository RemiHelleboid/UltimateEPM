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
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace uepm::amc {

namespace {

std::size_t particle_scattering_event_count(const particle_amc& particle, scattering_event event) {
    const std::size_t event_index = static_cast<std::size_t>(event);
    return particle.history().scattering_events()[event_index];
}

double particle_drift_velocity_along_field_m_per_s(const particle_amc& particle, const mesh::vector3& electric_field) {
    const double field_norm = electric_field.norm();

    if (field_norm <= 0.0) {
        return 0.0;
    }

    const mesh::vector3 field_direction = electric_field / field_norm;
    return std::abs(particle.state().velocity.dot(field_direction));
}

double sampled_time_after_warmup(double time_before_s, double time_after_s, double warmup_time_s, double final_time_s) {
    const double sample_start = std::max(time_before_s, warmup_time_s);
    const double sample_stop  = std::min(time_after_s, final_time_s);

    return std::max(0.0, sample_stop - sample_start);
}

}  // namespace

amc_transport_config bulk_amc_simulation::make_transport_config(const bulk_amc_simulation_config& cfg) {
    amc_transport_config transport_cfg;

    transport_cfg.m_carrier_type                  = cfg.m_carrier_type;
    transport_cfg.m_lattice_temperature           = cfg.m_lattice_temperature;
    transport_cfg.m_max_energy_eV                 = cfg.m_max_energy_eV;
    transport_cfg.m_self_scattering_safety_factor = cfg.m_self_scattering_safety_factor;
    transport_cfg.m_gamma_max_energy_samples      = cfg.m_gamma_max_energy_samples;
    transport_cfg.m_enable_impact_ionization      = cfg.m_enable_impact_ionization;

    transport_cfg.m_enable_impurity_scattering = cfg.m_enable_impurity_scattering;
    transport_cfg.m_background_impurity_density_cm_3 =
        cfg.m_enable_impurity_scattering ? cfg.m_impurity_density_cm_3 : 0.0;
    transport_cfg.m_impurity_scattering_model = cfg.m_impurity_scattering_model;

    return transport_cfg;
}

std::size_t bulk_amc_simulation::count_scattering_events(scattering_event event) const {
    std::size_t total = 0;

    const std::size_t event_index = static_cast<std::size_t>(event);

    for (const auto& particle : m_particles) {
        total += particle.history().scattering_events()[event_index];
    }

    return total;
}

double bulk_amc_simulation::average_drift_velocity_along_field_m_per_s() const {
    const double field_norm = m_cfg.m_electric_field.norm();

    if (field_norm <= 0.0 || m_particles.empty()) {
        return 0.0;
    }

    const mesh::vector3 field_direction = m_cfg.m_electric_field / field_norm;

    double velocity_sum = 0.0;

    for (const auto& particle : m_particles) {
        velocity_sum += particle.state().velocity.dot(field_direction);
    }

    const double mean_parallel_velocity = velocity_sum / static_cast<double>(m_particles.size());

    return std::abs(mean_parallel_velocity);
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
    fmt::print("Initialized {} particles\n", m_particles.size());
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
    if (m_cfg.m_time_step <= 0.0) {
        throw std::invalid_argument("time step must be > 0");
    }
    if (m_cfg.m_final_time <= 0.0) {
        throw std::invalid_argument("final time must be > 0");
    }
    if (m_cfg.m_warmup_fraction < 0.0 || m_cfg.m_warmup_fraction >= 1.0) {
        throw std::invalid_argument("warmup fraction must be in [0, 1)");
    }

    const double      dt           = m_cfg.m_time_step;
    const std::size_t n_steps      = static_cast<std::size_t>(std::ceil(m_cfg.m_final_time / dt));
    const std::size_t warmup_steps = static_cast<std::size_t>(m_cfg.m_warmup_fraction * static_cast<double>(n_steps));

    m_observables                                 = {};
    m_impact_ionization_coefficient_statistics    = {};
    std::size_t previous_impact_ionization_events = count_scattering_events(scattering_event::impact_ionization);

    for (std::size_t step = 0; step < n_steps; ++step) {
        for (auto& p : m_particles) {
            m_transport.drift_particle(p, m_cfg.m_electric_field, dt);
            m_transport.scatter_particle(p, dt);
            if (m_cfg.m_record_history) {
                p.record_state();
            }
        }
        const std::size_t current_impact_ionization_events =
            count_scattering_events(scattering_event::impact_ionization);
        const std::size_t new_impact_ionization_events =
            current_impact_ionization_events - previous_impact_ionization_events;
        previous_impact_ionization_events = current_impact_ionization_events;
        const bool collect_statistics     = step >= warmup_steps;

        if (collect_statistics) {
            accumulate_observables(dt);
            auto& stats = m_impact_ionization_coefficient_statistics;
            stats.m_events += new_impact_ionization_events;
            stats.m_carrier_time_s += static_cast<double>(m_particles.size()) * dt;
            stats.m_drift_velocity_time_integral_m_per_s_times_s += average_drift_velocity_along_field_m_per_s() * dt;
            stats.m_sampling_time_s += dt;
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
    fmt::print("Average velocity: ({:.6e}, {:.6e}, {:.6e}) m/s\n",
               avg_velocity.x(),
               avg_velocity.y(),
               avg_velocity.z());

    std::size_t total_acoustic_events               = 0;
    std::size_t total_intervalley_absorption_events = 0;
    std::size_t total_intervalley_emission_events   = 0;
    std::size_t total_impact_ionization_events      = 0;
    std::size_t total_self_scattering_events        = 0;

    for (const auto& p : m_particles) {
        const auto& events = p.history().scattering_events();
        total_acoustic_events += events[static_cast<std::size_t>(scattering_event::acoustic)];
        total_intervalley_absorption_events +=
            events[static_cast<std::size_t>(scattering_event::intervalley_absorption)];
        total_intervalley_emission_events += events[static_cast<std::size_t>(scattering_event::intervalley_emission)];
        total_impact_ionization_events += events[static_cast<std::size_t>(scattering_event::impact_ionization)];
        total_self_scattering_events += events[static_cast<std::size_t>(scattering_event::self_scattering)];
    }

    fmt::print("\nTotal acoustic events: {}\n", total_acoustic_events);
    fmt::print("Total intervalley absorption events: {}\n", total_intervalley_absorption_events);
    fmt::print("Total intervalley emission events: {}\n", total_intervalley_emission_events);
    fmt::print("Total impact ionization events: {}\n", total_impact_ionization_events);
    fmt::print("Total self-scattering events: {}\n\n", total_self_scattering_events);
    const auto& ii_stats = m_impact_ionization_coefficient_statistics;
    fmt::print("Impact ionization coefficient statistics:\n");
    fmt::print("  sampled II events: {}\n", ii_stats.m_events);
    fmt::print("  carrier-time: {:.6e} particle.s\n", ii_stats.m_carrier_time_s);
    fmt::print("  event rate per carrier: {:.6e} s^-1\n", ii_stats.event_rate_per_carrier_s_1());
    fmt::print("  drift velocity along field: {:.6e} m/s\n", ii_stats.average_drift_velocity_m_per_s());
    fmt::print("  ionization coefficient: {:.6e} cm^-1\n\n", ii_stats.ionization_coefficient_cm_1());

    m_observables.electric_field_V_per_m = m_cfg.m_electric_field.norm();

    if (m_observables.accumulated_time_s == 0.0) {
        throw std::runtime_error("no steady-state samples accumulated");
    }

    const double avg_vx_m_per_s = m_observables.weighted_velocity_x_m2_per_s2 / m_observables.accumulated_time_s;
    const double avg_energy_eV  = m_observables.weighted_kinetic_energy_eV_s / m_observables.accumulated_time_s;
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
    m_observables                              = {};
    m_impact_ionization_coefficient_statistics = {};
    const double warmup_time                   = m_cfg.m_warmup_fraction * m_cfg.m_final_time;
    std::size_t  number_threads                = 1;
    number_threads                             = static_cast<std::size_t>(std::max(number_threads, m_cfg.m_nb_threads));
    std::vector<amc_transport_kernel> thread_transports;
    thread_transports.reserve(number_threads);

    for (std::size_t i = 0; i < number_threads; ++i) {
        auto                 transport_cfg = make_transport_config(m_cfg);
        amc_transport_kernel transport(transport_cfg, static_cast<std::uint64_t>(1234 + 7919 * i));
        transport.initialize();
        thread_transports.push_back(std::move(transport));
    }

    double max_observed_total_rate = 0.0;

    double reduced_weighted_velocity_x = 0.0;
    double reduced_weighted_energy     = 0.0;
    double reduced_accumulated_time    = 0.0;

    std::size_t reduced_ii_events                     = 0;
    double      reduced_ii_carrier_time               = 0.0;
    double      reduced_ii_velocity_time_integral     = 0.0;
    double      reduced_ii_sampling_time              = 0.0;
    double      reduced_ii_ionization_coefficient_raw = 0.0;

    std::size_t reduced_impurity_events       = 0;
    double      reduced_impurity_carrier_time = 0.0;

#pragma omp parallel for if (m_cfg.m_nb_threads > 1) num_threads(m_cfg.m_nb_threads)          \
    reduction(max : max_observed_total_rate) reduction(+ : reduced_weighted_velocity_x,       \
                                                           reduced_weighted_energy,           \
                                                           reduced_accumulated_time,          \
                                                           reduced_ii_events,                 \
                                                           reduced_ii_carrier_time,           \
                                                           reduced_ii_velocity_time_integral, \
                                                           reduced_ii_sampling_time,          \
                                                           reduced_impurity_events,           \
                                                           reduced_impurity_carrier_time)
    for (std::int64_t idx_particle = 0; idx_particle < static_cast<std::int64_t>(m_particles.size()); ++idx_particle) {
        const int thread_id = omp_get_thread_num();
        auto&     transport = thread_transports[static_cast<std::size_t>(thread_id)];
        auto&     p         = m_particles[static_cast<std::size_t>(idx_particle)];

        std::size_t previous_particle_impact_ionization_events =
            particle_scattering_event_count(p, scattering_event::impact_ionization);

        std::size_t previous_particle_impurity_events = particle_scattering_event_count(p, scattering_event::impurity);

        while (p.state().time < m_cfg.m_final_time) {
            const double time_before_drift = p.state().time;

            const double tau            = transport.sample_free_flight_time();
            const double remaining_time = m_cfg.m_final_time - p.state().time;
            const double drift_time     = std::min(tau, remaining_time);
            transport.drift_particle(p, m_cfg.m_electric_field, drift_time);
            const double sampled_drift_time =
                sampled_time_after_warmup(time_before_drift, p.state().time, warmup_time, m_cfg.m_final_time);

            if (sampled_drift_time > 0.0) {
                reduced_weighted_velocity_x += p.state().velocity.x() * sampled_drift_time;
                reduced_weighted_energy += p.state().kinetic_energy * sampled_drift_time;
                reduced_accumulated_time += sampled_drift_time;
                reduced_ii_carrier_time += sampled_drift_time;
                reduced_ii_velocity_time_integral +=
                    particle_drift_velocity_along_field_m_per_s(p, m_cfg.m_electric_field) * sampled_drift_time;
                reduced_ii_sampling_time += sampled_drift_time;
                reduced_impurity_carrier_time += sampled_drift_time;
            }
            if (m_cfg.m_record_history) {
                p.record_state();
            }
            if (drift_time < tau) {
                break;
            }

            const double total_rate = transport.total_scattering_rate(p);
            max_observed_total_rate = std::max(max_observed_total_rate, total_rate);
            // Do not call ensure_gamma_max_covers() in parallel on shared m_transport.
            // This local transport is thread-private, so this is safe.
            transport.ensure_gamma_max_covers(total_rate);
            const double u = transport.uniform01();
            if (u < total_rate / transport.gamma_max()) {
                const auto channel = transport.select_scattering_channel(p);
                transport.apply_scattering_channel(p, channel);
            } else {
                p.increment_scattering_event_count();
                p.add_scattering_event(scattering_event::self_scattering);
            }
            const std::size_t current_particle_impact_ionization_events =
                particle_scattering_event_count(p, scattering_event::impact_ionization);
            const std::size_t new_particle_impact_ionization_events =
                current_particle_impact_ionization_events - previous_particle_impact_ionization_events;
            previous_particle_impact_ionization_events = current_particle_impact_ionization_events;
            const std::size_t current_particle_impurity_events =
                particle_scattering_event_count(p, scattering_event::impurity);
            const std::size_t new_particle_impurity_events =
                current_particle_impurity_events - previous_particle_impurity_events;
            previous_particle_impurity_events = current_particle_impurity_events;
            if (p.state().time >= warmup_time) {
                reduced_ii_events += new_particle_impact_ionization_events;
                reduced_impurity_events += new_particle_impurity_events;
            }
            if (m_cfg.m_record_history) {
                p.record_state();
            }
        }
        double particle_ii_coefficient_raw = p.compute_raw_impact_ionization_coefficient();
        fmt::print("Particle {} raw II coefficient: {:.6e} cm^-1\n", p.index(), particle_ii_coefficient_raw);
        reduced_ii_ionization_coefficient_raw += particle_ii_coefficient_raw;
    }

    m_observables.weighted_velocity_x_m2_per_s2 += reduced_weighted_velocity_x;
    m_observables.weighted_kinetic_energy_eV_s += reduced_weighted_energy;
    m_observables.accumulated_time_s += reduced_accumulated_time;

    m_impact_ionization_coefficient_statistics.m_events += reduced_ii_events;
    m_impact_ionization_coefficient_statistics.m_carrier_time_s += reduced_ii_carrier_time;
    m_impact_ionization_coefficient_statistics.m_drift_velocity_time_integral_m_per_s_times_s +=
        reduced_ii_velocity_time_integral;
    m_impact_ionization_coefficient_statistics.m_sampling_time_s += reduced_ii_sampling_time;
    m_impact_ionization_coefficient_statistics.m_raw_ii_coefficient_cm_1 += reduced_ii_ionization_coefficient_raw;

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
    fmt::print("Average velocity: ({:.6e}, {:.6e}, {:.6e}) m/s\n",
               avg_velocity.x(),
               avg_velocity.y(),
               avg_velocity.z());

    std::size_t total_acoustic_events               = 0;
    std::size_t total_intervalley_absorption_events = 0;
    std::size_t total_intervalley_emission_events   = 0;
    std::size_t total_impurity_events               = 0;
    std::size_t total_impact_ionization_events      = 0;
    std::size_t total_self_scattering_events        = 0;

    for (const auto& p : m_particles) {
        const auto& events = p.history().scattering_events();
        total_acoustic_events += events[static_cast<std::size_t>(scattering_event::acoustic)];
        total_intervalley_absorption_events +=
            events[static_cast<std::size_t>(scattering_event::intervalley_absorption)];
        total_intervalley_emission_events += events[static_cast<std::size_t>(scattering_event::intervalley_emission)];
        total_impact_ionization_events += events[static_cast<std::size_t>(scattering_event::impact_ionization)];
        total_impurity_events += events[static_cast<std::size_t>(scattering_event::impurity)];
        total_self_scattering_events += events[static_cast<std::size_t>(scattering_event::self_scattering)];
    }

    fmt::print("\nTotal acoustic events: {}\n", total_acoustic_events);
    fmt::print("Total intervalley absorption events: {}\n", total_intervalley_absorption_events);
    fmt::print("Total intervalley emission events: {}\n", total_intervalley_emission_events);
    fmt::print("Total impact ionization events: {}\n", total_impact_ionization_events);
    fmt::print("Total impurity events: {}\n", total_impurity_events);
    fmt::print("Total self-scattering events: {}\n\n", total_self_scattering_events);

    const auto& ii_stats = m_impact_ionization_coefficient_statistics;

    fmt::print("Impact ionization coefficient statistics:\n");
    fmt::print("  sampled II events: {}\n", ii_stats.m_events);
    fmt::print("  carrier-time: {:.6e} particle.s\n", ii_stats.m_carrier_time_s);
    fmt::print("  event rate per carrier: {:.6e} s^-1\n", ii_stats.event_rate_per_carrier_s_1());
    fmt::print("  drift velocity along field: {:.6e} m/s\n", ii_stats.average_drift_velocity_m_per_s());
    fmt::print("  ionization coefficient: {:.6e} cm^-1\n\n", ii_stats.ionization_coefficient_cm_1());
    fmt::print("  raw ionization coefficient (from particle data): {:.6e} cm^-1\n",
               ii_stats.m_raw_ii_coefficient_cm_1 / static_cast<double>(m_particles.size()));

    const double impurity_rate_per_carrier_s_1 =
        reduced_impurity_carrier_time > 0.0
            ? static_cast<double>(reduced_impurity_events) / reduced_impurity_carrier_time
            : 0.0;

    fmt::print("Impurity scattering statistics:\n");
    fmt::print("  sampled impurity events: {}\n", reduced_impurity_events);
    fmt::print("  carrier-time: {:.6e} particle.s\n", reduced_impurity_carrier_time);
    fmt::print("  event rate per carrier: {:.6e} s^-1\n\n", impurity_rate_per_carrier_s_1);

    m_observables.electric_field_V_per_m = m_cfg.m_electric_field.norm();

    if (m_observables.accumulated_time_s == 0.0) {
        throw std::runtime_error("no steady-state samples accumulated");
    }

    const double avg_vx_m_per_s = m_observables.weighted_velocity_x_m2_per_s2 / m_observables.accumulated_time_s;
    const double avg_energy_eV  = m_observables.weighted_kinetic_energy_eV_s / m_observables.accumulated_time_s;
    fmt::print("Steady-state average vx: {:.6e} m/s\n", avg_vx_m_per_s);
    fmt::print("Steady-state average energy: {:.6f} eV\n", avg_energy_eV);
}

void bulk_amc_simulation::export_observables_to_csv(const std::string& filename) const {
    std::ofstream file(filename, std::ios::app);
    if (!file.is_open()) {
        fmt::print(stderr, "Failed to open file for writing: {}\n", filename);
        return;
    }

    if (file.tellp() == 0) {
        file << "charge_C,"
                "temperature_K,"
                "electric_field_V_per_m,"
                "impurity_density_cm_3,"
                "mean_velocity_x_m_per_s,"
                "mean_kinetic_energy_eV,"
                "sample_count,"
                "impact_ionization_events,"
                "impact_ionization_rate_per_carrier_s_1,"
                "impact_ionization_drift_velocity_m_per_s,"
                "impact_ionization_coefficient_cm_1\n";
    }

    if (m_observables.accumulated_time_s <= 0.0) {
        fmt::print(stderr, "No observables accumulated; skipping CSV export.\n");
        return;
    }

    const auto& ii_stats = m_impact_ionization_coefficient_statistics;

    file << fmt::format("{},{},{},{},{},{},{},{},{},{},{}\n",
                        signed_charge_C(m_cfg.m_carrier_type),
                        m_cfg.m_lattice_temperature,
                        m_observables.electric_field_V_per_m,
                        m_cfg.m_impurity_density_cm_3,
                        m_observables.weighted_velocity_x_m2_per_s2 / m_observables.accumulated_time_s,
                        m_observables.weighted_kinetic_energy_eV_s / m_observables.accumulated_time_s,
                        m_observables.accumulated_time_s,
                        ii_stats.m_events,
                        ii_stats.event_rate_per_carrier_s_1(),
                        ii_stats.average_drift_velocity_m_per_s(),
                        ii_stats.ionization_coefficient_cm_1());
}

void bulk_amc_simulation::export_particles_history_to_csv(const std::string& prefix_name) const {
    for (const auto& p : m_particles) {
        const auto&       history  = p.history();
        const std::string filename = fmt::format("{}_particle_{}.csv", prefix_name, p.index());
        std::ofstream     file(filename);

        if (!file.is_open()) {
            fmt::print(stderr, "Failed to open file for writing: {}\n", filename);
            continue;
        }

        file << "time,"
                "position_x,position_y,position_z,"
                "local_k_x,local_k_y,local_k_z,"
                "velocity_x,velocity_y,velocity_z,"
                "kinetic_energy,gamma,valley_index\n";

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