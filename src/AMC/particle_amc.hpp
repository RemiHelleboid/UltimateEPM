/**
 * @file particle_amc.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-03-30
 * 
 * 
 */

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include "vector.hpp"
#include "physical_constants.hpp"
#include "scattering_events.hpp"

namespace uepm::amc {

enum class particle_type : std::int8_t { electron, hole };


using vector3 = uepm::mesh_bz::vector3;

struct particle_state {
    double      time = 0.0;
    vector3     position{};
    vector3     local_k{};
    vector3     velocity{};
    double      kinetic_energy = 0.0; // in eV
    double      gamma          = 0.0; // gamma = E for parabolic, gamma = E * (1 + alpha * E) for Kane (in eV)
    std::size_t valley_index   = 0;
};

struct particle_snapshot {
    double      time = 0.0;
    vector3     position{};
    vector3     local_k{};
    vector3     velocity{};
    double      kinetic_energy = 0.0;
    double      gamma          = 0.0;
    std::size_t valley_index   = 0;
};

class particle_history {
 public:
    static constexpr std::size_t num_scattering_channels = static_cast<std::size_t>(scattering_event::count);

 private:
    std::size_t                                      m_particle_index = 0;
    std::size_t                                      m_total_nb_steps = 0;
    std::vector<particle_snapshot>                   m_snapshots;
    std::array<std::size_t, num_scattering_channels> m_scattering_events{};

 public:
    particle_history() = default;

    explicit particle_history(std::size_t particle_index) : m_particle_index(particle_index) {}

    std::size_t particle_index() const noexcept { return m_particle_index; }

    std::size_t total_number_of_steps() const noexcept { return m_total_nb_steps; }

    std::size_t recorded_number_of_steps() const noexcept { return m_snapshots.size(); }

    const std::vector<particle_snapshot>& snapshots() const noexcept { return m_snapshots; }

    const std::array<std::size_t, num_scattering_channels>& scattering_events() const noexcept { return m_scattering_events; }

    void reserve(std::size_t n_steps) { m_snapshots.reserve(n_steps); }

    void clear() noexcept {
        m_total_nb_steps = 0;
        m_snapshots.clear();
        m_scattering_events.fill(0);
    }

    void set_particle_index(std::size_t particle_index) noexcept { m_particle_index = particle_index; }

    void increment_total_steps() noexcept { ++m_total_nb_steps; }

    void record(const particle_state& state) {
        m_snapshots.push_back(particle_snapshot{.time           = state.time,
                                                .position       = state.position,
                                                .local_k        = state.local_k,
                                                .velocity       = state.velocity,
                                                .kinetic_energy = state.kinetic_energy,
                                                .gamma          = state.gamma,
                                                .valley_index   = state.valley_index});
    }

    void add_event(scattering_event event) noexcept {
        const auto index = static_cast<std::size_t>(event);
        if (index < m_scattering_events.size()) {
            ++m_scattering_events[index];
        }
    }
};

class particle_amc {
 private:
    std::size_t      m_index  = 0;
    particle_type    m_type   = particle_type::electron;
    double           m_weight = 1.0;
    particle_state   m_state{};
    particle_history m_history{};

 public:
    particle_amc() = default;

    explicit particle_amc(std::size_t index, particle_type type = particle_type::electron, double weight = 1.0)
        : m_index(index),
          m_type(type),
          m_weight(weight),
          m_history(index) {
        if (weight <= 0.0) {
            throw std::invalid_argument("particle weight must be > 0");
        }
    }

    std::size_t index() const noexcept { return m_index; }

    particle_type type() const noexcept { return m_type; }

    double weight() const noexcept { return m_weight; }

    double signed_charge() const noexcept {
        constexpr double q = 1.602176634e-19;
        return (m_type == particle_type::electron) ? -q : q;
    }

    const particle_state& state() const noexcept { return m_state; }

    particle_state& state() noexcept { return m_state; }

    const particle_history& history() const noexcept { return m_history; }

    particle_history& history() noexcept { return m_history; }

    void set_index(std::size_t index) noexcept {
        m_index = index;
        m_history.set_particle_index(index);
    }

    void set_weight(double weight) {
        if (weight <= 0.0) {
            throw std::invalid_argument("particle weight must be > 0");
        }
        m_weight = weight;
    }

    void set_type(particle_type type) noexcept { m_type = type; }

    void set_state(const particle_state& state) noexcept { m_state = state; }

    void advance_time(double dt) {
        if (dt < 0.0) {
            throw std::invalid_argument("time step must be >= 0");
        }
        m_state.time += dt;
    }

    void translate(const vector3& dr) noexcept { m_state.position += dr; }

    void set_local_k(const vector3& local_k) noexcept { m_state.local_k = local_k; }

    void set_velocity(const vector3& velocity) noexcept { m_state.velocity = velocity; }

    void set_kinetic_energy(double energy) {
        if (energy < 0.0) {
            throw std::invalid_argument("kinetic energy must be >= 0");
        }
        m_state.kinetic_energy = energy;
    }

    void set_gamma(double gamma) {
        if (gamma < 0.0) {
            throw std::invalid_argument("gamma must be >= 0");
        }
        m_state.gamma = gamma;
    }

    void set_valley_index(std::size_t valley_index) noexcept { m_state.valley_index = valley_index; }

    void increment_total_step() noexcept { m_history.increment_total_steps(); }

    void record_state() { m_history.record(m_state); }

    void add_scattering_event(scattering_event event) noexcept { m_history.add_event(event); }

    void reset_history() noexcept { m_history.clear(); }

    void print_info() const;
};

}  // namespace uepm::amc