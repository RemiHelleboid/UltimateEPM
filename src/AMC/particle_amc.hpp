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
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "element.hpp"
#include "physical_constants.hpp"
#include "scattering_events.hpp"
#include "vector.hpp"

namespace uepm::amc {

enum class particle_type : std::int8_t { electron, hole };

std::string_view carrier_type_to_string(particle_type type);
double           signed_charge_C(particle_type type);

using vector3 = uepm::mesh::vector3;
using element = uepm::mesh::element;

struct particle_state {
    double      time                  = 0.0;
    double      lattice_temperature_K = 300.0;  // in K
    vector3     position{};
    vector3     local_k{};
    vector3     velocity{};
    vector3     electric_field{};                   // in V/m
    double      doping_concentration_cm_3   = 0.0;  // signed or net doping
    double      impurity_concentration_cm_3 = 0.0;  // positive scattering-center density
    double      kinetic_energy              = 0.0;  // in eV
    double      gamma        = 0.0;  // gamma = E for parabolic, gamma = E * (1 + alpha * E) for Kane (in eV)
    std::size_t valley_index = 0;

    /**
     * @brief Pointer to the element containing the particle. This is used to avoid having to search for the containing
     * element at each step, which can be costly. It is updated after each scattering event and after each free flight.
     *
     */
    element* m_containing_element = nullptr;

    /**
     * @brief Flag to indicate if the particle has crossed a contact during its free flight.
     *
     */
    bool m_crossed_contact = false;

    vector3 previous_position{};
};

struct particle_snapshot {
    double      time = 0.0;
    vector3     position{};
    vector3     local_k{};
    vector3     velocity{};
    double      kinetic_energy      = 0.0;
    double      gamma               = 0.0;
    std::size_t valley_index        = 0;
    double      electric_field_norm = 0.0;
};

class particle_history {
 public:
    static constexpr std::size_t num_scattering_channels = static_cast<std::size_t>(scattering_event::count);

 private:
    std::size_t                                      m_particle_index          = 0;
    std::size_t                                      m_total_scattering_events = 0;
    std::vector<particle_snapshot>                   m_snapshots;
    std::array<std::size_t, num_scattering_channels> m_scattering_events{};
    std::unordered_map<std::string, std::size_t>     m_transition_events;

 public:
    particle_history() = default;

    explicit particle_history(std::size_t particle_index) : m_particle_index(particle_index) {}
    std::size_t particle_index() const noexcept { return m_particle_index; }
    std::size_t total_scattering_event_count() const noexcept { return m_total_scattering_events; }
    std::size_t recorded_number_of_steps() const noexcept { return m_snapshots.size(); }
    const std::vector<particle_snapshot>&                   snapshots() const noexcept { return m_snapshots; }
    const std::array<std::size_t, num_scattering_channels>& scattering_events() const noexcept {
        return m_scattering_events;
    }
    const std::unordered_map<std::string, std::size_t>& transition_events() const noexcept {
        return m_transition_events;
    }
    void reserve(std::size_t n_steps) { m_snapshots.reserve(n_steps); }
    void clear() noexcept {
        m_total_scattering_events = 0;
        m_snapshots.clear();
        m_scattering_events.fill(0);
        m_transition_events.clear();
    }

    void set_particle_index(std::size_t particle_index) noexcept { m_particle_index = particle_index; }
    void increment_scattering_event_count() noexcept { ++m_total_scattering_events; }
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
    void add_transition_event(std::string_view transition_name) {
        if (!transition_name.empty()) {
            ++m_transition_events[std::string(transition_name)];
        }
    }
    void export_trajectory_as_csv(const std::string& filename) const;
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

    particle_amc(std::size_t index, particle_type type, const particle_state& initial_state, double weight = 1.0)
        : m_index(index),
          m_type(type),
          m_state(initial_state),
          m_weight(weight),
          m_history(index) {
        if (weight <= 0.0) {
            throw std::invalid_argument("particle weight must be > 0");
        }
    }

    std::size_t   index() const noexcept { return m_index; }
    particle_type type() const noexcept { return m_type; }
    double        weight() const noexcept { return m_weight; }
    double        get_signed_charge() const noexcept {
        constexpr double q = 1.602176634e-19;
        return (m_type == particle_type::electron) ? -q : q;
    }
    const particle_state&   state() const noexcept { return m_state; }
    particle_state&         state() noexcept { return m_state; }
    const particle_history& history() const noexcept { return m_history; }
    particle_history&       history() noexcept { return m_history; }
    void                    set_index(std::size_t index) noexcept {
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
    void advance_time(double dt) { m_state.time += dt; }

    void set_position(const vector3& position) noexcept { m_state.position = position; }
    void translate(const vector3& dr) noexcept { m_state.position += dr; }
    void set_local_k(const vector3& local_k) noexcept { m_state.local_k = local_k; }
    void set_velocity(const vector3& velocity) noexcept { m_state.velocity = velocity; }
    void set_kinetic_energy(double energy) { m_state.kinetic_energy = energy; }
    void set_gamma(double gamma) { m_state.gamma = gamma; }

    double get_lattice_temperature() const noexcept { return m_state.lattice_temperature_K; }
    void   set_lattice_temperature(double temperature_K) noexcept { m_state.lattice_temperature_K = temperature_K; }

    void     set_crossed_contact(bool crossed) noexcept { m_state.m_crossed_contact = crossed; }
    void     set_valley_index(std::size_t valley_index) noexcept { m_state.valley_index = valley_index; }
    void     set_containing_element(mesh::element* element) noexcept { m_state.m_containing_element = element; }
    element* get_containing_element() const noexcept { return m_state.m_containing_element; }
    void     increment_scattering_event_count() noexcept { m_history.increment_scattering_event_count(); }
    void     record_state() { m_history.record(m_state); }
    void     add_scattering_event(scattering_event event) noexcept { m_history.add_event(event); }
    void     add_transition_event(std::string_view transition_name) { m_history.add_transition_event(transition_name); }
    void     reset_history() noexcept { m_history.clear(); }
    void     set_data_from_device(int m_dimension);
    void     print_info() const;
    double   compute_raw_impact_ionization_coefficient() const;
    void export_trajectory_as_csv(const std::string& filename) const { m_history.export_trajectory_as_csv(filename); }
};

}  // namespace uepm::amc
