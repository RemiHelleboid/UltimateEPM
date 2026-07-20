/**
 * @author remzerrr (remi.helleboid@gmail.com)
 * @file admc_carrier.hpp
 * @brief
 * @version 0.1
 * @date 2026-06-15
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "vector.hpp"

namespace uepm::ADMC {

using vector3 = uepm::common::vector3;

enum class carrier_type : std::int8_t { electron, hole };

constexpr double carrier_charge_sign(carrier_type type) {
    switch (type) {
        case carrier_type::electron:
            return -1.0;
        case carrier_type::hole:
            return 1.0;
    }
    throw std::invalid_argument("invalid ADMC carrier type");
}

constexpr std::string_view carrier_type_name(carrier_type type) {
    switch (type) {
        case carrier_type::electron:
            return "electron";
        case carrier_type::hole:
            return "hole";
    }
    return "unknown";
}

struct admc_particle_state {
    double  time_s = 0.0;
    vector3 position_m{};
    vector3 previous_position_m{};
    vector3 electric_field_V_per_m{};
    vector3 drift_velocity_m_per_s{};
    vector3 total_velocity_m_per_s{};
    double  doping_concentration_cm_3 = 0.0;
    double  lattice_temperature_K     = 300.0;
    double  mobility_m2_per_V_s       = 0.0;
    double  diffusion_m2_per_s        = 0.0;
};

struct admc_particle_snapshot {
    double  time_s = 0.0;
    vector3 position_m{};
    vector3 electric_field_V_per_m{};
    vector3 drift_velocity_m_per_s{};
    vector3 total_velocity_m_per_s{};
    double  doping_concentration_cm_3 = 0.0;
    double  lattice_temperature_K     = 300.0;
    double  mobility_m2_per_V_s       = 0.0;
    double  diffusion_m2_per_s        = 0.0;
};

class admc_particle_history {
 public:
    admc_particle_history() = default;
    explicit admc_particle_history(std::size_t particle_index) : m_particle_index(particle_index) {}

    std::size_t                                particle_index() const noexcept { return m_particle_index; }
    std::size_t                                recorded_number_of_steps() const noexcept { return m_snapshots.size(); }
    const std::vector<admc_particle_snapshot>& snapshots() const noexcept { return m_snapshots; }

    void reserve(std::size_t number_of_snapshots) { m_snapshots.reserve(number_of_snapshots); }
    void clear() noexcept { m_snapshots.clear(); }
    void record(const admc_particle_state& state) {
        m_snapshots.push_back(admc_particle_snapshot{
            .time_s                    = state.time_s,
            .position_m                = state.position_m,
            .electric_field_V_per_m    = state.electric_field_V_per_m,
            .drift_velocity_m_per_s    = state.drift_velocity_m_per_s,
            .total_velocity_m_per_s    = state.total_velocity_m_per_s,
            .doping_concentration_cm_3 = state.doping_concentration_cm_3,
            .lattice_temperature_K     = state.lattice_temperature_K,
            .mobility_m2_per_V_s       = state.mobility_m2_per_V_s,
            .diffusion_m2_per_s        = state.diffusion_m2_per_s,
        });
    }

 private:
    std::size_t                         m_particle_index = 0;
    std::vector<admc_particle_snapshot> m_snapshots;
};

class admc_particle {
 public:
    admc_particle(std::size_t index, carrier_type type, vector3 initial_position_m = {})
        : m_index(index),
          m_type(type),
          m_history(index) {
        m_state.position_m          = initial_position_m;
        m_state.previous_position_m = initial_position_m;
    }

    std::size_t                  index() const noexcept { return m_index; }
    carrier_type                 type() const noexcept { return m_type; }
    const admc_particle_state&   state() const noexcept { return m_state; }
    admc_particle_state&         state() noexcept { return m_state; }
    const admc_particle_history& history() const noexcept { return m_history; }
    admc_particle_history&       history() noexcept { return m_history; }
    void                         record_state() { m_history.record(m_state); }

 private:
    std::size_t           m_index;
    carrier_type          m_type;
    admc_particle_state   m_state{};
    admc_particle_history m_history{};
};

}  // namespace uepm::ADMC
