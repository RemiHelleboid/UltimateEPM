/**
 * @file scattering_channels.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-04-16
 * 
 * 
 */

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string_view>

namespace uepm::amc {

struct intervalley_phonon_branch;

enum class scattering_mechanism : std::uint8_t { acoustic, intervalley, impurity, impact_ionization };

enum class intervalley_process : std::uint8_t { none, absorption, emission };

struct scattering_channel {
    scattering_mechanism             mechanism         = scattering_mechanism::acoustic;
    double                           rate_s_1          = 0.0;
    double                           final_energy_eV   = 0.0;
    std::size_t                      destination_index = 0;
    const intervalley_phonon_branch* branch            = nullptr;
    intervalley_process              process           = intervalley_process::none;
    std::string_view                 transition_name{};
};

class scattering_channel_list {
 public:
    static constexpr std::size_t capacity = 16;

    using const_iterator = std::array<scattering_channel, capacity>::const_iterator;

    void push_back(const scattering_channel& channel) {
        if (m_size == capacity) {
            throw std::length_error("too many AMC scattering channels");
        }
        m_channels[m_size++] = channel;
    }

    [[nodiscard]] bool empty() const noexcept { return m_size == 0; }
    [[nodiscard]] std::size_t size() const noexcept { return m_size; }
    [[nodiscard]] const scattering_channel& front() const { return m_channels.front(); }
    [[nodiscard]] const_iterator begin() const noexcept { return m_channels.begin(); }
    [[nodiscard]] const_iterator end() const noexcept {
        return m_channels.begin() + static_cast<std::ptrdiff_t>(m_size);
    }

 private:
    std::array<scattering_channel, capacity> m_channels{};
    std::size_t                              m_size = 0;
};
}  // namespace uepm::amc
