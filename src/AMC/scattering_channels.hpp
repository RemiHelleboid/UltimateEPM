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

#include <cstddef>
#include <cstdint>

namespace uepm::amc {

struct intervalley_phonon_branch;

enum class scattering_mechanism : std::uint8_t { acoustic, intervalley };

enum class intervalley_process : std::uint8_t { none, absorption, emission };

struct scattering_channel {
    scattering_mechanism             mechanism         = scattering_mechanism::acoustic;
    double                           rate_s_1          = 0.0;
    double                           final_energy_eV   = 0.0;
    std::size_t                      destination_index = 0;
    const intervalley_phonon_branch* branch            = nullptr;
    intervalley_process              process           = intervalley_process::none;
};
}  // namespace uepm::amc