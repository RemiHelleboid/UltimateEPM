/**
 * @file scattering_events.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-04-16
 * 
 * 
 */

#pragma once

#include <cstddef>

namespace uepm::amc {

enum class scattering_event : std::size_t { acoustic = 0, intervalley_absorption, intervalley_emission, self_scattering, count };

}  // namespace uepm::amc