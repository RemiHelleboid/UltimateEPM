/**
 * @file scattering_event.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-04-01
 *
 *
 */

#pragma once

#include "particle_amc.hpp"
#include "valley_model.hpp"

namespace uepm::amc {

enum class scattering_type : std::size_t { acoustic, intervalley_g, intervalley_f, self_scattering };

enum class scattering_event : std::size_t { acoustic = 0, intervalley_absorption, intervalley_emission, self_scattering, count };
}  // namespace uepm::amc