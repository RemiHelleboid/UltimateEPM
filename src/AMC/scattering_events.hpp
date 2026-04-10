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

struct scattering_result {
    scattering_type type             = scattering_type::self_scattering;
    std::size_t     new_valley_index = 0;
    double          new_energy       = 0.0;
    vector3         new_local_k{};
    bool            real_scattering = false;
};

} // namespace uepm::amc