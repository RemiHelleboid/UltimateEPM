/**
 * @file amc_material_model.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-26
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <cstddef>
#include <vector>

#include "valley_model.hpp"

namespace uepm::amc {

struct impact_ionization_parameters {
    double m_threshold_eV  = 1.1;
    double m_prefactor_s_1 = 1.0e11;
    double m_exponent      = 4.6;
};

struct carrier_impact_ionization_parameters {
    impact_ionization_parameters m_electron;
    impact_ionization_parameters m_hole;
};
struct hole_optical_transition {
    const char* name                           = "";
    std::size_t initial_band                   = 0;
    std::size_t final_band                     = 0;
    double      phonon_energy_eV               = 0.0;
    double      deformation_potential_eV_per_m = 0.0;
    double      overlap_factor                 = 1.0;
};

std::vector<valley_model>            make_silicon_delta_valleys();
std::vector<valley_model>            make_silicon_hole_bands();
std::vector<hole_optical_transition> make_silicon_hole_optical_transitions();
carrier_impact_ionization_parameters make_silicon_impact_ionization_parameters();

}  // namespace uepm::amc
