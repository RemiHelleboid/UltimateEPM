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
#include <string>
#include <vector>

#include "intervalley_phonon.hpp"
#include "valley_model.hpp"

namespace uepm::amc {

struct dielectric_properties {
    double epsilon_r = 1.0;
};

struct acoustic_scattering_parameters {
    double mass_density_kg_per_m3       = 0.0;
    double sound_velocity_m_per_s       = 0.0;
    double deformation_potential_eV     = 0.0;
    double overlap_factor               = 1.0;
};

struct impact_ionization_parameters {
    double m_threshold_eV  = 1.1;
    double m_prefactor_s_1 = 1.0e11;
    double m_exponent      = 4.6;
};

struct impurity_mobility_parameters {
    double m_mu0_cm2_per_V_s    = 1400.0;
    double m_mu_min_cm2_per_V_s = 50.0;
    double m_n_ref_cm_3         = 1.0e17;
    double m_alpha              = 0.7;
};

struct carrier_impurity_mobility_parameters {
    impurity_mobility_parameters m_electron;
    impurity_mobility_parameters m_hole;
};

struct carrier_impact_ionization_parameters {
    impact_ionization_parameters m_electron;
    impact_ionization_parameters m_hole;
};

struct hole_optical_transition {
    std::string name;
    std::size_t initial_band                   = 0;
    std::size_t final_band                     = 0;
    double      phonon_energy_eV               = 0.0;
    double      deformation_potential_eV_per_m = 0.0;
    double      overlap_factor                 = 1.0;
};

struct amc_material_model {
    std::string m_name;
    std::string m_symbol;

    dielectric_properties          m_dielectric;
    acoustic_scattering_parameters m_electron_acoustic;
    acoustic_scattering_parameters m_hole_acoustic;

    std::vector<valley_model>              m_electron_valleys;
    std::vector<intervalley_phonon_branch> m_electron_intervalley_transitions;
    std::vector<valley_model>              m_hole_bands;
    std::vector<hole_optical_transition>   m_hole_optical_transitions;

    carrier_impurity_mobility_parameters m_impurity_mobility;
    carrier_impact_ionization_parameters m_impact_ionization;

    void validate() const;
};

dielectric_properties                make_silicon_dielectric_properties();
carrier_impurity_mobility_parameters make_silicon_impurity_mobility_parameters();
std::vector<valley_model>            make_silicon_delta_valleys();
std::vector<valley_model>            make_silicon_hole_bands();
std::vector<hole_optical_transition> make_silicon_hole_optical_transitions();
carrier_impact_ionization_parameters make_silicon_impact_ionization_parameters();
amc_material_model                   make_silicon_amc_material_model();

}  // namespace uepm::amc
