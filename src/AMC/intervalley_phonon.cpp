/**
 * @file intervalley_phonon.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-04-13
 *
 *
 */

#include "intervalley_phonon.hpp"

namespace uepm::amc {

std::vector<intervalley_phonon_branch> make_silicon_intervalley_phonon_branches() {
    std::vector<intervalley_phonon_branch> branches;
    branches.reserve(5);

    // Values inspired by the silicon MC model summarized in Aubry-Fortuna et al.
    // Energies are in eV here.
    // D0 values are converted from eV/cm to eV/m.
    // The paper lists:
    // g1 = 11.4 meV, g2 = 18.8 meV, g3 = 63.2 meV,
    // f1 = 21.9 meV, f2 = 46.3 meV.
    // with g1, g2, f1 treated as first-order,
    // and g3, f2 treated as zeroth-order. :contentReference[oaicite:1]{index=1}

    branches.push_back(intervalley_phonon_branch{.m_name                    = "g1_TA",
                                                 .m_family                  = intervalley_family::g,
                                                 .m_order                   = intervalley_order::first,
                                                 .m_phonon_energy_eV        = 11.4e-3,
                                                 .m_deformation_potential_0 = 0.0,
                                                 .m_deformation_potential_1 = 3.0,
                                                 .m_final_valley_count      = 1});

    branches.push_back(intervalley_phonon_branch{.m_name                    = "g2_LA",
                                                 .m_family                  = intervalley_family::g,
                                                 .m_order                   = intervalley_order::first,
                                                 .m_phonon_energy_eV        = 18.8e-3,
                                                 .m_deformation_potential_0 = 0.0,
                                                 .m_deformation_potential_1 = 3.0,
                                                 .m_final_valley_count      = 1});

    branches.push_back(intervalley_phonon_branch{.m_name                    = "g3_LO",
                                                 .m_family                  = intervalley_family::g,
                                                 .m_order                   = intervalley_order::zeroth,
                                                 .m_phonon_energy_eV        = 63.2e-3,
                                                 .m_deformation_potential_0 = 3.4e10,
                                                 .m_deformation_potential_1 = 0.0,
                                                 .m_final_valley_count      = 1});

    branches.push_back(intervalley_phonon_branch{.m_name                    = "f1_TA",
                                                 .m_family                  = intervalley_family::f,
                                                 .m_order                   = intervalley_order::first,
                                                 .m_phonon_energy_eV        = 21.9e-3,
                                                 .m_deformation_potential_0 = 0.0,
                                                 .m_deformation_potential_1 = 3.0,
                                                 .m_final_valley_count      = 4});

    branches.push_back(intervalley_phonon_branch{.m_name                    = "f2_LA",
                                                 .m_family                  = intervalley_family::f,
                                                 .m_order                   = intervalley_order::zeroth,
                                                 .m_phonon_energy_eV        = 46.3e-3,
                                                 .m_deformation_potential_0 = 3.4e10,
                                                 .m_deformation_potential_1 = 0.0,
                                                 .m_final_valley_count      = 4});

    return branches;
}

}  // namespace uepm::amc