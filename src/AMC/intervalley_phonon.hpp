/**
 * @file intervalley_phonon.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-04-13
 * 
 * 
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace uepm::amc {

enum class intervalley_family : std::uint8_t { f, g };

enum class intervalley_order : std::uint8_t { zeroth, first };

struct intervalley_phonon_branch {
    std::string        m_name;
    intervalley_family m_family = intervalley_family::f;
    intervalley_order  m_order  = intervalley_order::zeroth;

    double m_phonon_energy_eV = 0.0;

    // Use only one depending on the order:
    // zeroth-order -> m_deformation_potential_0
    // first-order  -> m_deformation_potential_1
    double m_deformation_potential_0 = 0.0;  // eV / m
    double m_deformation_potential_1 = 0.0;  // eV

    std::size_t m_final_valley_count = 0;

    bool is_zeroth_order() const noexcept { return m_order == intervalley_order::zeroth; }

    bool is_first_order() const noexcept { return m_order == intervalley_order::first; }
};

std::vector<intervalley_phonon_branch> make_silicon_intervalley_phonon_branches();

}  // namespace uepm::amc