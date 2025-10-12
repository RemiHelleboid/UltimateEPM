/**
 * @file physical_constants.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief  Common physical constants used in the project.
 * @version 0.1
 * @date 2025-10-12
 *
 *
 */

#pragma once

#include <cmath>
#include <numbers>

namespace uepm {

namespace Constants {

constexpr double h   = 6.62607015e-34;   // J·s (exact)
constexpr double k_B = 1.380649e-23;     // J/K (exact)
constexpr double q_e = 1.602176634e-19;  // C = J/eV (exact)
constexpr double N_A = 6.02214076e23;    // 1/mol (exact)
constexpr double c   = 299792458.0;      // m/s (exact)

// === Derived fundamentals ===
constexpr double pi       = std::numbers::pi_v<double>;
constexpr double h_bar    = h / (2.0 * pi);   // J·s
constexpr double eV_to_J  = q_e;              // J/eV
constexpr double h_bar_eV = h_bar / eV_to_J;  // eV·s
constexpr double k_b_eV   = k_B / eV_to_J;    // eV/K

constexpr double eps_0 = 8.8541878128e-12;  // F/m
constexpr double m_e   = 9.1093837015e-31;  // kg

constexpr double bohr_radius    = 5.29177210903e-11;    // m
constexpr double Hartree_to_J   = 4.3597447222071e-18;  // J
constexpr double Hartree_to_eV  = Hartree_to_J / eV_to_J;
constexpr double Ryd_to_eV      = Hartree_to_eV / 2.0;  // ~13.605693122994
constexpr double Ryd_to_Hartree = 0.5;

// conversions
constexpr double angstrom_to_m   = 1e-10;              // m/Å
constexpr double eV_to_cm_inv    = 8065.54429;         // cm^-1/eV
constexpr double eV_to_J_per_mol = eV_to_J * N_A;      // J/mol per eV
constexpr double eV_to_kg        = eV_to_J / (c * c);  // kg/eV (~1.78266192e-36)

// Back-conversions
constexpr double eV_to_Hartree = 1.0 / Hartree_to_eV;
constexpr double eV_to_Ryd     = 1.0 / Ryd_to_eV;

}  // namespace Constants
}  // namespace uepm