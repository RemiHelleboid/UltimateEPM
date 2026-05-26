/**
 * @file amc_scattering_model.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-05-26
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#include "amc_scattering_model.hpp"

namespace uepm::amc {

double bose_einstein_occupation(double phonon_energy_eV, double temperature_K) {
    if (temperature_K <= 0.0) {
        return 0.0;
    }

    const double x = phonon_energy_eV / (uepm::constants::k_b_eV * temperature_K);
    return 1.0 / (std::exp(x) - 1.0);
}

double acoustic_scattering_rate_silicon(const valley_model& valley, double energy_eV, double temperature_K) {
    if (temperature_K <= 0.0) {
        return 0.0;
    }

    constexpr double rho_si  = 2.329e3;  // kg / m^3
    constexpr double u_l     = 9.0e3;    // m / s
    constexpr double u_t     = 5.4e3;    // m / s
    constexpr double u_avg   = (u_l + 2.0 * u_t) / 3.0;
    constexpr double D_ac_eV = 6.6;  // eV

    const double energy_clamped_eV = std::max(energy_eV, 1.0e-9);
    const double energy_J          = energy_clamped_eV * uepm::constants::eV_to_J;

    const double alpha_per_J = valley.non_parabolicity() / uepm::constants::eV_to_J;

    const double mt = valley.transverse_effective_mass();
    const double ml = valley.longitudinal_effective_mass();
    const double mD = std::cbrt(mt * mt * ml);

    const double D_ac_J = D_ac_eV * uepm::constants::eV_to_J;

    const double gamma_factor = std::sqrt(energy_J * (1.0 + alpha_per_J * energy_J)) * (1.0 + 2.0 * alpha_per_J * energy_J);

    const double prefactor = std::sqrt(2.0) * uepm::constants::k_B * temperature_K * std::pow(mD, 1.5) * (D_ac_J * D_ac_J) /
                             (uepm::constants::pi * std::pow(uepm::constants::h_bar, 4) * rho_si * u_avg * u_avg);

    return prefactor * gamma_factor;
}

double acoustic_scattering_rate_silicon_holes(const valley_model& band, double energy_eV, double temperature_K) {
    if (temperature_K <= 0.0) {
        return 0.0;
    }

    constexpr double rho_si         = 2.329e3;  // kg / m^3
    constexpr double v_s            = 6.6e3;    // m / s
    constexpr double D_ac_eV        = 5.5;      // eV
    constexpr double overlap_factor = 0.5;      // A + B/3 for intra-band holes

    const double energy_clamped_eV = std::max(energy_eV, 1.0e-9);
    const double energy_J          = energy_clamped_eV * uepm::constants::eV_to_J;

    const double alpha_per_J = band.non_parabolicity() / uepm::constants::eV_to_J;

    const double mt = band.transverse_effective_mass();
    const double ml = band.longitudinal_effective_mass();
    const double mD = std::cbrt(mt * mt * ml);

    const double D_ac_J = D_ac_eV * uepm::constants::eV_to_J;

    const double gamma_factor = std::sqrt(energy_J * (1.0 + alpha_per_J * energy_J)) * (1.0 + 2.0 * alpha_per_J * energy_J);

    const double prefactor = std::sqrt(2.0) * uepm::constants::k_B * temperature_K * std::pow(mD, 1.5) * (D_ac_J * D_ac_J) *
                             overlap_factor / (uepm::constants::pi * std::pow(uepm::constants::h_bar, 4) * rho_si * v_s * v_s);

    return prefactor * gamma_factor;
}

double silicon_density_of_states_mass(const valley_model& valley) {
    const double mt = valley.transverse_effective_mass();
    const double ml = valley.longitudinal_effective_mass();
    return std::cbrt(mt * mt * ml);
}

double silicon_nonparabolicity_per_joule(const valley_model& valley) { return valley.non_parabolicity() / uepm::constants::eV_to_J; }

double intervalley_zeroth_order_rate(const valley_model&              valley,
                                     const intervalley_phonon_branch& branch,
                                     double                           initial_energy_eV,
                                     bool                             absorption,
                                     double                           temperature_K) {
    constexpr double rho_si = 2.329e3;  // kg / m^3

    const double phonon_energy_eV = branch.m_phonon_energy_eV;
    const double final_energy_eV  = absorption ? (initial_energy_eV + phonon_energy_eV) : (initial_energy_eV - phonon_energy_eV);

    if (final_energy_eV < 0.0) {
        return 0.0;
    }

    const double Nop           = bose_einstein_occupation(phonon_energy_eV, temperature_K);
    const double phonon_factor = absorption ? Nop : (Nop + 1.0);

    const double alpha_per_J = silicon_nonparabolicity_per_joule(valley);
    const double mD          = silicon_density_of_states_mass(valley);

    const double final_energy_J  = final_energy_eV * uepm::constants::eV_to_J;
    const double phonon_energy_J = phonon_energy_eV * uepm::constants::eV_to_J;

    const double D0_J_per_m = branch.m_deformation_potential_0 * uepm::constants::eV_to_J;

    const double gamma_final = final_energy_J * (1.0 + alpha_per_J * final_energy_J);

    const double prefactor = std::sqrt(2.0) * static_cast<double>(branch.m_final_valley_count) * std::pow(mD, 1.5) *
                             (D0_J_per_m * D0_J_per_m) /
                             (uepm::constants::pi * rho_si * uepm::constants::h_bar * uepm::constants::h_bar * phonon_energy_J);

    return prefactor * phonon_factor * std::sqrt(std::max(0.0, gamma_final)) * (1.0 + 2.0 * alpha_per_J * final_energy_J);
}

double intervalley_first_order_rate(const valley_model&              valley,
                                    const intervalley_phonon_branch& branch,
                                    double                           initial_energy_eV,
                                    bool                             absorption,
                                    double                           temperature_K) {
    constexpr double rho_si = 2.329e3;  // kg / m^3

    const double phonon_energy_eV = branch.m_phonon_energy_eV;
    const double final_energy_eV  = absorption ? (initial_energy_eV + phonon_energy_eV) : (initial_energy_eV - phonon_energy_eV);

    if (final_energy_eV < 0.0) {
        return 0.0;
    }

    const double Nop           = bose_einstein_occupation(phonon_energy_eV, temperature_K);
    const double phonon_factor = absorption ? Nop : (Nop + 1.0);

    const double alpha_per_J = silicon_nonparabolicity_per_joule(valley);
    const double mD          = silicon_density_of_states_mass(valley);

    const double initial_energy_J = initial_energy_eV * uepm::constants::eV_to_J;
    const double final_energy_J   = final_energy_eV * uepm::constants::eV_to_J;
    const double phonon_energy_J  = phonon_energy_eV * uepm::constants::eV_to_J;

    const double D1_J = branch.m_deformation_potential_1 * uepm::constants::eV_to_J;

    const double gamma_initial = initial_energy_J * (1.0 + alpha_per_J * initial_energy_J);
    const double gamma_final   = final_energy_J * (1.0 + alpha_per_J * final_energy_J);

    const double prefactor = std::sqrt(2.0) * static_cast<double>(branch.m_final_valley_count) * std::pow(mD, 2.5) * (D1_J * D1_J) /
                             (uepm::constants::pi * rho_si * std::pow(uepm::constants::h_bar, 4) * phonon_energy_J);

    return prefactor * phonon_factor * std::sqrt(std::max(0.0, gamma_final)) * (1.0 + 2.0 * alpha_per_J * final_energy_J) *
           (gamma_final + gamma_initial);
}

double intervalley_scattering_rate(const valley_model&              valley,
                                   const intervalley_phonon_branch& branch,
                                   double                           initial_energy_eV,
                                   bool                             absorption,
                                   double                           temperature_K) {
    if (branch.is_zeroth_order()) {
        return intervalley_zeroth_order_rate(valley, branch, initial_energy_eV, absorption, temperature_K);
    }

    return intervalley_first_order_rate(valley, branch, initial_energy_eV, absorption, temperature_K);
}

double optical_scattering_rate_silicon_holes(const valley_model&            final_band,
                                             const hole_optical_transition& transition,
                                             double                         initial_energy_eV,
                                             bool                           absorption,
                                             double                         temperature_K) {
    constexpr double rho_si = 2.329e3;  // kg / m^3

    const double phonon_energy_eV = transition.phonon_energy_eV;
    const double final_energy_eV  = absorption ? (initial_energy_eV + phonon_energy_eV) : (initial_energy_eV - phonon_energy_eV);

    if (final_energy_eV < 0.0) {
        return 0.0;
    }

    const double n_op          = bose_einstein_occupation(phonon_energy_eV, temperature_K);
    const double phonon_factor = absorption ? n_op : (n_op + 1.0);

    const double alpha_per_J = final_band.non_parabolicity() / uepm::constants::eV_to_J;
    const double mt          = final_band.transverse_effective_mass();
    const double ml          = final_band.longitudinal_effective_mass();
    const double mD          = std::cbrt(mt * mt * ml);

    const double final_energy_J  = final_energy_eV * uepm::constants::eV_to_J;
    const double phonon_energy_J = phonon_energy_eV * uepm::constants::eV_to_J;
    const double dop_J_per_m     = transition.deformation_potential_eV_per_m * uepm::constants::eV_to_J;

    const double gamma_final = final_energy_J * (1.0 + alpha_per_J * final_energy_J);

    const double prefactor = std::sqrt(2.0) * std::pow(mD, 1.5) * (dop_J_per_m * dop_J_per_m) * transition.overlap_factor /
                             (uepm::constants::pi * rho_si * uepm::constants::h_bar * uepm::constants::h_bar * phonon_energy_J);

    return prefactor * phonon_factor * std::sqrt(std::max(0.0, gamma_final)) * (1.0 + 2.0 * alpha_per_J * final_energy_J);
}

} // namespace uepm::amc