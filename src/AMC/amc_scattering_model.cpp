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

#include <cmath>
#include <stdexcept>

#include "physical_constants.hpp"

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
    // constexpr double D_ac_eV = 6.6;  // eV (Valeur de Aubry-Fortuna et al. pour les électrons, ajustée pour mieux correspondre aux données expérimentales de mobilité)
    constexpr double D_ac_eV = 9.0;  // eV (Valeur dans l'HDR de Phillipe)

    const double energy_clamped_eV = std::max(energy_eV, 1.0e-9);
    const double energy_J          = energy_clamped_eV * uepm::constants::eV_to_J;

    const double alpha_per_J = valley.non_parabolicity() / uepm::constants::eV_to_J;

    const double mt = valley.transverse_effective_mass();
    const double ml = valley.longitudinal_effective_mass();
    const double mD = std::cbrt(mt * mt * ml);

    const double D_ac_J = D_ac_eV * uepm::constants::eV_to_J;

    const double gamma_factor =
        std::sqrt(energy_J * (1.0 + alpha_per_J * energy_J)) * (1.0 + 2.0 * alpha_per_J * energy_J);

    const double prefactor = std::sqrt(2.0) * uepm::constants::k_B * temperature_K * std::pow(mD, 1.5) *
                             (D_ac_J * D_ac_J) /
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

    const double gamma_factor =
        std::sqrt(energy_J * (1.0 + alpha_per_J * energy_J)) * (1.0 + 2.0 * alpha_per_J * energy_J);

    const double prefactor = std::sqrt(2.0) * uepm::constants::k_B * temperature_K * std::pow(mD, 1.5) *
                             (D_ac_J * D_ac_J) * overlap_factor /
                             (uepm::constants::pi * std::pow(uepm::constants::h_bar, 4) * rho_si * v_s * v_s);

    return prefactor * gamma_factor;
}

double silicon_density_of_states_mass(const valley_model& valley) {
    const double mt = valley.transverse_effective_mass();
    const double ml = valley.longitudinal_effective_mass();
    return std::cbrt(mt * mt * ml);
}

double silicon_nonparabolicity_per_joule(const valley_model& valley) {
    return valley.non_parabolicity() / uepm::constants::eV_to_J;
}

double intervalley_zeroth_order_rate(const valley_model&              valley,
                                     const intervalley_phonon_branch& branch,
                                     double                           initial_energy_eV,
                                     bool                             absorption,
                                     double                           temperature_K) {
    constexpr double rho_si = 2.329e3;  // kg / m^3

    const double phonon_energy_eV = branch.m_phonon_energy_eV;
    const double final_energy_eV =
        absorption ? (initial_energy_eV + phonon_energy_eV) : (initial_energy_eV - phonon_energy_eV);

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

    const double prefactor =
        std::sqrt(2.0) * static_cast<double>(branch.m_final_valley_count) * std::pow(mD, 1.5) *
        (D0_J_per_m * D0_J_per_m) /
        (uepm::constants::pi * rho_si * uepm::constants::h_bar * uepm::constants::h_bar * phonon_energy_J);

    return prefactor * phonon_factor * std::sqrt(std::max(0.0, gamma_final)) *
           (1.0 + 2.0 * alpha_per_J * final_energy_J);
}

double intervalley_first_order_rate(const valley_model&              valley,
                                    const intervalley_phonon_branch& branch,
                                    double                           initial_energy_eV,
                                    bool                             absorption,
                                    double                           temperature_K) {
    constexpr double rho_si = 2.329e3;  // kg / m^3

    const double phonon_energy_eV = branch.m_phonon_energy_eV;
    const double final_energy_eV =
        absorption ? (initial_energy_eV + phonon_energy_eV) : (initial_energy_eV - phonon_energy_eV);

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

    const double prefactor = std::sqrt(2.0) * static_cast<double>(branch.m_final_valley_count) * std::pow(mD, 2.5) *
                             (D1_J * D1_J) /
                             (uepm::constants::pi * rho_si * std::pow(uepm::constants::h_bar, 4) * phonon_energy_J);

    return prefactor * phonon_factor * std::sqrt(std::max(0.0, gamma_final)) *
           (1.0 + 2.0 * alpha_per_J * final_energy_J) * (gamma_final + gamma_initial);
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
    const double final_energy_eV =
        absorption ? (initial_energy_eV + phonon_energy_eV) : (initial_energy_eV - phonon_energy_eV);

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

    const double prefactor =
        std::sqrt(2.0) * std::pow(mD, 1.5) * (dop_J_per_m * dop_J_per_m) * transition.overlap_factor /
        (uepm::constants::pi * rho_si * uepm::constants::h_bar * uepm::constants::h_bar * phonon_energy_J);

    return prefactor * phonon_factor * std::sqrt(std::max(0.0, gamma_final)) *
           (1.0 + 2.0 * alpha_per_J * final_energy_J);
}

double caughey_thomas_mobility_cm2_per_V_s(double                              impurity_density_cm_3,
                                           const impurity_mobility_parameters& parameters) {
    if (impurity_density_cm_3 <= 0.0) {
        return parameters.m_mu0_cm2_per_V_s;
    }

    const double ratio = impurity_density_cm_3 / parameters.m_n_ref_cm_3;

    return parameters.m_mu_min_cm2_per_V_s + (parameters.m_mu0_cm2_per_V_s - parameters.m_mu_min_cm2_per_V_s) /
                                                 (1.0 + std::pow(ratio, parameters.m_alpha));
}

double impurity_momentum_relaxation_rate_silicon(const valley_model&                         band_or_valley,
                                                 particle_type                               carrier_type,
                                                 double                                      impurity_density_cm_3,
                                                 const carrier_impurity_mobility_parameters& parameters) {
    if (impurity_density_cm_3 <= 0.0) {
        return 0.0;
    }

    const auto& mobility_parameters =
        (carrier_type == particle_type::electron) ? parameters.m_electron : parameters.m_hole;

    const double mu_doped_cm2_per_V_s = caughey_thomas_mobility_cm2_per_V_s(impurity_density_cm_3, mobility_parameters);
    const double mu_lattice_cm2_per_V_s = mobility_parameters.m_mu0_cm2_per_V_s;
    const double inv_mu_impurity        = 1.0 / mu_doped_cm2_per_V_s - 1.0 / mu_lattice_cm2_per_V_s;

    if (inv_mu_impurity <= 0.0) {
        return 0.0;
    }

    const double mu_impurity_cm2_per_V_s = 1.0 / inv_mu_impurity;
    const double mu_impurity_m2_per_V_s  = mu_impurity_cm2_per_V_s * 1.0e-4;
    const double mt                      = band_or_valley.transverse_effective_mass();
    const double ml                      = band_or_valley.longitudinal_effective_mass();

    // Scalar approximation for the momentum relaxation mass.
    const double m_eff = std::cbrt(mt * mt * ml);

    return uepm::constants::q_e / (m_eff * mu_impurity_m2_per_V_s);
}

double screened_coulomb_impurity_momentum_relaxation_rate_silicon(const valley_model& band_or_valley,
                                                                  particle_type       carrier_type,
                                                                  double              energy_eV,
                                                                  double              impurity_density_cm_3,
                                                                  double              screening_density_cm_3,
                                                                  double              temperature_K) {
    if (energy_eV <= 0.0) {
        return 0.0;
    }

    if (impurity_density_cm_3 <= 0.0) {
        return 0.0;
    }

    if (screening_density_cm_3 <= 0.0) {
        screening_density_cm_3 = impurity_density_cm_3;
    }

    if (temperature_K <= 0.0) {
        throw std::invalid_argument("temperature must be positive for screened Coulomb impurity scattering");
    }

    constexpr double q_C               = uepm::constants::q_e;
    constexpr double k_B_J_per_K       = uepm::constants::k_B;
    constexpr double epsilon0_F_per_m  = uepm::constants::eps_0;
    constexpr double epsilon_r_silicon = 11.7;
    constexpr double pi                = uepm::constants::pi;

    const double epsilon_s_F_per_m = epsilon_r_silicon * epsilon0_F_per_m;

    const double impurity_density_m_3  = impurity_density_cm_3 * 1.0e6;
    const double screening_density_m_3 = screening_density_cm_3 * 1.0e6;

    const double kBT_J = k_B_J_per_K * temperature_K;
    const double E_J   = energy_eV * q_C;

    const double mt = band_or_valley.transverse_effective_mass();
    const double ml = band_or_valley.longitudinal_effective_mass();

    const double m_d = std::cbrt(mt * mt * ml);

    const double alpha_eV_inv = carrier_type == particle_type::electron ? 0.5 : 0.0;

    const double alpha_J_inv = alpha_eV_inv / q_C;

    const double gamma_J = E_J * (1.0 + alpha_J_inv * E_J);

    if (gamma_J <= 0.0) {
        return 0.0;
    }

    const double hbar = uepm::constants::h_bar;

    const double k2 = 2.0 * m_d * gamma_J / (hbar * hbar);

    if (k2 <= 0.0) {
        return 0.0;
    }

    const double q_screen2 = q_C * q_C * screening_density_m_3 / (epsilon_s_F_per_m * kBT_J);

    if (q_screen2 <= 0.0) {
        return 0.0;
    }

    /*
     * Momentum-relaxation angular integral:
     *
     * I = ∫_{-1}^{1} (1 - cosθ) d(cosθ)
     *     / (q_s² + 2 k² (1 - cosθ))²
     *
     * With beta = 4 k² / q_s²:
     *
     * I = [ln(1 + beta) - beta / (1 + beta)] / (4 k⁴)
     *
     * This is the Brooks-Herring / screened-Coulomb momentum-relaxation form.
     */
    const double beta = 4.0 * k2 / q_screen2;

    if (beta <= 0.0) {
        return 0.0;
    }

    const double angular_factor = std::log1p(beta) - beta / (1.0 + beta);

    if (angular_factor <= 0.0 || !std::isfinite(angular_factor)) {
        return 0.0;
    }

    const double angular_integral = angular_factor / (4.0 * k2 * k2);

    const double nonparabolic_factor = (1.0 + 2.0 * alpha_J_inv * E_J) * std::sqrt(gamma_J);

    const double prefactor = std::sqrt(2.0) * std::pow(q_C, 4) * std::pow(m_d, 1.5) * impurity_density_m_3 /
                             (pi * std::pow(hbar, 4) * epsilon_s_F_per_m * epsilon_s_F_per_m);

    const double rate = prefactor * nonparabolic_factor * angular_integral;

    if (!std::isfinite(rate) || rate < 0.0) {
        return 0.0;
    }

    return rate;
}

}  // namespace uepm::amc