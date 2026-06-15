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

#include <array>
#include <cmath>
#include <stdexcept>

#include "integrals.hpp"
#include "physical_constants.hpp"

namespace uepm::amc {

namespace {

double screening_function_by_quadrature(double xi) {
    const double xi2 = xi * xi;
    return uepm::integrate::integrate_gauss_legendre_32([xi2](double s) { return std::exp(-xi2 * (1.0 - s * s)); },
                                                        0.0,
                                                        1.0);
}

const std::array<double, 4097>& screening_function_table() {
    static const std::array<double, 4097> table = [] {
        std::array<double, 4097> values{};
        constexpr double         max_xi = 8.0;
        constexpr double         step   = max_xi / static_cast<double>(values.size() - 1);

        for (std::size_t i = 0; i < values.size(); ++i) {
            values[i] = screening_function_by_quadrature(static_cast<double>(i) * step);
        }
        return values;
    }();
    return table;
}

}  // namespace

double bose_einstein_occupation(double phonon_energy_eV, double temperature_K) {
    if (temperature_K <= 0.0) {
        return 0.0;
    }
    const double x = phonon_energy_eV / (uepm::constants::k_b_eV * temperature_K);
    return 1.0 / std::expm1(x);
}

double acoustic_scattering_rate(const valley_model&                   valley,
                                const acoustic_scattering_parameters& parameters,
                                double                                energy_eV,
                                double                                temperature_K) {
    if (temperature_K <= 0.0) {
        return 0.0;
    }

    const double energy_clamped_eV = std::max(energy_eV, 1.0e-9);
    const double energy_J          = energy_clamped_eV * uepm::constants::eV_to_J;

    const double alpha_per_J = valley.non_parabolicity() / uepm::constants::eV_to_J;

    const double mt = valley.transverse_effective_mass();
    const double ml = valley.longitudinal_effective_mass();
    const double mD = std::cbrt(mt * mt * ml);

    const double D_ac_J = parameters.deformation_potential_eV * uepm::constants::eV_to_J;

    const double gamma_factor =
        std::sqrt(energy_J * (1.0 + alpha_per_J * energy_J)) * (1.0 + 2.0 * alpha_per_J * energy_J);

    const double prefactor =
        std::sqrt(2.0) * uepm::constants::k_B * temperature_K * std::pow(mD, 1.5) * (D_ac_J * D_ac_J) *
        parameters.overlap_factor /
        (uepm::constants::pi * std::pow(uepm::constants::h_bar, 4) * parameters.mass_density_kg_per_m3 *
         parameters.sound_velocity_m_per_s * parameters.sound_velocity_m_per_s);

    return prefactor * gamma_factor;
}

double intervalley_zeroth_order_rate(const valley_model&              valley,
                                     const intervalley_phonon_branch& branch,
                                     double                           mass_density_kg_per_m3,
                                     double                           initial_energy_eV,
                                     bool                             absorption,
                                     double                           temperature_K) {
    const double phonon_energy_eV = branch.m_phonon_energy_eV;
    const double final_energy_eV =
        absorption ? (initial_energy_eV + phonon_energy_eV) : (initial_energy_eV - phonon_energy_eV);

    if (final_energy_eV < 0.0) {
        return 0.0;
    }

    const double Nop           = bose_einstein_occupation(phonon_energy_eV, temperature_K);
    const double phonon_factor = absorption ? Nop : (Nop + 1.0);

    const double alpha_per_J = valley.non_parabolicity() / uepm::constants::eV_to_J;
    const double mD          = valley.density_of_states_effective_mass();

    const double final_energy_J  = final_energy_eV * uepm::constants::eV_to_J;
    const double phonon_energy_J = phonon_energy_eV * uepm::constants::eV_to_J;

    const double D0_J_per_m = branch.m_deformation_potential_0 * uepm::constants::eV_to_J;

    const double gamma_final = final_energy_J * (1.0 + alpha_per_J * final_energy_J);

    const double prefactor = std::sqrt(2.0) * static_cast<double>(branch.m_final_valley_count) * std::pow(mD, 1.5) *
                             (D0_J_per_m * D0_J_per_m) /
                             (uepm::constants::pi * mass_density_kg_per_m3 * uepm::constants::h_bar *
                              uepm::constants::h_bar * phonon_energy_J);

    return prefactor * phonon_factor * std::sqrt(std::max(0.0, gamma_final)) *
           (1.0 + 2.0 * alpha_per_J * final_energy_J);
}

double intervalley_first_order_rate(const valley_model&              valley,
                                    const intervalley_phonon_branch& branch,
                                    double                           mass_density_kg_per_m3,
                                    double                           initial_energy_eV,
                                    bool                             absorption,
                                    double                           temperature_K) {
    const double phonon_energy_eV = branch.m_phonon_energy_eV;
    const double final_energy_eV =
        absorption ? (initial_energy_eV + phonon_energy_eV) : (initial_energy_eV - phonon_energy_eV);

    if (final_energy_eV < 0.0) {
        return 0.0;
    }

    const double Nop           = bose_einstein_occupation(phonon_energy_eV, temperature_K);
    const double phonon_factor = absorption ? Nop : (Nop + 1.0);

    const double alpha_per_J = valley.non_parabolicity() / uepm::constants::eV_to_J;
    const double mD          = valley.density_of_states_effective_mass();

    const double initial_energy_J = initial_energy_eV * uepm::constants::eV_to_J;
    const double final_energy_J   = final_energy_eV * uepm::constants::eV_to_J;
    const double phonon_energy_J  = phonon_energy_eV * uepm::constants::eV_to_J;

    const double D1_J = branch.m_deformation_potential_1 * uepm::constants::eV_to_J;

    const double gamma_initial = initial_energy_J * (1.0 + alpha_per_J * initial_energy_J);
    const double gamma_final   = final_energy_J * (1.0 + alpha_per_J * final_energy_J);

    const double prefactor =
        std::sqrt(2.0) * static_cast<double>(branch.m_final_valley_count) * std::pow(mD, 2.5) * (D1_J * D1_J) /
        (uepm::constants::pi * mass_density_kg_per_m3 * std::pow(uepm::constants::h_bar, 4) * phonon_energy_J);

    return prefactor * phonon_factor * std::sqrt(std::max(0.0, gamma_final)) *
           (1.0 + 2.0 * alpha_per_J * final_energy_J) * (gamma_final + gamma_initial);
}

double intervalley_scattering_rate(const valley_model&              valley,
                                   const intervalley_phonon_branch& branch,
                                   double                           mass_density_kg_per_m3,
                                   double                           initial_energy_eV,
                                   bool                             absorption,
                                   double                           temperature_K) {
    if (branch.is_zeroth_order()) {
        return intervalley_zeroth_order_rate(valley,
                                             branch,
                                             mass_density_kg_per_m3,
                                             initial_energy_eV,
                                             absorption,
                                             temperature_K);
    }

    return intervalley_first_order_rate(valley,
                                        branch,
                                        mass_density_kg_per_m3,
                                        initial_energy_eV,
                                        absorption,
                                        temperature_K);
}

double optical_scattering_rate_holes(const valley_model&            final_band,
                                     const hole_optical_transition& transition,
                                     double                         mass_density_kg_per_m3,
                                     double                         initial_energy_eV,
                                     bool                           absorption,
                                     double                         temperature_K) {
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

    const double prefactor = std::sqrt(2.0) * std::pow(mD, 1.5) * (dop_J_per_m * dop_J_per_m) *
                             transition.overlap_factor /
                             (uepm::constants::pi * mass_density_kg_per_m3 * uepm::constants::h_bar *
                              uepm::constants::h_bar * phonon_energy_J);

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

double impurity_momentum_relaxation_rate(const valley_model&                         band_or_valley,
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

double impurity_screening_function(double xi) {
    const double x = std::abs(xi);
    if (!std::isfinite(x)) {
        return 0.0;
    }
    if (x < 1.0e-4) {
        const double x2 = x * x;
        return 1.0 - (2.0 / 3.0) * x2 + (4.0 / 15.0) * x2 * x2 - (8.0 / 105.0) * x2 * x2 * x2;
    }
    if (x >= 8.0) {
        const double inverse_x2 = 1.0 / (x * x);
        return 0.5 * inverse_x2 + 0.25 * inverse_x2 * inverse_x2 + 0.375 * inverse_x2 * inverse_x2 * inverse_x2 +
               0.9375 * inverse_x2 * inverse_x2 * inverse_x2 * inverse_x2;
    }

    constexpr double max_xi       = 8.0;
    const auto&      table        = screening_function_table();
    const double     scaled_index = x * static_cast<double>(table.size() - 1) / max_xi;
    const auto       lower_index  = static_cast<std::size_t>(scaled_index);
    const double     fraction     = scaled_index - static_cast<double>(lower_index);

    return table[lower_index] + fraction * (table[lower_index + 1] - table[lower_index]);
}

double analytic_brooks_herring_momentum_integral(double k2, double q_screen2) {
    if (k2 <= 0.0 || q_screen2 <= 0.0) {
        return 0.0;
    }
    const double beta           = 4.0 * k2 / q_screen2;
    const double angular_factor = std::log1p(beta) - beta / (1.0 + beta);
    if (angular_factor <= 0.0 || !std::isfinite(angular_factor)) {
        return 0.0;
    }
    return angular_factor / (4.0 * k2 * k2);
}

double full_screening_momentum_integral(double k2, double q_screen2, double gamma_J, double kBT_J) {
    if (k2 <= 0.0 || q_screen2 <= 0.0 || gamma_J <= 0.0 || kBT_J <= 0.0) {
        return 0.0;
    }

    // Resolve the forward-scattering region with u = u0 * (exp(y) - 1),
    // where u0 is the natural Debye screening scale.
    const double u0       = q_screen2 / (2.0 * k2);
    const double y_max    = std::log1p(2.0 / u0);
    const double integral = uepm::integrate::integrate_gauss_legendre_32(
        [=](double y) {
            const double exp_y       = std::exp(y);
            const double u           = u0 * (exp_y - 1.0);
            const double du_dy       = u0 * exp_y;
            const double xi          = std::sqrt(gamma_J * u / (2.0 * kBT_J));
            const double denominator = 2.0 * k2 * u + q_screen2 * impurity_screening_function(xi);
            return u * du_dy / (denominator * denominator);
        },
        0.0,
        y_max);

    return std::isfinite(integral) && integral > 0.0 ? integral : 0.0;
}

double screened_coulomb_impurity_momentum_relaxation_rate(const valley_model&      band_or_valley,
                                                          double                   relative_permittivity,
                                                          double                   energy_eV,
                                                          double                   impurity_density_cm_3,
                                                          double                   screening_density_cm_3,
                                                          double                   temperature_K,
                                                          impurity_screening_model screening_model) {
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
    constexpr double q_C                   = uepm::constants::q_e;
    constexpr double k_B_J_per_K           = uepm::constants::k_B;
    constexpr double epsilon0_F_per_m      = uepm::constants::eps_0;
    constexpr double pi                    = uepm::constants::pi;
    const double     epsilon_s_F_per_m     = relative_permittivity * epsilon0_F_per_m;
    const double     impurity_density_m_3  = impurity_density_cm_3 * 1.0e6;
    const double     screening_density_m_3 = screening_density_cm_3 * 1.0e6;
    const double     kBT_J                 = k_B_J_per_K * temperature_K;
    const double     E_J                   = energy_eV * q_C;
    const double     mt                    = band_or_valley.transverse_effective_mass();
    const double     ml                    = band_or_valley.longitudinal_effective_mass();
    const double     m_d                   = std::cbrt(mt * mt * ml);
    const double     alpha_eV_inv          = band_or_valley.non_parabolicity();
    const double     alpha_J_inv           = alpha_eV_inv / q_C;
    const double     gamma_J               = E_J * (1.0 + alpha_J_inv * E_J);
    if (gamma_J <= 0.0) {
        return 0.0;
    }
    const double hbar = uepm::constants::h_bar;
    const double k2   = 2.0 * m_d * gamma_J / (hbar * hbar);

    if (k2 <= 0.0) {
        return 0.0;
    }
    const double q_screen2 = q_C * q_C * screening_density_m_3 / (epsilon_s_F_per_m * kBT_J);
    if (q_screen2 <= 0.0) {
        return 0.0;
    }
    const double angular_integral = screening_model == impurity_screening_model::debye_analytic
                                        ? analytic_brooks_herring_momentum_integral(k2, q_screen2)
                                        : full_screening_momentum_integral(k2, q_screen2, gamma_J, kBT_J);
    if (angular_integral <= 0.0) {
        return 0.0;
    }
    const double nonparabolic_factor = (1.0 + 2.0 * alpha_J_inv * E_J) * std::sqrt(gamma_J);
    const double prefactor           = std::sqrt(2.0) * std::pow(q_C, 4) * std::pow(m_d, 1.5) * impurity_density_m_3 /
                             (pi * std::pow(hbar, 4) * epsilon_s_F_per_m * epsilon_s_F_per_m);

    const double rate = prefactor * nonparabolic_factor * angular_integral;

    if (!std::isfinite(rate) || rate < 0.0) {
        return 0.0;
    }

    return rate;
}

}  // namespace uepm::amc
