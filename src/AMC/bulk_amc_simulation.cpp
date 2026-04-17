/**
 * @file bulk_amc_simulation.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-04-10
 *
 *
 */

#include "bulk_amc_simulation.hpp"

#include <fmt/core.h>
#include <fmt/format.h>

#include <array>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>

#include "intervalley_phonon.hpp"

namespace uepm::amc {

constexpr double m0 = uepm::constants::m_e;

// Local z -> global z
valley_model::mat3 rotation_local_z_to_global_z() { return valley_model::identity_matrix(); }

// Local z -> global x
valley_model::mat3 rotation_local_z_to_global_x() { return {{{{0.0, 0.0, 1.0}}, {{0.0, 1.0, 0.0}}, {{-1.0, 0.0, 0.0}}}}; }

// Local z -> global y
valley_model::mat3 rotation_local_z_to_global_y() { return {{{{1.0, 0.0, 0.0}}, {{0.0, 0.0, 1.0}}, {{0.0, -1.0, 0.0}}}}; }

std::vector<valley_model> make_silicon_delta_valleys() {
    std::vector<valley_model> valleys;
    valleys.reserve(6);

    valley_model::parameters p;
    p.transverse_effective_mass   = 0.19 * m0;
    p.longitudinal_effective_mass = 0.916 * m0;

    // Keep this simple for now.
    // Use parabolic first; switch to Kane once the loop works.
    p.dispersion = valley_model::band_type::parabolic;

    // If you want Kane immediately, set:
    // p.dispersion = valley_model::band_type::kane;
    // p.non_parabolicity = ...; // in 1/eV or 1/J depending on your convention

    p.non_parabolicity        = 0.0;
    p.energy_offset           = 0.0;
    p.phonon_reference_energy = 0.0;
    p.degeneracy              = 1;

    p.name     = "Delta_x_plus";
    p.rotation = rotation_local_z_to_global_x();
    valleys.emplace_back(p);

    p.name     = "Delta_x_minus";
    p.rotation = rotation_local_z_to_global_x();
    valleys.emplace_back(p);

    p.name     = "Delta_y_plus";
    p.rotation = rotation_local_z_to_global_y();
    valleys.emplace_back(p);

    p.name     = "Delta_y_minus";
    p.rotation = rotation_local_z_to_global_y();
    valleys.emplace_back(p);

    p.name     = "Delta_z_plus";
    p.rotation = rotation_local_z_to_global_z();
    valleys.emplace_back(p);

    p.name     = "Delta_z_minus";
    p.rotation = rotation_local_z_to_global_z();
    valleys.emplace_back(p);

    return valleys;
}

std::vector<valley_model> make_silicon_hole_bands() {
    std::vector<valley_model> bands;
    bands.reserve(2);

    valley_model::parameters p;
    p.dispersion              = valley_model::band_type::parabolic;
    p.non_parabolicity        = 0.0;
    p.energy_offset           = 0.0;
    p.phonon_reference_energy = 0.0;
    p.degeneracy              = 1;
    p.rotation                = valley_model::identity_matrix();

    // Heavy-hole band: isotropic parabolic approximation
    p.name                        = "heavy_hole";
    p.transverse_effective_mass   = 0.87 * m0;
    p.longitudinal_effective_mass = 0.87 * m0;
    bands.emplace_back(p);

    // Light-hole band: isotropic parabolic approximation
    p.name                        = "light_hole";
    p.transverse_effective_mass   = 0.24 * m0;
    p.longitudinal_effective_mass = 0.24 * m0;
    bands.emplace_back(p);

    return bands;
}

double sample_thermal_energy_eV(double temperature_K, std::mt19937_64& rng) {
    const double                    kT_eV = uepm::constants::k_b_eV * temperature_K;
    std::gamma_distribution<double> dist(1.5, kT_eV);
    return dist(rng);
}

std::size_t valley_axis(std::size_t valley_index) {
    if (valley_index < 2) {
        return 0;
    }
    if (valley_index < 4) {
        return 1;
    }
    if (valley_index < 6) {
        return 2;
    }
    throw std::out_of_range("invalid valley index in valley_axis");
}

std::size_t draw_g_destination_valley(std::size_t current_valley_index) {
    switch (current_valley_index) {
        case 0:
            return 1;
        case 1:
            return 0;
        case 2:
            return 3;
        case 3:
            return 2;
        case 4:
            return 5;
        case 5:
            return 4;
        default:
            throw std::out_of_range("invalid valley index in draw_g_destination_valley");
    }
}

std::size_t draw_f_destination_valley(std::size_t current_valley_index, std::mt19937_64& rng) {
    const std::size_t current_axis = valley_axis(current_valley_index);

    std::array<std::size_t, 4> candidates{};
    std::size_t                count = 0;

    for (std::size_t valley_index = 0; valley_index < 6; ++valley_index) {
        if (valley_axis(valley_index) != current_axis) {
            candidates[count++] = valley_index;
        }
    }

    if (count != 4) {
        throw std::runtime_error("unexpected number of f-type destination valleys");
    }

    std::uniform_int_distribution<std::size_t> dist(0, 3);
    return candidates[dist(rng)];
}

std::size_t draw_intervalley_destination_valley(const intervalley_phonon_branch& branch,
                                                std::size_t                      current_valley_index,
                                                std::mt19937_64&                 rng) {
    switch (branch.m_family) {
        case intervalley_family::g:
            return draw_g_destination_valley(current_valley_index);
        case intervalley_family::f:
            return draw_f_destination_valley(current_valley_index, rng);
        default:
            throw std::runtime_error("unknown intervalley family");
    }
}

std::vector<hole_optical_transition> make_silicon_hole_optical_transitions() {
    constexpr double eV_per_cm_to_eV_per_m = 100.0;

    const double dop_eV_per_m = 8.0e8 * eV_per_cm_to_eV_per_m;
    const double eop_eV       = 63.0e-3;

    std::vector<hole_optical_transition> transitions;
    transitions.reserve(4);

    transitions.push_back({"hh_to_hh", 0, 0, eop_eV, dop_eV_per_m, 0.5});
    transitions.push_back({"hh_to_lh", 0, 1, eop_eV, dop_eV_per_m, 1.0});
    transitions.push_back({"lh_to_hh", 1, 0, eop_eV, dop_eV_per_m, 1.0});
    transitions.push_back({"lh_to_lh", 1, 1, eop_eV, dop_eV_per_m, 0.5});

    return transitions;
}

void bulk_amc_simulation::initialize() {
    if (m_cfg.m_carrier_type == particle_type::electron) {
        m_valleys              = make_silicon_delta_valleys();
        m_intervalley_branches = make_silicon_intervalley_phonon_branches();
        m_hole_optical_transitions.clear();
    } else {
        m_valleys = make_silicon_hole_bands();
        m_intervalley_branches.clear();
        m_hole_optical_transitions = make_silicon_hole_optical_transitions();
    }

    m_particles.clear();
    m_particles.reserve(m_cfg.m_number_of_particles);

    for (std::size_t i = 0; i < m_cfg.m_number_of_particles; ++i) {
        particle_amc p{i, m_cfg.m_carrier_type, 1.0};

        p.state().time         = 0.0;
        p.state().position     = vector3{0.0, 0.0, 0.0};
        p.state().valley_index = i % m_valleys.size();

        const auto& valley = m_valleys[p.state().valley_index];

        const double energy_eV   = sample_thermal_energy_eV(m_cfg.m_lattice_temperature, m_rng);
        p.state().local_k        = valley.draw_random_k_valley_at_energy(energy_eV, m_rng);
        p.state().gamma          = valley.gamma_from_k_valley(p.state().local_k);
        p.state().kinetic_energy = valley.kinetic_energy_from_gamma(p.state().gamma);
        p.state().velocity       = valley.velocity_from_k_valley(p.state().local_k);

        if (m_cfg.m_record_history) {
            p.record_state();
        }

        m_particles.push_back(std::move(p));
    }

    m_gamma_max_s_1 = compute_max_self_scattering_rate(m_cfg.m_max_energy_eV, m_cfg.m_gamma_max_energy_samples);
    fmt::print("Computed max self-scattering rate: {:.6e} s^-1\n", m_gamma_max_s_1);

    fmt::print("Computed max self-scattering rate: {:.6e} s^-1\n", m_gamma_max_s_1);

    std::cout << "Initialized " << m_particles.size() << " particles\n";
}

void bulk_amc_simulation::drift_particle(particle_amc& p, double dt) {
    if (dt < 0.0) {
        throw std::invalid_argument("drift time step must be non-negative");
    }

    const auto valley_index = p.state().valley_index;
    if (valley_index >= m_valleys.size()) {
        throw std::out_of_range("invalid valley index in drift_particle");
    }

    const auto&   valley                = m_valleys[valley_index];
    const vector3 electric_field_valley = valley.to_valley_frame(m_cfg.m_electric_field);

    const double  prefactor    = p.signed_charge() / uepm::constants::h_bar;
    const vector3 old_velocity = p.state().velocity;

    p.state().local_k += electric_field_valley * (prefactor * dt);

    p.state().gamma          = valley.gamma_from_k_valley(p.state().local_k);
    p.state().kinetic_energy = valley.kinetic_energy_from_gamma(p.state().gamma);
    p.state().velocity       = valley.velocity_from_k_valley(p.state().local_k);

    const vector3 avg_velocity = 0.5 * (old_velocity + p.state().velocity);
    p.state().position += avg_velocity * dt;
    p.state().time += dt;
}

double bose_einstein_occupation(double phonon_energy_eV, double temperature_K) {
    if (temperature_K <= 0.0) {
        return 0.0;
    }

    const double x = phonon_energy_eV / (uepm::constants::k_b_eV * temperature_K);
    return 1.0 / (std::exp(x) - 1.0);
}

double bulk_amc_simulation::sample_free_flight_time() {
    if (m_gamma_max_s_1 <= 0.0) {
        throw std::invalid_argument("max self-scattering rate must be > 0");
    }

    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    double                                 u = 0.0;

    do {
        u = unif01(m_rng);
    } while (u <= 0.0);

    return -std::log(u) / m_gamma_max_s_1;
}

scattering_channel bulk_amc_simulation::select_scattering_channel(const particle_amc& p) {
    const auto channels = build_scattering_channels(p);

    double total_rate = 0.0;
    for (const auto& channel : channels) {
        total_rate += channel.rate_s_1;
    }

    if (total_rate <= 0.0) {
        throw std::runtime_error("cannot select scattering channel with zero total rate");
    }

    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    const double                           r_select = unif01(m_rng) * total_rate;

    double cumulative = 0.0;
    for (const auto& channel : channels) {
        cumulative += channel.rate_s_1;
        if (r_select < cumulative) {
            return channel;
        }
    }

    throw std::runtime_error("failed to select a real scattering channel");
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

// double bose_einstein_occupation(double phonon_energy_eV, double temperature_K) {
//     if (temperature_K <= 0.0) {
//         return 0.0;
//     }

//     const double x = phonon_energy_eV / (uepm::constants::k_b_eV * temperature_K);
//     return 1.0 / (std::exp(x) - 1.0);
// }

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

std::vector<scattering_channel> bulk_amc_simulation::build_scattering_channels(const particle_amc& p) const {
    const auto current_band_index = p.state().valley_index;
    if (current_band_index >= m_valleys.size()) {
        throw std::out_of_range("invalid band/valley index in build_scattering_channels");
    }

    const auto&  current_band = m_valleys[current_band_index];
    const double energy_eV    = p.state().kinetic_energy;

    std::vector<scattering_channel> channels;

    if (p.type() == particle_type::hole) {
        channels.reserve(1 + 2 * m_hole_optical_transitions.size());

        const double acoustic_rate = acoustic_scattering_rate_silicon_holes(current_band, energy_eV, m_cfg.m_lattice_temperature);

        if (acoustic_rate > 0.0) {
            channels.push_back(scattering_channel{.mechanism         = scattering_mechanism::acoustic,
                                                  .rate_s_1          = acoustic_rate,
                                                  .final_energy_eV   = energy_eV,
                                                  .destination_index = current_band_index,
                                                  .branch            = nullptr,
                                                  .process           = intervalley_process::none});
        }

        for (const auto& transition : m_hole_optical_transitions) {
            if (transition.initial_band != current_band_index) {
                continue;
            }

            const auto& final_band = m_valleys[transition.final_band];

            const double rate_abs =
                optical_scattering_rate_silicon_holes(final_band, transition, energy_eV, true, m_cfg.m_lattice_temperature);

            if (rate_abs > 0.0) {
                channels.push_back(scattering_channel{.mechanism         = scattering_mechanism::intervalley,
                                                      .rate_s_1          = rate_abs,
                                                      .final_energy_eV   = energy_eV + transition.phonon_energy_eV,
                                                      .destination_index = transition.final_band,
                                                      .branch            = nullptr,
                                                      .process           = intervalley_process::absorption});
            }

            const double final_energy_emission_eV = energy_eV - transition.phonon_energy_eV;
            const double rate_em =
                optical_scattering_rate_silicon_holes(final_band, transition, energy_eV, false, m_cfg.m_lattice_temperature);

            if (rate_em > 0.0 && final_energy_emission_eV >= 0.0) {
                channels.push_back(scattering_channel{.mechanism         = scattering_mechanism::intervalley,
                                                      .rate_s_1          = rate_em,
                                                      .final_energy_eV   = final_energy_emission_eV,
                                                      .destination_index = transition.final_band,
                                                      .branch            = nullptr,
                                                      .process           = intervalley_process::emission});
            }
        }

        return channels;
    }

    channels.reserve(1 + 2 * m_intervalley_branches.size());

    const double acoustic_rate = acoustic_scattering_rate_silicon(current_band, energy_eV, m_cfg.m_lattice_temperature);

    if (acoustic_rate > 0.0) {
        channels.push_back(scattering_channel{.mechanism       = scattering_mechanism::acoustic,
                                              .rate_s_1        = acoustic_rate,
                                              .final_energy_eV = energy_eV,
                                              .branch          = nullptr,
                                              .process         = intervalley_process::none});
    }

    for (const auto& branch : m_intervalley_branches) {
        const double rate_abs = intervalley_scattering_rate(current_band, branch, energy_eV, true, m_cfg.m_lattice_temperature);

        if (rate_abs > 0.0) {
            channels.push_back(scattering_channel{.mechanism       = scattering_mechanism::intervalley,
                                                  .rate_s_1        = rate_abs,
                                                  .final_energy_eV = energy_eV + branch.m_phonon_energy_eV,
                                                  .branch          = &branch,
                                                  .process         = intervalley_process::absorption});
        }

        const double rate_em = intervalley_scattering_rate(current_band, branch, energy_eV, false, m_cfg.m_lattice_temperature);

        const double final_energy_emission_eV = energy_eV - branch.m_phonon_energy_eV;
        if (rate_em > 0.0 && final_energy_emission_eV >= 0.0) {
            channels.push_back(scattering_channel{.mechanism       = scattering_mechanism::intervalley,
                                                  .rate_s_1        = rate_em,
                                                  .final_energy_eV = final_energy_emission_eV,
                                                  .branch          = &branch,
                                                  .process         = intervalley_process::emission});
        }
    }

    return channels;
}

double bulk_amc_simulation::total_scattering_rate(const particle_amc& p) const {
    const auto channels = build_scattering_channels(p);

    double total_rate = 0.0;
    for (const auto& channel : channels) {
        total_rate += channel.rate_s_1;
    }

    return total_rate;
}

double bulk_amc_simulation::total_scattering_rate_for_energy(std::size_t band_or_valley_index, double energy_eV) const {
    if (band_or_valley_index >= m_valleys.size()) {
        throw std::out_of_range("invalid band/valley index in total_scattering_rate_for_energy");
    }

    const auto& band_or_valley = m_valleys[band_or_valley_index];
    double      total_rate     = 0.0;

    if (m_cfg.m_carrier_type == particle_type::hole) {
        total_rate += acoustic_scattering_rate_silicon_holes(band_or_valley, energy_eV, m_cfg.m_lattice_temperature);

        for (const auto& transition : m_hole_optical_transitions) {
            if (transition.initial_band != band_or_valley_index) {
                continue;
            }

            const auto& final_band = m_valleys[transition.final_band];

            total_rate += optical_scattering_rate_silicon_holes(final_band, transition, energy_eV, true, m_cfg.m_lattice_temperature);

            total_rate += optical_scattering_rate_silicon_holes(final_band, transition, energy_eV, false, m_cfg.m_lattice_temperature);
        }

        return total_rate;
    }

    total_rate += acoustic_scattering_rate_silicon(band_or_valley, energy_eV, m_cfg.m_lattice_temperature);

    for (const auto& branch : m_intervalley_branches) {
        total_rate += intervalley_scattering_rate(band_or_valley, branch, energy_eV, true, m_cfg.m_lattice_temperature);
        total_rate += intervalley_scattering_rate(band_or_valley, branch, energy_eV, false, m_cfg.m_lattice_temperature);
    }

    return total_rate;
}


double bulk_amc_simulation::compute_max_self_scattering_rate(double max_energy_eV, std::size_t n_samples) const {
    if (max_energy_eV <= 0.0) {
        throw std::invalid_argument("max energy must be > 0");
    }

    if (n_samples < 2) {
        throw std::invalid_argument("number of gamma-max energy samples must be >= 2");
    }

    double gamma_max = 0.0;

    for (std::size_t valley_index = 0; valley_index < m_valleys.size(); ++valley_index) {
        for (std::size_t i = 0; i < n_samples; ++i) {
            const double x         = static_cast<double>(i) / static_cast<double>(n_samples - 1);
            const double energy_eV = x * max_energy_eV;

            const double total_rate = total_scattering_rate_for_energy(valley_index, energy_eV);

            gamma_max = std::max(gamma_max, total_rate);
        }
    }

    return gamma_max * m_cfg.m_self_scattering_safety_factor;
}
void bulk_amc_simulation::apply_scattering_channel(particle_amc& p, const scattering_channel& channel) {
    if (channel.rate_s_1 < 0.0) {
        throw std::invalid_argument("negative scattering channel rate");
    }

    switch (channel.mechanism) {
        case scattering_mechanism::acoustic: {
            const auto band_or_valley_index = p.state().valley_index;
            if (band_or_valley_index >= m_valleys.size()) {
                throw std::out_of_range("invalid band/valley index in apply_scattering_channel acoustic");
            }

            const auto& band_or_valley = m_valleys[band_or_valley_index];

            p.state().local_k        = band_or_valley.draw_random_k_valley_at_energy(channel.final_energy_eV, m_rng);
            p.state().gamma          = band_or_valley.gamma_from_k_valley(p.state().local_k);
            p.state().kinetic_energy = band_or_valley.kinetic_energy_from_gamma(p.state().gamma);
            p.state().velocity       = band_or_valley.velocity_from_k_valley(p.state().local_k);

            p.increment_scattering_event_count();
            p.add_scattering_event(scattering_event::acoustic);
            return;
        }

        case scattering_mechanism::intervalley: {
            if (p.type() == particle_type::hole) {
                if (channel.destination_index >= m_valleys.size()) {
                    throw std::out_of_range("invalid destination band in apply_scattering_channel for holes");
                }

                const auto& dst_band = m_valleys[channel.destination_index];

                p.state().valley_index   = channel.destination_index;
                p.state().local_k        = dst_band.draw_random_k_valley_at_energy(channel.final_energy_eV, m_rng);
                p.state().gamma          = dst_band.gamma_from_k_valley(p.state().local_k);
                p.state().kinetic_energy = dst_band.kinetic_energy_from_gamma(p.state().gamma);
                p.state().velocity       = dst_band.velocity_from_k_valley(p.state().local_k);

                p.increment_scattering_event_count();

                if (channel.process == intervalley_process::absorption) {
                    p.add_scattering_event(scattering_event::intervalley_absorption);
                } else if (channel.process == intervalley_process::emission) {
                    p.add_scattering_event(scattering_event::intervalley_emission);
                } else {
                    throw std::runtime_error("hole optical channel missing absorption/emission tag");
                }

                return;
            }

            if (channel.branch == nullptr) {
                throw std::runtime_error("electron intervalley channel missing branch");
            }

            const auto current_valley_index = p.state().valley_index;
            if (current_valley_index >= m_valleys.size()) {
                throw std::out_of_range("invalid current valley index in apply_scattering_channel for electrons");
            }

            const std::size_t destination_valley = draw_intervalley_destination_valley(*channel.branch, current_valley_index, m_rng);

            if (destination_valley >= m_valleys.size()) {
                throw std::out_of_range("invalid destination valley in apply_scattering_channel for electrons");
            }

            const auto& dst_valley = m_valleys[destination_valley];

            p.state().valley_index   = destination_valley;
            p.state().local_k        = dst_valley.draw_random_k_valley_at_energy(channel.final_energy_eV, m_rng);
            p.state().gamma          = dst_valley.gamma_from_k_valley(p.state().local_k);
            p.state().kinetic_energy = dst_valley.kinetic_energy_from_gamma(p.state().gamma);
            p.state().velocity       = dst_valley.velocity_from_k_valley(p.state().local_k);

            p.increment_scattering_event_count();

            if (channel.process == intervalley_process::absorption) {
                p.add_scattering_event(scattering_event::intervalley_absorption);
            } else if (channel.process == intervalley_process::emission) {
                p.add_scattering_event(scattering_event::intervalley_emission);
            } else {
                throw std::runtime_error("electron intervalley channel missing absorption/emission tag");
            }

            return;
        }
    }

    throw std::runtime_error("unknown scattering mechanism");
}

void bulk_amc_simulation::scatter_particle(particle_amc& p, double dt) {
    if (dt < 0.0) {
        throw std::invalid_argument("scatter time step must be non-negative");
    }

    const auto channels = build_scattering_channels(p);

    double total_rate = 0.0;
    for (const auto& channel : channels) {
        total_rate += channel.rate_s_1;
    }

    if (total_rate <= 0.0) {
        return;
    }

    const double scatter_probability = 1.0 - std::exp(-total_rate * dt);

    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    if (unif01(m_rng) >= scatter_probability) {
        return;
    }

    const double r_select = unif01(m_rng) * total_rate;

    double cumulative = 0.0;
    for (const auto& channel : channels) {
        cumulative += channel.rate_s_1;
        if (r_select < cumulative) {
            apply_scattering_channel(p, channel);
            return;
        }
    }

    throw std::runtime_error("failed to select a scattering channel");
}

void bulk_amc_simulation::accumulate_observables(double dt) {
    if (dt <= 0.0) {
        return;
    }

    for (const auto& p : m_particles) {
        m_observables.weighted_velocity_x_m2_per_s2 += p.state().velocity.x() * dt;
        m_observables.weighted_kinetic_energy_eV_s += p.state().kinetic_energy * dt;
    }

    m_observables.accumulated_time_s += static_cast<double>(m_particles.size()) * dt;
}

void bulk_amc_simulation::accumulate_particle_observables(const particle_amc& p, double dt) {
    if (dt <= 0.0) {
        return;
    }

    m_observables.weighted_velocity_x_m2_per_s2 += p.state().velocity.x() * dt;
    m_observables.weighted_kinetic_energy_eV_s += p.state().kinetic_energy * dt;
    m_observables.accumulated_time_s += dt;
}

void bulk_amc_simulation::run() {
    if (m_particles.empty()) {
        throw std::runtime_error("simulation not initialized");
    }

    const double      dt           = m_cfg.m_time_step;
    const std::size_t n_steps      = static_cast<std::size_t>(std::ceil(m_cfg.m_final_time / dt));
    const std::size_t warmup_steps = static_cast<std::size_t>(m_cfg.m_warmup_fraction * n_steps);

    m_observables = {};

    for (std::size_t step = 0; step < n_steps; ++step) {
        for (auto& p : m_particles) {
            drift_particle(p, dt);
            scatter_particle(p, dt);

            if (m_cfg.m_record_history) {
                p.record_state();
            }
        }

        if (step >= warmup_steps) {
            accumulate_observables(dt);
        }
    }

    fmt::print("Completed {} steps of {} particles\n", n_steps, m_particles.size());

    double  avg_energy = 0.0;
    vector3 avg_velocity{0.0, 0.0, 0.0};

    for (const auto& p : m_particles) {
        avg_energy += p.state().kinetic_energy;
        avg_velocity += p.state().velocity;
    }

    avg_energy /= static_cast<double>(m_particles.size());
    avg_velocity /= static_cast<double>(m_particles.size());

    fmt::print("Average kinetic energy: {:.6f} eV\n", avg_energy);
    fmt::print("Average velocity: ({:.6e}, {:.6e}, {:.6e}) m/s\n", avg_velocity.x(), avg_velocity.y(), avg_velocity.z());

    std::size_t total_acoustic_events               = 0;
    std::size_t total_intervalley_absorption_events = 0;
    std::size_t total_intervalley_emission_events   = 0;
    std::size_t total_self_scattering_events        = 0;

    for (const auto& p : m_particles) {
        const auto& events = p.history().scattering_events();

        total_acoustic_events += events[static_cast<std::size_t>(scattering_event::acoustic)];
        total_intervalley_absorption_events += events[static_cast<std::size_t>(scattering_event::intervalley_absorption)];
        total_intervalley_emission_events += events[static_cast<std::size_t>(scattering_event::intervalley_emission)];
        total_self_scattering_events += events[static_cast<std::size_t>(scattering_event::self_scattering)];
    }

    fmt::print("Total acoustic events: {}\n", total_acoustic_events);
    fmt::print("Total intervalley absorption events: {}\n", total_intervalley_absorption_events);
    fmt::print("Total intervalley emission events: {}\n", total_intervalley_emission_events);
    fmt::print("Total self-scattering events: {}\n", total_self_scattering_events);

    m_observables.electric_field_V_per_m = m_cfg.m_electric_field.norm();

    if (m_observables.accumulated_time_s == 0.0) {
        throw std::runtime_error("no steady-state samples accumulated");
    }

    const double avg_vx_m_per_s = m_observables.weighted_velocity_x_m2_per_s2 / m_observables.accumulated_time_s;

    const double avg_energy_eV = m_observables.weighted_kinetic_energy_eV_s / m_observables.accumulated_time_s;

    fmt::print("Steady-state average vx: {:.6e} m/s\n", avg_vx_m_per_s);
    fmt::print("Steady-state average energy: {:.6f} eV\n", avg_energy_eV);
}

void bulk_amc_simulation::run_self_scattering_emc() {
    if (m_particles.empty()) {
        throw std::runtime_error("simulation not initialized");
    }

    if (m_cfg.m_final_time <= 0.0) {
        throw std::invalid_argument("final time must be > 0");
    }

    if (m_cfg.m_warmup_fraction < 0.0 || m_cfg.m_warmup_fraction >= 1.0) {
        throw std::invalid_argument("warmup fraction must be in [0, 1)");
    }

    if (m_gamma_max_s_1 <= 0.0) {
        throw std::invalid_argument("max self-scattering rate must be > 0");
    }

    m_observables = {};

    const double warmup_time             = m_cfg.m_warmup_fraction * m_cfg.m_final_time;
    double       max_observed_total_rate = 0.0;

    for (auto& p : m_particles) {
        while (p.state().time < m_cfg.m_final_time) {
            const double tau            = sample_free_flight_time();
            const double remaining_time = m_cfg.m_final_time - p.state().time;
            const double drift_time     = std::min(tau, remaining_time);

            drift_particle(p, drift_time);

            if (m_cfg.m_record_history) {
                p.record_state();
            }

            if (p.state().time >= warmup_time) {
                accumulate_particle_observables(p, drift_time);
            }

            if (drift_time < tau) {
                break;
            }

            const double total_rate = total_scattering_rate(p);
            max_observed_total_rate = std::max(max_observed_total_rate, total_rate);

            if (total_rate > m_gamma_max_s_1) {
                throw std::runtime_error("total scattering rate exceeds max self-scattering rate");
            }

            std::uniform_real_distribution<double> unif01(0.0, 0.0 + 1.0);
            const double                           u = unif01(m_rng);

            if (u < total_rate / m_gamma_max_s_1) {
                const auto channel = select_scattering_channel(p);
                apply_scattering_channel(p, channel);
            } else {
                p.increment_scattering_event_count();
                p.add_scattering_event(scattering_event::self_scattering);
            }

            if (m_cfg.m_record_history) {
                p.record_state();
            }
        }
    }

    fmt::print("Completed self-scattering EMC run with {} particles\n", m_particles.size());
    fmt::print("Maximum observed total scattering rate: {:.6e} s^-1\n", max_observed_total_rate);

    double  avg_energy = 0.0;
    vector3 avg_velocity{0.0, 0.0, 0.0};

    for (const auto& p : m_particles) {
        avg_energy += p.state().kinetic_energy;
        avg_velocity += p.state().velocity;
    }

    avg_energy /= static_cast<double>(m_particles.size());
    avg_velocity /= static_cast<double>(m_particles.size());

    fmt::print("Average kinetic energy: {:.6f} eV\n", avg_energy);
    fmt::print("Average velocity: ({:.6e}, {:.6e}, {:.6e}) m/s\n", avg_velocity.x(), avg_velocity.y(), avg_velocity.z());

    std::size_t total_acoustic_events               = 0;
    std::size_t total_intervalley_absorption_events = 0;
    std::size_t total_intervalley_emission_events   = 0;
    std::size_t total_self_scattering_events        = 0;

    for (const auto& p : m_particles) {
        const auto& events = p.history().scattering_events();

        total_acoustic_events += events[static_cast<std::size_t>(scattering_event::acoustic)];
        total_intervalley_absorption_events += events[static_cast<std::size_t>(scattering_event::intervalley_absorption)];
        total_intervalley_emission_events += events[static_cast<std::size_t>(scattering_event::intervalley_emission)];
        total_self_scattering_events += events[static_cast<std::size_t>(scattering_event::self_scattering)];
    }

    fmt::print("Total acoustic events: {}\n", total_acoustic_events);
    fmt::print("Total intervalley absorption events: {}\n", total_intervalley_absorption_events);
    fmt::print("Total intervalley emission events: {}\n", total_intervalley_emission_events);
    fmt::print("Total self-scattering events: {}\n", total_self_scattering_events);

    m_observables.electric_field_V_per_m = m_cfg.m_electric_field.norm();

    if (m_observables.accumulated_time_s == 0.0) {
        throw std::runtime_error("no steady-state samples accumulated");
    }

    const double avg_vx_m_per_s = m_observables.weighted_velocity_x_m2_per_s2 / m_observables.accumulated_time_s;

    const double avg_energy_eV = m_observables.weighted_kinetic_energy_eV_s / m_observables.accumulated_time_s;

    fmt::print("Steady-state average vx: {:.6e} m/s\n", avg_vx_m_per_s);
    fmt::print("Steady-state average energy: {:.6f} eV\n", avg_energy_eV);
}

/**
 * @brief Append the current observables to a CSV file. If the file does not exist, it will be created with a header. If it already exists,
 * a new line will be appended with the current observables values.
 *
 * @param filename
 */
void bulk_amc_simulation::export_observables_to_csv(const std::string& filename) const {
    std::ofstream file(filename, std::ios::app);
    if (!file.is_open()) {
        fmt::print(stderr, "Failed to open file for writing: {}\n", filename);
        return;
    }
    // Check if the file is empty to write the header
    if (file.tellp() == 0) {
        file << "electric_field_V_per_m,mean_velocity_x_m_per_s,mean_kinetic_energy_eV,sample_count\n";
    }
    file << fmt::format("{},{},{},{}\n",
                        m_observables.electric_field_V_per_m,
                        m_observables.weighted_velocity_x_m2_per_s2 / m_observables.accumulated_time_s,
                        m_observables.weighted_kinetic_energy_eV_s / m_observables.accumulated_time_s,
                        m_observables.accumulated_time_s);
}

void bulk_amc_simulation::export_particles_history_to_csv(const std::string& prefix_name) const {
    // Implementation for exporting particle history to CSV
    for (const auto& p : m_particles) {
        const auto&       history  = p.history();
        const std::string filename = fmt::format("{}_particle_{}.csv", prefix_name, p.index());
        std::ofstream     file(filename);
        if (!file.is_open()) {
            fmt::print(stderr, "Failed to open file for writing: {}\n", filename);
            continue;
        }

        // Write CSV header
        file << "time,position_x,position_y,position_z,local_k_x,local_k_y,local_k_z,velocity_x,velocity_y,velocity_z,kinetic_energy,gamma,"
                "valley_index\n";

        // Write particle history
        for (const auto& snapshot : history.snapshots()) {
            file << fmt::format("{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
                                snapshot.time,
                                snapshot.position.x(),
                                snapshot.position.y(),
                                snapshot.position.z(),
                                snapshot.local_k.x(),
                                snapshot.local_k.y(),
                                snapshot.local_k.z(),
                                snapshot.velocity.x(),
                                snapshot.velocity.y(),
                                snapshot.velocity.z(),
                                snapshot.kinetic_energy,
                                snapshot.gamma,
                                snapshot.valley_index);
        }
    }
}

}  // namespace uepm::amc