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
void bulk_amc_simulation::initialize() {
    m_valleys              = make_silicon_delta_valleys();
    m_intervalley_branches = make_silicon_intervalley_phonon_branches();

    // std::cout << "Initialized " << m_valleys.size() << " valleys\n";
    // // for (std::size_t i = 0; i < m_valleys.size(); ++i) {
    // //     std::cout << "  valley[" << i << "] = " << m_valleys[i].name() << '\n';
    // // }

    // // std::cout << "Initialized " << m_intervalley_branches.size() << " intervalley phonon branches\n";
    // // for (const auto& branch : m_intervalley_branches) {
    // //     std::cout << "  branch " << branch.m_name << " family=" << (branch.m_family == intervalley_family::f ? "f" : "g")
    // //               << " final_valleys=" << branch.m_final_valley_count << '\n';
    // // }



    m_particles.clear();
    m_particles.reserve(m_cfg.m_number_of_particles);

    for (std::size_t i = 0; i < m_cfg.m_number_of_particles; ++i) {
        particle_amc p{i, particle_type::electron, 1.0};

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

struct sampled_intervalley_event {
    const intervalley_phonon_branch* branch             = nullptr;
    bool                             absorption         = false;
    double                           rate               = 0.0;
    std::size_t                      destination_valley = 0;
    double                           final_energy_eV    = 0.0;
};

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

std::size_t draw_random_different_valley(std::size_t current_valley, std::size_t number_of_valleys, std::mt19937_64& rng) {
    if (number_of_valleys < 2) {
        throw std::invalid_argument("need at least two valleys for intervalley scattering");
    }
    if (current_valley >= number_of_valleys) {
        throw std::out_of_range("current valley index out of range");
    }

    std::uniform_int_distribution<std::size_t> dist(0, number_of_valleys - 2);
    std::size_t                                draw = dist(rng);

    if (draw >= current_valley) {
        ++draw;
    }
    return draw;
}

void bulk_amc_simulation::scatter_particle(particle_amc& p, double dt) {
    if (dt < 0.0) {
        throw std::invalid_argument("scatter time step must be non-negative");
    }

    const auto current_valley_index = p.state().valley_index;
    if (current_valley_index >= m_valleys.size()) {
        throw std::out_of_range("invalid valley index in scatter_particle");
    }

    const double energy_eV = p.state().kinetic_energy;

    struct sampled_intervalley_event {
        const intervalley_phonon_branch* branch             = nullptr;
        bool                             absorption         = false;
        double                           rate               = 0.0;
        std::size_t                      destination_valley = 0;
        double                           final_energy_eV    = 0.0;
    };

    const auto& current_valley = m_valleys[current_valley_index];
    const double rate_ac = acoustic_scattering_rate_silicon(current_valley, energy_eV, m_cfg.m_lattice_temperature);


    std::vector<sampled_intervalley_event> intervalley_events;
    intervalley_events.reserve(2 * m_intervalley_branches.size());

    double total_intervalley_rate = 0.0;

    for (const auto& branch : m_intervalley_branches) {
        {
            const double rate_abs = intervalley_scattering_rate(current_valley, branch, energy_eV, true, m_cfg.m_lattice_temperature);

            if (rate_abs > 0.0) {
                intervalley_events.push_back(sampled_intervalley_event{
                    .branch             = &branch,
                    .absorption         = true,
                    .rate               = rate_abs,
                    .destination_valley = draw_intervalley_destination_valley(branch, current_valley_index, m_rng),
                    .final_energy_eV    = energy_eV + branch.m_phonon_energy_eV});
                total_intervalley_rate += rate_abs;
            }
        }

        {
            const double rate_em = intervalley_scattering_rate(current_valley, branch, energy_eV, false, m_cfg.m_lattice_temperature);

            if (rate_em > 0.0) {
                const double final_energy_eV = energy_eV - branch.m_phonon_energy_eV;
                if (final_energy_eV >= 0.0) {
                    intervalley_events.push_back(sampled_intervalley_event{
                        .branch             = &branch,
                        .absorption         = false,
                        .rate               = rate_em,
                        .destination_valley = draw_intervalley_destination_valley(branch, current_valley_index, m_rng),
                        .final_energy_eV    = final_energy_eV});
                    total_intervalley_rate += rate_em;
                }
            }
        }
    }

    const double total_rate = rate_ac + total_intervalley_rate;
    if (total_rate <= 0.0) {
        return;
    }

    const double scatter_probability = 1.0 - std::exp(-total_rate * dt);

    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    if (unif01(m_rng) >= scatter_probability) {
        return;
    }

    const double r_select = unif01(m_rng) * total_rate;

    if (r_select < rate_ac) {
        const auto& valley = m_valleys[current_valley_index];

        p.state().local_k        = valley.draw_random_k_valley_at_energy(energy_eV, m_rng);
        p.state().gamma          = valley.gamma_from_k_valley(p.state().local_k);
        p.state().kinetic_energy = valley.kinetic_energy_from_gamma(p.state().gamma);
        p.state().velocity       = valley.velocity_from_k_valley(p.state().local_k);

        p.increment_total_step();
        p.add_scattering_event(scattering_event::acoustic);
        return;
    }

    double cumulative = rate_ac;
    for (const auto& event : intervalley_events) {
        cumulative += event.rate;
        if (r_select < cumulative) {
            const auto& dst_valley = m_valleys[event.destination_valley];

            p.state().valley_index   = event.destination_valley;
            p.state().local_k        = dst_valley.draw_random_k_valley_at_energy(event.final_energy_eV, m_rng);
            p.state().gamma          = dst_valley.gamma_from_k_valley(p.state().local_k);
            p.state().kinetic_energy = dst_valley.kinetic_energy_from_gamma(p.state().gamma);
            p.state().velocity       = dst_valley.velocity_from_k_valley(p.state().local_k);

            p.increment_total_step();
            p.add_scattering_event(scattering_event::intervalley);
            return;
        }
    }

    throw std::runtime_error("failed to select a scattering event");
}

// double bulk_amc_simulation::acoustic_scattering_rate_silicon(double energy_eV, double temperature_K) const {
//     const double energy_clamped = std::max(energy_eV, 1.0e-6);

//     // Simple fitted placeholder coefficient for now.
//     // Units chosen so the result is in s^-1 when E is in eV and T in K.
//     constexpr double C_ac = 1.0e12;

//     return C_ac * (temperature_K / 300.0) * std::sqrt(energy_clamped);
// }

void bulk_amc_simulation::run() {
    if (m_particles.empty()) {
        throw std::runtime_error("simulation not initialized");
    }

    constexpr double  dt      = 5.0e-15;
    const std::size_t n_steps = static_cast<std::size_t>(std::ceil(m_cfg.m_final_time / dt));

    for (std::size_t step = 0; step < n_steps; ++step) {
        for (auto& p : m_particles) {
            drift_particle(p, dt);
            scatter_particle(p, dt);

            if (m_cfg.m_record_history) {
                p.record_state();
            }
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


    std::size_t total_acoustic_events    = 0;
    std::size_t total_intervalley_events = 0;

    for (const auto& p : m_particles) {
        const auto& events = p.history().scattering_events();
        total_acoustic_events += events[static_cast<std::size_t>(scattering_event::acoustic)];
        total_intervalley_events += events[static_cast<std::size_t>(scattering_event::intervalley)];
    }

    fmt::print("Total acoustic events: {}\n", total_acoustic_events);
    fmt::print("Total intervalley events: {}\n", total_intervalley_events);
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