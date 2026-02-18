/**
 * @file single_part_fbmc.cpp
 * @author
 * @brief Single-particle FBMC with null-collision (Option A)
 * @version 0.1
 * @date 2025-09-19
 */

#include "single_part_fbmc.hpp"

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <vector>

#include "keldysh_impactio.hpp"
#include "physical_constants.hpp"

namespace uepm::fbmc {

Single_particle_simulation::Single_particle_simulation(uepm::mesh_bz::ElectronPhonon* ptr_mesh_bz,
                                                       const Bulk_environment&        bulk_env,
                                                       const Simulation_parameters&   sim_params,
                                                       std::size_t                    nb_particles)
    : m_ptr_mesh_bz(ptr_mesh_bz),
      m_bulk_env(bulk_env),
      m_sim_params(sim_params),
      m_nb_particles(nb_particles),
      m_time(0.0) {
    m_list_particle.reserve(m_nb_particles);
    for (std::size_t i = 0; i < m_nb_particles; ++i) {
        particle particle(i, particle_type::electron, m_ptr_mesh_bz);
        particle.set_position({0.0, 0.0, 0.0});
        m_list_particle.emplace_back(std::move(particle));
    }

    const double thermal_energy = (3.0 / 2.0) * m_bulk_env.m_temperature * uepm::constants::k_b_eV;
    std::cout << "Thermal energy at " << m_bulk_env.m_temperature << " K: " << thermal_energy << " eV\n";

    for (auto& particle : m_list_particle) {
        particle.set_energy(thermal_energy);
        particle.set_velocity({0.0, 0.0, 0.0});
    }

    // Draw an initial k-vector at that energy (band 0 here)
    fmt::print("Drawing initial k-vector at thermal energy...\n");
#pragma omp parallel for
    for (auto& particle : m_list_particle) {
        vector3 initial_k = m_ptr_mesh_bz->draw_random_k_point_at_energy(thermal_energy, 0, particle.get_random_generator());
        particle.set_k_vector(initial_k);
        std::cout << "Initial k-vector (drawn at thermal energy): " << initial_k << "\n";
    }

    // Attach containing tetra
    for (auto& particle : m_list_particle) {
        uepm::mesh_bz::Tetra* containing_tetra = m_ptr_mesh_bz->find_tetra_at_location(particle.get_k_vector());
        if (containing_tetra == nullptr) {
            throw std::runtime_error("Initial k-point is out of the Brillouin zone mesh.");
        }
        particle.set_containing_bz_mesh_tetra(containing_tetra);
        // Physical total rate at current k
        std::array<double, 8> rates = particle.interpolate_phonon_scattering_rate_at_location(particle.get_k_vector());
        const double          Gamma = std::accumulate(rates.begin(), rates.end(), 0.0);
        particle.set_gamma(Gamma);

        // Sync velocity & energy
        particle.update_group_velocity();
        particle.update_energy();

        const double init_energy_true = containing_tetra->interpolate_energy_at_band(particle.get_k_vector(), 0);

        fmt::print("Initial energy (set):  {:.10f} eV\n", thermal_energy);
        fmt::print("Initial energy (true): {:.10f} eV\n", init_energy_true);
    }
}


// void Single_particle_simulation::run_simulation() {
//     // Upper bound for null-collision
//     double p_gamma_elph = m_ptr_mesh_bz->compute_P_Gamma() * 1.2;
//     if (!(p_gamma_elph > 0.0) || !std::isfinite(p_gamma_elph)) {
//         fmt::print("Invalid p_gamma (<=0 or NaN)\n");
//         return;
//     }
//     KeldyshImpactIonization keldysh_impactio;
//     double                  max_energy     = 8.0;  // Example value, replace with actual max energy
//     double                  p_gamma_impact = keldysh_impactio.compute_rate(max_energy);

//     double p_gamma = p_gamma_elph + p_gamma_impact;
//     if (!(p_gamma > 0.0) || !std::isfinite(p_gamma)) {
//         fmt::print("Invalid total p_gamma (<=0 or NaN)\n");
//         return;
//     }
//     fmt::print("Using p_gamma = {:.3e} 1/s (el-ph: {:.3e}, impact ionization: {:.3e})\n", p_gamma, p_gamma_elph, p_gamma_impact);

//     const double WarningGammaThreshold = 1e10;  // 1/s
//     const double T_end                 = m_sim_params.m_simulation_time;

// // Parallelize over particles; each runs to T_end; each particle are independent (it's like multiple single-particle simulations)
// #pragma omp parallel for schedule(dynamic) num_threads(m_sim_params.m_nb_openmp_threads)
//     for (std::size_t idx = 0; idx < m_list_particle.size(); ++idx) {
//         fmt::print("Running simulation for particle {}\n", idx);
//         auto& particle = m_list_particle[idx];

//         std::uniform_real_distribution<double> U01(0.0, 1.0);

//         while (particle.get_time() <= T_end) {
//             if (particle.get_index() == 0 && particle.get_iter() % 1000 == 0) {
//                 fmt::print("Particle {} at time {:.3e} / {:.3e} s, iteration {}, energy {:.4f} eV\n",
//                            particle.get_index(),
//                            particle.get_time(),
//                            T_end,
//                            particle.get_iter(),
//                            particle.get_energy());
//             }

//             // Draw free flight with upper bound p_gamma
//             particle.draw_free_flight_time(p_gamma);
//             const double dt = particle.get_current_free_flight_time();
//             // Advance particle time (add a setter if needed)
//             // particle.set_time(particle.get_time() + dt);

//             // Drift
//             particle.update_k_vector(m_bulk_env.m_electric_field);

//             // Re-attach tetra; fold if needed
//             uepm::mesh_bz::Tetra* containing_tetra = m_ptr_mesh_bz->find_tetra_at_location(particle.get_k_vector());
//             if (containing_tetra == nullptr) {
//                 const vector3 folded_k = m_ptr_mesh_bz->retrieve_k_inside_mesh_geometry(particle.get_k_vector());
//                 particle.set_k_vector(folded_k);
//                 containing_tetra = m_ptr_mesh_bz->find_tetra_at_location(particle.get_k_vector());
//                 if (containing_tetra == nullptr) {
//                     throw std::runtime_error("Particle k-point is out of the Brillouin zone mesh after folding.");
//                 }
//             }
//             particle.set_containing_bz_mesh_tetra(containing_tetra);

//             // Update E & v after drift
//             particle.update_energy();
//             particle.update_group_velocity();

//             // Physical total rate at current k
//             std::array<double, 8> rates_elph = particle.interpolate_phonon_scattering_rate_at_location(particle.get_k_vector());
//             std::array<double, 8> rates_tetra_elph =
//                 containing_tetra->interpolate_phonon_scattering_rate_at_location(particle.get_k_vector(), particle.get_band_index());
//             const double Gamma_elph           = std::accumulate(rates_elph.begin(), rates_elph.end(), 0.0) + p_gamma_impact;
//             double       sum_rates_tetra = std::accumulate(rates_tetra_elph.begin(), rates_tetra_elph.end(), 0.0);
//             double       energy_tetra    = containing_tetra->interpolate_energy_at_band(particle.get_k_vector(), particle.get_band_index());
//             double       energy_particle = particle.get_energy();
//             double       rate_impactio   = keldysh_impactio.compute_rate(energy_particle);
//             double       Gamma           = Gamma_elph + rate_impactio;
//             // Debug checks
//             if (std::fabs(energy_tetra - particle.get_energy()) > 1e-6) {
//                 fmt::print("Warning: Energy mismatch: particle = {:.6e} eV, tetra = {:.6e} eV\n", particle.get_energy(), energy_tetra);
//             }
//             if (std::fabs(Gamma - sum_rates_tetra) > 1e-6) {
//                 fmt::print("Warning: Gamma mismatch: Gamma = {:.6e}, sum_rates_tetra = {:.6e}\n", Gamma, sum_rates_tetra);
//             }
//             particle.set_gamma(Gamma);

//             if (!(Gamma > 0.0) || !std::isfinite(Gamma)) {
//                 // Null event only this step
//                 continue;
//             }
//             if (Gamma < WarningGammaThreshold) {
//                 fmt::print("Warning: Gamma < {} for particle {} at time {} s\n",
//                            WarningGammaThreshold,
//                            particle.get_index(),
//                            particle.get_time());
//             }

//             // Null-collision acceptance
//             double accept = Gamma / p_gamma;
//             if (accept > 1.0) {
//                 accept = 1.0;  // clamp if bound violated
//             }

//             if (U01(particle.get_random_generator()) > accept) {
//                 particle.add_scattering_event_to_history(9);
//                 continue;
//             }
//             // to gain place, we update the history only after a real event
//             particle.update_history();

//             // Real event: pick channel
//             const double rsel = U01(particle.get_random_generator()) * Gamma;
//             double       cum  = 0.0;
//             std::size_t  event_idx   = 0;
//             for (; event_idx < rates_elph.size(); ++event_idx) {
//                 cum += rates_elph[event_idx];
//                 if (rsel <= cum) {
//                     break;
//                 }
//             }
//             if (event_idx == rates_elph.size()) {
//                 // Impact ionization event
//                 event_idx = 8;  // index 8 for impact ionization
//                 particle.select_final_state_after_impact_ionization();
//                 fmt::print("Particle {} underwent impact ionization at time {:.3e} s, energy {:.3f} eV\n",
//                            particle.get_index(),
//                            particle.get_time(),
//                            particle.get_energy());
//             } else {
//                 // Phonon scattering event
//                 // event_idx is already set correctly
//                 particle.select_final_state_after_phonon_scattering(event_idx);
//             }
//             particle.add_scattering_event_to_history(event_idx);




//         }  // while particle
//         particle.update_history();  // final update
//                                     // Print summary for this particle
//         particle.print_history_summary();

//     }  // omp parallel for
// }

void Single_particle_simulation::run_simulation() {
    // Upper bound for null-collision (el-ph)
    double p_gamma_elph = m_ptr_mesh_bz->compute_P_Gamma() * 1.2;
    if (!(p_gamma_elph > 0.0) || !std::isfinite(p_gamma_elph)) {
        fmt::print("Invalid p_gamma_elph (<=0 or NaN)\n");
        return;
    }

    // Upper bound contribution for impact ionization (global, based on an assumed max energy)
    KeldyshImpactIonization keldysh_impactio;
    double                  max_energy     = 8.0;  // must be a true upper bound of particle energy
    double                  p_gamma_impact = keldysh_impactio.compute_rate(max_energy) * 1.2;
    if (!(p_gamma_impact >= 0.0) || !std::isfinite(p_gamma_impact)) {
        fmt::print("Invalid p_gamma_impact (NaN)\n");
        return;
    }

    const double p_gamma = p_gamma_elph + p_gamma_impact;
    if (!(p_gamma > 0.0) || !std::isfinite(p_gamma)) {
        fmt::print("Invalid total p_gamma (<=0 or NaN)\n");
        return;
    }

    fmt::print("Using p_gamma = {:.3e} 1/s (el-ph: {:.3e}, impact bound: {:.3e})\n", p_gamma, p_gamma_elph, p_gamma_impact);

    const double WarningGammaThreshold = 1e10;  // 1/s
    const double T_end                 = m_sim_params.m_simulation_time;

#pragma omp parallel for schedule(dynamic) num_threads(m_sim_params.m_nb_openmp_threads)
    for (std::size_t idx = 0; idx < m_list_particle.size(); ++idx) {
#pragma omp critical
        fmt::print("Running simulation for particle {}\n", idx);

        auto&                                  particle = m_list_particle[idx];
        std::uniform_real_distribution<double> U01(0.0, 1.0);

        while (particle.get_time() <= T_end) {
            if (particle.get_index() == 0 && particle.get_iter() % 1000 == 0) {
#pragma omp critical
                fmt::print("Particle {} at time {:.3e} / {:.3e} s, iteration {}, energy {:.4f} eV\n",
                           particle.get_index(),
                           particle.get_time(),
                           T_end,
                           particle.get_iter(),
                           particle.get_energy());
            }

            // Draw free flight with global upper bound p_gamma
            particle.draw_free_flight_time(p_gamma);
            const double dt = particle.get_current_free_flight_time();

            // Advance particle time (required if draw_free_flight_time does not do it internally)
            particle.set_time(particle.get_time() + dt);

            // Drift
            particle.update_k_vector(m_bulk_env.m_electric_field);

            // Re-attach tetra; fold if needed
            uepm::mesh_bz::Tetra* containing_tetra = m_ptr_mesh_bz->find_tetra_at_location(particle.get_k_vector());
            if (containing_tetra == nullptr) {
                const vector3 folded_k = m_ptr_mesh_bz->retrieve_k_inside_mesh_geometry(particle.get_k_vector());
                particle.set_k_vector(folded_k);
                containing_tetra = m_ptr_mesh_bz->find_tetra_at_location(particle.get_k_vector());
                if (containing_tetra == nullptr) {
                    throw std::runtime_error("Particle k-point is out of the Brillouin zone mesh after folding.");
                }
            }
            particle.set_containing_bz_mesh_tetra(containing_tetra);

            // Update E & v after drift
            particle.update_energy();
            particle.update_group_velocity();

            // Physical rates at current k
            std::array<double, 8> rates_elph = particle.interpolate_phonon_scattering_rate_at_location(particle.get_k_vector());
            const double          Gamma_elph = std::accumulate(rates_elph.begin(), rates_elph.end(), 0.0);

            const double energy_particle = particle.get_energy();
            const double rate_impactio   = keldysh_impactio.compute_rate(energy_particle);
            const double Gamma           = Gamma_elph + rate_impactio;

            // Debug checks (optional)
            const double energy_tetra = containing_tetra->interpolate_energy_at_band(particle.get_k_vector(), particle.get_band_index());
            if (std::fabs(energy_tetra - energy_particle) > 1e-6) {
#pragma omp critical
                fmt::print("Warning: Energy mismatch: particle = {:.6e} eV, tetra = {:.6e} eV\n", energy_particle, energy_tetra);
            }

            particle.set_gamma(Gamma);

            if (!(Gamma > 0.0) || !std::isfinite(Gamma)) {
                // Null event only this step
                continue;
            }

            if (Gamma < WarningGammaThreshold) {
#pragma omp critical
                fmt::print("Warning: Gamma < {} for particle {} at time {} s\n",
                           WarningGammaThreshold,
                           particle.get_index(),
                           particle.get_time());
            }

            // Null-collision acceptance
            double accept = Gamma / p_gamma;

            // Bound violation check: null-collision method requires Gamma <= p_gamma
            // Allow a tiny epsilon for floating rounding.
            constexpr double eps = 1e-12;
            if (accept > 1.0 + eps) {
                throw std::runtime_error("Null-collision bound violated: Gamma > p_gamma. Increase bounds.");
            }
            if (accept > 1.0) {
                accept = 1.0;
            }

            if (U01(particle.get_random_generator()) > accept) {
                particle.add_scattering_event_to_history(9);  // null/self-scatter
                continue;
            }

            // Real event
            particle.update_history();

            const double rsel = U01(particle.get_random_generator()) * Gamma;

            // Select phonon channel or impact ionization
            const double sum_elph = Gamma_elph;

            if (rsel > sum_elph) {
                // Impact ionization event
                const std::size_t event_idx = 8;
                particle.select_final_state_after_impact_ionization();
#pragma omp critical
                fmt::print("Particle {} underwent impact ionization at time {:.3e} s, energy {:.3f} eV\n",
                           particle.get_index(),
                           particle.get_time(),
                           particle.get_energy());
                particle.add_scattering_event_to_history(event_idx);
            } else {
                // Phonon scattering event
                double      cum       = 0.0;
                std::size_t event_idx = 0;
                for (; event_idx < rates_elph.size(); ++event_idx) {
                    cum += rates_elph[event_idx];
                    if (rsel <= cum) {
                        break;
                    }
                }
                if (event_idx >= rates_elph.size()) {
                    // Guard against numerical edge cases (rsel very close to sum_elph)
                    event_idx = rates_elph.size() - 1;
                }

                particle.select_final_state_after_phonon_scattering(event_idx);
                particle.add_scattering_event_to_history(event_idx);
            }
        }  // while particle

        particle.update_history();  // final update
        // particle.print_history_summary();
    }  // omp parallel for
}

void Single_particle_simulation::extract_stats_and_export(const std::string& filename) {
    double mean_energy = 0.0;
    double mean_velocity_norm = 0.0;
    double ionization_coeff = 0.0;
    for (const auto& particle : m_list_particle) {
        mean_energy += particle.compute_mean_energy();
        mean_velocity_norm += particle.extract_global_average_velocity().norm();
        ionization_coeff += particle.extract_impact_ionization_coeff();
    }
    mean_energy /= static_cast<double>(m_list_particle.size());
    mean_velocity_norm /= static_cast<double>(m_list_particle.size());
    ionization_coeff /= static_cast<double>(m_list_particle.size());

    fmt::print("Extracted stats across {} particles:\n", m_list_particle.size());
    fmt::print("  Mean energy: {:.6f} eV\n", mean_energy);
    fmt::print("  Mean velocity norm: {:.6e} m/s\n", mean_velocity_norm);
    fmt::print("  Mean impact ionization coefficient: {:.3e} 1/s\n", ionization_coeff);

    // Export to CSV
    std::ofstream ofs(filename);
    if (!ofs) {
        fmt::print("Error: Could not open file {} for writing\n", filename);
        return;
    }
    ofs << "mean_energy_eV,mean_velocity_norm_m_per_s,mean_ionization_coeff_1_per_s\n";
    ofs << fmt::format("{:.6f},{:.6e},{:.3e}\n", mean_energy, mean_velocity_norm, ionization_coeff);
    fmt::print("Exported extracted stats to {}\n", filename);
}

void Single_particle_simulation::export_history(const std::string& filename) {
    for (auto& particle : m_list_particle) {
        particle.print_history_summary();
        std::string filename_particle = fmt::format("{}_particle_{}.csv", filename, particle.get_index());
        particle.export_history_to_csv(filename_particle);
        fmt::print("Exported history of particle {} to {}\n", particle.get_index(), filename_particle);
    }
}

}  // namespace uepm::fbmc
