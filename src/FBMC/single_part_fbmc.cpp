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
    // Initialize the particles
    m_list_particle.reserve(m_nb_particles);
    for (std::size_t i = 0; i < m_nb_particles; ++i) {
        particle particle(i, particle_type::electron, m_ptr_mesh_bz);
        particle.set_position({0.0, 0.0, 0.0});
        m_list_particle.emplace_back(std::move(particle));
    }

    const double thermal_energy = (3.0 / 2.0) * m_bulk_env.m_temperature * uepm::constants::k_b_eV;
    std::cout << "Thermal energy at " << m_bulk_env.m_temperature << " K: " << thermal_energy << " eV\n";

    // for (std::size_t i = 0; i < m_nb_particles; ++i) {
    //     m_particle[i].set_energy(thermal_energy);
    //     m_particle[i].set_velocity({0.0, 0.0, 0.0});
    // }
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

        // Sync velocity & energy
        particle.update_group_velocity();
        particle.update_energy();

        const double init_energy_true = containing_tetra->interpolate_energy_at_band(particle.get_k_vector(), 0);

        fmt::print("Initial energy (set):  {:.10f} eV\n", thermal_energy);
        fmt::print("Initial energy (true): {:.10f} eV\n", init_energy_true);
    }
}
void Single_particle_simulation::run_simulation() {
    // Upper bound for null-collision
    double p_gamma = m_ptr_mesh_bz->compute_P_Gamma() * 1.2;
    if (!(p_gamma > 0.0) || !std::isfinite(p_gamma)) {
        fmt::print("Invalid p_gamma (<=0 or NaN)\n");
        return;
    }
    for (auto& p : m_list_particle) {
        p.set_gamma(p_gamma);
    }

    // Global tallies (optional)
    int    nb_foldings_total = 0;
    double max_time_reached  = 0.0;

    const double T_end = m_sim_params.m_simulation_time;

// Parallelize over particles; each runs to T_end
#pragma omp parallel for reduction(+ : nb_foldings_total) reduction(max : max_time_reached) schedule(dynamic)
    for (std::size_t idx = 0; idx < m_list_particle.size(); ++idx) {
        fmt::print("Running simulation for particle {}\n", idx);
        auto& particle = m_list_particle[idx];

        std::uniform_real_distribution<double> U01(0.0, 1.0);
        int                                    nb_foldings_local = 0;

        while (particle.get_time() < T_end) {
            // Save state
            particle.update_history();

            // Draw free flight with upper bound p_gamma
            particle.draw_free_flight_time(p_gamma);
            const double dt = particle.get_current_free_flight_time();
            // Advance particle time (add a setter if needed)
            // particle.set_time(particle.get_time() + dt);

            // Drift
            particle.update_k_vector(m_bulk_env.m_electric_field);

            // Re-attach tetra; fold if needed
            uepm::mesh_bz::Tetra* containing_tetra = m_ptr_mesh_bz->find_tetra_at_location(particle.get_k_vector());
            if (containing_tetra == nullptr) {
                const vector3 folded_k = m_ptr_mesh_bz->retrieve_k_inside_mesh_geometry(particle.get_k_vector());
                particle.set_k_vector(folded_k);
                containing_tetra = m_ptr_mesh_bz->find_tetra_at_location(particle.get_k_vector());
                if (containing_tetra == nullptr) {
                    // Mark and stop this particle; keep others running
                    // (optional) set a flag on the particle
                    break;
                } else {
                    ++nb_foldings_local;
                }
            }
            particle.set_containing_bz_mesh_tetra(containing_tetra);

            // Update E & v after drift
            particle.update_energy();
            particle.update_group_velocity();

            // Physical total rate at current k
            std::array<double, 8> rates = particle.interpolate_phonon_scattering_rate_at_location(particle.get_k_vector());
            const double          Gamma = std::accumulate(rates.begin(), rates.end(), 0.0);
            particle.set_gamma(Gamma);

            if (!(Gamma > 0.0) || !std::isfinite(Gamma)) {
                // Null event only this step
                continue;
            }

            // Null-collision acceptance
            double accept = Gamma / p_gamma;
            if (accept > 1.0) {
                accept = 1.0;  // clamp if bound violated
            }

            if (U01(particle.get_random_generator()) > accept) {
                // self-scatter: do nothing
                continue;
            }

            // Real event: pick channel
            const double rsel = U01(particle.get_random_generator()) * Gamma;
            double       cum  = 0.0;
            std::size_t  ev   = 0;
            for (; ev < rates.size(); ++ev) {
                cum += rates[ev];
                if (rsel <= cum) {
                    break;
                }
            }
            if (ev >= rates.size()) {
                ev = rates.size() - 1;
            }

            particle.select_final_state_after_phonon_scattering(ev);
        }  // while particle

        nb_foldings_total += nb_foldings_local;
        if (particle.get_time() > max_time_reached) {
            max_time_reached = particle.get_time();
        }
    }  // omp parallel for

    std::cout << "Simulation finished. max_time = " << max_time_reached << " s, total foldings = " << nb_foldings_total << "\n";
}

void Single_particle_simulation::export_history(const std::string& filename) {
    for (auto& particle : m_list_particle) {
        std::string   filename_particle = fmt::format("{}_particle_{}.csv", filename, particle.get_index());
        std::ofstream history_file(filename_particle);
        if (!history_file.is_open()) {
            std::cerr << "Could not open history file for writing.\n";
            return;
        }
        // Prefer your constants over M_PI for portability
        constexpr double lattice_constant     = 5.43e-10;  // meters (Si)
        const double     normalization_factor = 2.0 * uepm::constants::pi / lattice_constant;

        const auto history = particle.get_history();

        history_file << "time,gamma,x,y,z,kx,ky,kz,vx,vy,vz,energy\n";
        for (std::size_t step = 0; step < history.get_number_of_steps(); ++step) {
            history_file << history.m_time_history[step] << "," << history.m_gammas[step] << "," << history.m_positions[step].x() << ","
                         << history.m_positions[step].y() << "," << history.m_positions[step].z() << ","
                         << history.m_k_vectors[step].x() / normalization_factor << ","
                         << history.m_k_vectors[step].y() / normalization_factor << ","
                         << history.m_k_vectors[step].z() / normalization_factor << "," << history.m_velocities[step].x() << ","
                         << history.m_velocities[step].y() << "," << history.m_velocities[step].z() << "," << history.m_energies[step]
                         << "\n";
        }
        history_file.close();
        std::cout << "Particle history exported to " << filename << "\n";
    }
}

}  // namespace uepm::fbmc
