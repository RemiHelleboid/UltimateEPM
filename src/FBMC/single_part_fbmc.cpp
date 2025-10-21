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
                                                       const Simulation_parameters&   sim_params)
    : m_ptr_mesh_bz(ptr_mesh_bz),
      m_bulk_env(bulk_env),
      m_sim_params(sim_params),
      m_time(0.0)  // ✅ initialize the member (fixes the shadowing bug)
{
    // Initialize the particle
    const std::size_t index = 0;
    m_particle              = particle(index, particle_type::electron, m_ptr_mesh_bz);

    m_particle.set_position({0.0, 0.0, 0.0});

    // Quick-and-simple initial energy: k_B T (note: not an equilibrium draw)
    const double thermal_energy = (3.0 / 2.0) * m_bulk_env.m_temperature * uepm::constants::k_b_eV;
    m_particle.set_energy(thermal_energy);
    m_particle.set_velocity({0.0, 0.0, 0.0});

    std::cout << "Thermal energy at " << m_bulk_env.m_temperature << " K: " << thermal_energy << " eV\n";

    // Draw an initial k-vector at that energy (band 0 here)
    vector3 initial_k = m_ptr_mesh_bz->draw_random_k_point_at_energy(thermal_energy, 0, m_particle.get_random_generator());
    m_particle.set_k_vector(initial_k);

    std::cout << "Initial k-vector (drawn at thermal energy): " << initial_k << "\n";

    // Attach containing tetra
    uepm::mesh_bz::Tetra* containing_tetra = m_ptr_mesh_bz->find_tetra_at_location(m_particle.get_k_vector());
    if (containing_tetra == nullptr) {
        throw std::runtime_error("Initial k-point is out of the Brillouin zone mesh.");
    }
    m_particle.set_containing_bz_mesh_tetra(containing_tetra);

    // Sync velocity & energy
    m_particle.update_group_velocity();
    m_particle.update_energy();

    const double init_energy_true = containing_tetra->interpolate_energy_at_band(m_particle.get_k_vector(), 0);

    fmt::print("Initial energy (set):  {:.10f} eV\n", thermal_energy);
    fmt::print("Initial energy (true): {:.10f} eV\n", init_energy_true);
}

void Single_particle_simulation::run_simulation() {
    // Global upper bound p_gamma for null-collision
    double p_gamma = m_ptr_mesh_bz->compute_P_Gamma();
    // Get some margins
    p_gamma *= 1.2;
    fmt::print("Using p_gamma = {:.4e} 1/s for null-collision upper bound.\n", p_gamma);

    if (!(p_gamma > 0.0) || !std::isfinite(p_gamma)) {
        fmt::print("Invalid p_gamma (<=0 or NaN).\n");
        return;
    }

    m_particle.reserve_history(10000);
    int nb_foldings = 0;

    // Uniform [0,1] for selections
    std::uniform_real_distribution<double> U01(0.0, 1.0);

    while (m_time < m_sim_params.m_simulation_time) {

        // Save current state to history
        m_particle.update_history();

        // ---- Draw free-flight time with upper bound p_gamma (null-collision scheme) ----
        m_particle.draw_free_flight_time(p_gamma);
        const double dt = m_particle.get_current_free_flight_time();
        m_time += dt;

        m_particle.update_k_vector(m_bulk_env.m_electric_field);

        // ---- Keep k inside the mesh/BZ (fold if needed), reattach tetra ----
        uepm::mesh_bz::Tetra* containing_tetra = m_ptr_mesh_bz->find_tetra_at_location(m_particle.get_k_vector());
        if (containing_tetra == nullptr) {
            const bool is_inside_bz = m_ptr_mesh_bz->inside_ws_bcc(m_particle.get_k_vector());
            if (is_inside_bz) {
                std::cout << "\nWARNING: Particle inside WS cell but outside mesh.\n";
            } else {
                std::cout << "\nINFO: Particle outside WS cell and mesh.\n";
            }

            // Fold back into 1st BZ
            const vector3 folded_k = m_ptr_mesh_bz->retrieve_k_inside_mesh_geometry(m_particle.get_k_vector());
            m_particle.set_k_vector(folded_k);
            containing_tetra = m_ptr_mesh_bz->find_tetra_at_location(m_particle.get_k_vector());
            if (containing_tetra == nullptr) {
                std::cerr << "\nParticle still outside the Brillouin zone mesh after folding back.\n";
                break;
            } else {
                std::cout << "\nFolded back to k = " << folded_k.x() << " " << folded_k.y() << " " << folded_k.z() << "\n";
                ++nb_foldings;
            }
        }
        m_particle.set_containing_bz_mesh_tetra(containing_tetra);

        // ---- Update energy & velocity after drift ----
        m_particle.update_energy();
        m_particle.update_group_velocity();

        // ---- Compute physical total rate at current k ----
        std::array<double, 8> scattering_rates = m_particle.interpolate_phonon_scattering_rate_at_location(m_particle.get_k_vector());

        const double Gamma = std::accumulate(scattering_rates.begin(), scattering_rates.end(), 0.0);
        m_particle.set_gamma(Gamma);

        if (!(Gamma > 0.0) || !std::isfinite(Gamma)) {
            // No physical channels at this state → null collision only (do nothing)
            continue;
        }
        if (Gamma < 1e10) {
            fmt::print("\n[WARN] Very low total scattering rate Gamma = {:.4e} 1/s at energy {:.4e} eV in band {} at k = ({:.4e}, {:.4e}, {:.4e}).\n",
                       Gamma,
                       m_particle.get_energy(),
                       m_particle.get_band_index(),
                       m_particle.get_k_vector().x(),
                       m_particle.get_k_vector().y(),
                       m_particle.get_k_vector().z());
        }

        // ---- Null-collision acceptance: accept real event with prob = Gamma / p_gamma ----
        double accept = Gamma / p_gamma;
        if (accept > 1.0) {
            // Bound violated somewhere; clamp (optional: log under verbose)
            std::cout << "\n[WARN] Gamma > p_gamma (" << Gamma << " > " << p_gamma << "). Clamping accept=1.\n";

            accept = 1.0;
        }

        if (U01(m_particle.get_random_generator()) > accept) {
            // ---- Null event (self scattering): keep state, go to next loop ----
            continue;
        }

        // ---- Real event: pick channel with weights ~ gamma_i/Gamma ----
        const double rsel             = U01(m_particle.get_random_generator()) * Gamma;
        double       cum              = 0.0;
        std::size_t  scattering_event = 0;
        for (; scattering_event < scattering_rates.size(); ++scattering_event) {
            cum += scattering_rates[scattering_event];
            if (rsel <= cum) {
                break;
            }
        }
        if (scattering_event >= scattering_rates.size()) {
            // Numerical corner case: fallback to last channel
            scattering_event = scattering_rates.size() - 1;
        }

        // Your helper returns (band index, tetra index). We'll keep it.
        m_particle.select_final_state_after_phonon_scattering(scattering_event);

        fmt::print("\r time: {:.6e} / {:.6e} s - energy: {:.4e} eV", m_time, m_sim_params.m_simulation_time, m_particle.get_energy());
        std::cout.flush();
    }

    std::cout << "\nSimulation finished at time " << m_time << " s\n";
    std::cout << "Number of foldings: " << nb_foldings << "\n";
}

void Single_particle_simulation::export_history(const std::string& filename) {
    std::ofstream history_file(filename);
    if (!history_file.is_open()) {
        std::cerr << "Could not open history file for writing.\n";
        return;
    }
    // Prefer your constants over M_PI for portability
    constexpr double lattice_constant     = 5.43e-10;  // meters (Si)
    const double     normalization_factor = 2.0 * uepm::constants::pi / lattice_constant;

    const auto history = m_particle.get_history();

    history_file << "time,gamma,x,y,z,kx,ky,kz,vx,vy,vz,energy\n";
    for (std::size_t step = 0; step < history.get_number_of_steps(); ++step) {
        history_file << history.m_time_history[step] << "," << history.m_gammas[step] << "," << history.m_positions[step].x() << "," << history.m_positions[step].y() << ","
                     << history.m_positions[step].z() << "," << history.m_k_vectors[step].x() / normalization_factor << ","
                     << history.m_k_vectors[step].y() / normalization_factor << "," << history.m_k_vectors[step].z() / normalization_factor
                     << "," << history.m_velocities[step].x() << "," << history.m_velocities[step].y() << ","
                     << history.m_velocities[step].z() << "," << history.m_energies[step] << "\n";
    }
    history_file.close();
    std::cout << "Particle history exported to " << filename << "\n";
}

}  // namespace uepm::fbmc
