/**
 * @file single_part_fbmc.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2025-09-19
 * 
 * @copyright Copyright (c) 2025
 * 
 */


#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <vector>

#include "Constants.hpp"
#include "single_part_fbmc.hpp"
namespace uepm::fbmc {

Single_particle_simulation::Single_particle_simulation(uepm::mesh_bz::ElectronPhonon*     ptr_mesh_bz,
                                                       const Bulk_environment&      bulk_env,
                                                       const Simulation_parameters& sim_params)
    : m_ptr_mesh_bz(ptr_mesh_bz),
      m_bulk_env(bulk_env),
      m_sim_params(sim_params) {
    double m_time = 0.0;
    // Initialize the particle
    const std::size_t index = 0;
    m_particle              = particle(index, particle_type::electron, m_ptr_mesh_bz);
    m_particle.set_position({0.0, 0.0, 0.0});
    const double thermal_energy = bulk_env.m_temperature * uepm::pseudopotential::Constants::k_b_eV;
    m_particle.set_energy(thermal_energy);
    m_particle.set_velocity({0.0, 0.0, 0.0});
    std::cout << "Thermal energy at " << bulk_env.m_temperature << " K: " << thermal_energy << " eV" << std::endl;
    vector3 initial_k = m_ptr_mesh_bz->draw_random_k_point_at_energy(thermal_energy, 0, m_particle.get_random_generator());
    m_particle.set_k_vector(initial_k);
    std::cout << "Initial k-vector (drawn at thermal energy): " << initial_k << std::endl;
    uepm::mesh_bz::Tetra* containing_tetra = m_ptr_mesh_bz->find_tetra_at_location(m_particle.get_k_vector());
    if (containing_tetra == nullptr) {
        throw std::runtime_error("Initial k-point is out of the Brillouin zone mesh.");
    }
    m_particle.set_containing_bz_mesh_tetra(containing_tetra);
    m_particle.update_group_velocity();
    m_particle.update_energy();

    double init_energy_true = containing_tetra->interpolate_energy_at_band(m_particle.get_k_vector(), 0);

    std::cout << "Initial k-vector: " << m_particle.get_k_vector().x() << " " << m_particle.get_k_vector().y() << " "
              << m_particle.get_k_vector().z() << std::endl;
    std::cout << "Initial energy (set): " << m_particle.get_energy() << " eV" << std::endl;
    std::cout << "Initial energy (true): " << init_energy_true << " eV" << std::endl;
}

void Single_particle_simulation::run_simulation() {
    // Freefligts time in a file
    std::ofstream free_flight_time_file("free_flight_times.txt");
    if (!free_flight_time_file.is_open()) {
        std::cerr << "Could not open free flight time file for writing." << std::endl;
        return;
    }

    // Compute total scattering rate
    const double p_gamma = m_ptr_mesh_bz->compute_P_Gamma();
    std::cout << "Total scattering rate: " << p_gamma << " 1/s" << std::endl;
    if (p_gamma <= 0.0) {
        std::cerr << "Total scattering rate is zero or negative." << std::endl;
        return;
    }

    int nb_foldings = 0;

    while (m_time < m_sim_params.m_simulation_time) {
        // std::cout << "Time: " << m_time << " s" << std::endl;
        m_particle.update_history();

        // Draw free flight time
        m_particle.draw_free_flight_time(p_gamma);
        free_flight_time_file << m_particle.get_current_free_flight_time() << std::endl;

        // Update time
        m_time += m_particle.get_current_free_flight_time();

        // Update k-vector
        m_particle.update_k_vector(m_bulk_env.m_electric_field);
        // std::cout << "New k-vector: " << m_particle.get_k_vector().x() << " " << m_particle.get_k_vector().y() << " "
        //           << m_particle.get_k_vector().z() << std::endl;

        // Find containing tetrahedron
        uepm::mesh_bz::Tetra* containing_tetra = m_ptr_mesh_bz->find_tetra_at_location(m_particle.get_k_vector());
        if (containing_tetra == nullptr) {
            // TEST
            bool is_inside_bz = m_ptr_mesh_bz->inside_ws_bcc(m_particle.get_k_vector());
            if (is_inside_bz) {
                std::cout << "ERROR : Particle is inside the Wigner-Seitz cell but outside the mesh." << std::endl;
            } else {
                std::cout << "GOOD : Particle is outside the Wigner-Seitz cell and the mesh." << std::endl;
            }

            std::cerr << "Particle exited the Brillouin zone mesh ("
                      << "k = " << m_particle.get_k_vector().x() << " " << m_particle.get_k_vector().y() << " "
                      << m_particle.get_k_vector().z() << "). Folding back..." << std::endl;
            // Fold back into the first Brillouin zone
            vector3 folded_k = m_ptr_mesh_bz->fold_ws_bcc(m_particle.get_k_vector());
            m_particle.set_k_vector(folded_k);
            containing_tetra = m_ptr_mesh_bz->find_tetra_at_location(m_particle.get_k_vector());
            if (containing_tetra == nullptr) {
                std::cerr << "Particle still outside the Brillouin zone mesh after folding back." << std::endl;
                break;
            } else {
                std::cout << "Particle folded back into the Brillouin zone mesh at k = " << folded_k.x() << " " << folded_k.y() << " "
                          << folded_k.z() << std::endl;
                nb_foldings++;
            }
        }
        m_particle.set_containing_bz_mesh_tetra(containing_tetra);

        m_particle.update_energy();

        // Update group velocity
        m_particle.update_group_velocity();

        // Update position
        vector3 new_position = m_particle.get_position() + m_particle.get_velocity() * m_particle.get_current_free_flight_time();
        m_particle.set_position(new_position);

        // Scatter the particle
        std::array<double, 8> scattering_rates      = m_particle.interpolate_phonon_scattering_rate_at_location(m_particle.get_k_vector());
        double                total_scattering_rate = std::accumulate(scattering_rates.begin(), scattering_rates.end(), 0.0);
        if (total_scattering_rate <= 0.0) {
            std::cerr << "Total scattering rate is zero or negative during scattering." << std::endl;
            break;
        }
        std::uniform_real_distribution<double> distribution(0.0, total_scattering_rate);
        double                                 random_value     = distribution(m_particle.get_random_generator());
        double                                 cumulative_rate  = 0.0;
        std::size_t                            scattering_event = 0;
        bool                                   event_found      = false;
        for (std::size_t idx_mode = 0; idx_mode < scattering_rates.size(); ++idx_mode) {
            cumulative_rate += scattering_rates[idx_mode];
            if (random_value <= cumulative_rate) {
                scattering_event = idx_mode;
                event_found      = true;
                break;
            }
        }
        // std::cout << "Scattering event: " << scattering_event << " at energy " << m_particle.get_energy() << " eV" << std::endl;
        if (!event_found) {
            std::cout << "Self-scattering event." << std::endl;
            continue;
        }
    }
    // auto indx_final_band_idx_final_tetra = 
    std::cout << "Simulation finished at time " << m_time << " s" << std::endl;
    std::cout << "Number of foldings: " << nb_foldings << std::endl;
    free_flight_time_file.close();
}

void Single_particle_simulation::export_history(const std::string& filename) {
    std::ofstream history_file(filename);
    if (!history_file.is_open()) {
        std::cerr << "Could not open history file for writing." << std::endl;
        return;
    }
    constexpr double lattice_constant     = 5.43e-10;  // in meters (for silicon)
    constexpr double normalization_factor = 2.0 * M_PI / lattice_constant;

    // Write particle history to file
    auto history = m_particle.get_history();
    history_file << "time,x,y,z,kx,ky,kz,vx,vy,vz,energy\n";
    for (std::size_t step = 0; step < history.get_number_of_steps(); ++step) {
        history_file << history.m_time_history[step] << "," << history.m_positions[step].x() << "," << history.m_positions[step].y() << ","
                     << history.m_positions[step].z() << "," << history.m_k_vectors[step].x() / normalization_factor << ","
                     << history.m_k_vectors[step].y() / normalization_factor << "," << history.m_k_vectors[step].z() / normalization_factor
                     << "," << history.m_velocities[step].x() << "," << history.m_velocities[step].y() << ","
                     << history.m_velocities[step].z() << "," << history.m_energies[step] << "\n";
    }
    history_file.close();
    std::cout << "Particle history exported to " << filename << std::endl;
}

}  // namespace uepm::fbmc