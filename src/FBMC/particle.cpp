/**
 * @file particle.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-09-10
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "particle.hpp"

#include <random>

#include "physical_constants.hpp"
#include "statistics_functions.hpp"
namespace uepm::fbmc {

particle::particle(std::size_t index, particle_type arg_particle_type, uepm::mesh_bz::ElectronPhonon* ptr_mesh_bz)
    : m_index(index),
      m_type(arg_particle_type),
      m_mesh_bz(ptr_mesh_bz) {}

/**
 * @brief Draw a new free flight time for the particle.
 *
 * @param p_gamma
 */

void particle::draw_free_flight_time(double p_gamma) {
    if (!(p_gamma > 0.0) || !std::isfinite(p_gamma)) {
        throw std::runtime_error("draw_free_flight_time: invalid p_gamma");
    }
    double u                   = m_random_distribution(m_random_generator);
    m_current_free_flight_time = -std::log(u) / p_gamma;
    m_time += m_current_free_flight_time;
    m_iter += 1;
    m_history.m_total_nb_steps += 1;
    update_position();
}

/**
 * @brief Update the k-vector of the particle based on the electric field.
 *
 * @param v_electric_field The electric field vector.
 */
void particle::update_k_vector(const vector3& v_electric_field) {
    m_k_vector += (get_charge_sign() * uepm::constants::q_e / uepm::constants::h_bar) * v_electric_field * m_current_free_flight_time;
}

void particle::update_group_velocity() {
    m_velocity = m_containing_bz_mesh_tetra->interpolate_gradient_energy_at_band(m_k_vector, m_band_index);
    m_velocity *= (1.0 / uepm::constants::h_bar_eV);
}

void particle::update_position() { m_position += m_velocity * m_current_free_flight_time; }

std::array<double, 8> particle::interpolate_phonon_scattering_rate_at_location(const vector3& location) {
    return m_containing_bz_mesh_tetra->interpolate_phonon_scattering_rate_at_location(location, m_band_index);
}

void particle::update_energy() { m_energy = m_containing_bz_mesh_tetra->interpolate_energy_at_band(m_k_vector, m_band_index); }

void particle::select_final_state_after_phonon_scattering(std::size_t idx_phonon_branch) {
    uepm::mesh_bz::SelectedFinalState Sf =
        m_mesh_bz->select_electron_phonon_final_state(m_band_index, m_k_vector, idx_phonon_branch, m_random_generator);
    m_k_vector                 = Sf.k_final;
    m_energy                   = Sf.E_final_eV;
    m_containing_bz_mesh_tetra = Sf.ptr_final_tetra;
    m_band_index               = Sf.idx_final_band;

    update_group_velocity();
}

double particle::compute_mean_energy() const {
    double weighted_sum = 0.0;
    double total_weight = 0.0;
    for (std::size_t i = 1; i < m_history.get_number_of_steps(); ++i) {
        double energy = m_history.m_energies[i];
        double weight = m_history.m_time_history[i] - m_history.m_time_history[i - 1];
        weighted_sum += energy * weight;
        total_weight += weight;
    }
    return (total_weight > 0.0) ? (weighted_sum / total_weight) : 0.0;
}

void particle::print_history_summary() const {
    std::size_t            total_events      = m_history.m_total_nb_steps;
    std::size_t            total_real_events = m_history.get_number_of_steps();
    std::size_t            nb_events_type    = m_history.m_scattering_events.size();
    std::array<double, 10> event_fractions   = {0.0};
    for (std::size_t i = 0; i < nb_events_type; ++i) {
        event_fractions[i] = 100.0 * static_cast<double>(m_history.m_scattering_events[i]) / static_cast<double>(total_real_events);
    }
    fmt::print("Particle {} history summary:\n", m_index);
    fmt::print("  Total recorded events: {}\n", total_events);
    for (std::size_t i = 0; i < nb_events_type; ++i) {
        fmt::print("\tEvent type {}: fraction = {:.4f}\n", i, event_fractions[i]);
    }
    fmt::print("\tEvent type 8: fraction = {:.2f}\n", event_fractions[8]);
    fmt::print("\tEvent type 9: fraction = {:.2f}\n", event_fractions[9]);

    double mean_energy = compute_mean_energy();
    fmt::print("  Mean energy: {:.6f} eV\n", mean_energy);
}

void particle::export_history_to_csv(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("export_history_to_csv: cannot open file for writing");
    }
    file << "time,gamma,x,y,z,kx,ky,kz,vx,vy,vz,energy,band_occupation\n";
    for (std::size_t i = 0; i < m_history.get_number_of_steps(); ++i) {
        file << m_history.m_time_history[i] << "," << m_history.m_gammas[i] << "," << m_history.m_positions[i].x() << ","
             << m_history.m_positions[i].y() << "," << m_history.m_positions[i].z() << "," << m_history.m_k_vectors[i].x() << ","
             << m_history.m_k_vectors[i].y() << "," << m_history.m_k_vectors[i].z() << "," << m_history.m_velocities[i].x() << ","
             << m_history.m_velocities[i].y() << "," << m_history.m_velocities[i].z() << "," << m_history.m_energies[i] << ","
             << m_history.m_band_occupations[i] << "\n";
    }
    file.close();
}

}  // namespace uepm::fbmc