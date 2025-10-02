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

#include "Constants.hpp"
namespace fbmc {

particle::particle(std::size_t index, particle_type arg_particle_type, bz_mesh::ElectronPhonon* ptr_mesh_bz)
    : m_index(index),
      m_type(arg_particle_type),
      m_mesh_bz(ptr_mesh_bz) {}

/**
 * @brief Draw a new free flight time for the particle.
 *
 * @param p_gamma
 */
void particle::draw_free_flight_time(double p_gamma) {
    m_current_free_flight_time = -(1.0 / p_gamma) * std::log(m_random_distribution(m_random_generator));
    m_time += m_current_free_flight_time;
}

/**
 * @brief Update the k-vector of the particle based on the electric field.
 *
 * @param v_electric_field The electric field vector.
 */
void particle::update_k_vector(const vector3& v_electric_field) {
    m_k_vector += (get_charge_sign() * EmpiricalPseudopotential::Constants::q_e / EmpiricalPseudopotential::Constants::h_bar) *
                  v_electric_field * m_current_free_flight_time;
}

void particle::update_group_velocity() {
    m_velocity = m_containing_bz_mesh_tetra->get_gradient_energy_at_band(m_band_index);
    m_velocity *= (1.0 / EmpiricalPseudopotential::Constants::h_bar_eV);
}

std::array<double, 8> particle::interpolate_phonon_scattering_rate_at_location(const vector3& location) {
    return m_mesh_bz->interpolate_phonon_scattering_rate_at_location(location, m_band_index);
}

void particle::update_energy() { m_energy = m_containing_bz_mesh_tetra->interpolate_energy_at_band(m_k_vector, m_band_index); }

void compute_post_phonon_scattering_state() {}

}  // namespace fbmc