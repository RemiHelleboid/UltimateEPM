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
    double u = m_random_distribution(m_random_generator);  // [0,1)
    // Map to (0,1] robustly to avoid log(0):
    if (u <= 0.0) {
        u = std::numeric_limits<double>::min();  // needs <limits>
    }
    if (u >= 1.0) {
        u = std::nextafter(1.0, 0.0);  // needs <cmath>
    }
    m_current_free_flight_time = -std::log(u) / p_gamma;
    m_time += m_current_free_flight_time;
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
    m_velocity = m_containing_bz_mesh_tetra->get_gradient_energy_at_band(m_band_index);
    m_velocity *= (1.0 / uepm::constants::h_bar_eV);
}

std::array<double, 8> particle::interpolate_phonon_scattering_rate_at_location(const vector3& location) {
    return m_mesh_bz->interpolate_phonon_scattering_rate_at_location(location, m_band_index);
}

void particle::update_energy() { m_energy = m_containing_bz_mesh_tetra->interpolate_energy_at_band(m_k_vector, m_band_index); }

void particle::select_final_state_after_phonon_scattering(std::size_t idx_phonon_branch) {
    uepm::mesh_bz::SelectedFinalState Sf =
        m_mesh_bz->select_electron_phonon_final_state(m_band_index, m_k_vector, idx_phonon_branch, m_random_generator);
    m_k_vector                 = Sf.k_final;
    m_energy                   = Sf.E_final_eV;
    m_containing_bz_mesh_tetra = Sf.tetra_ptr;

    update_group_velocity();
}

}  // namespace uepm::fbmc