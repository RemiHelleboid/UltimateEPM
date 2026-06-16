/**
 * @file particle.hpp
 * @author your name (you@domain.com)
 * @brief Particle simulation class for the FBMC (Full Band Monte Carlo) method.
 * @version 0.1
 * @date 2025-09-19
 *
 * @copyright Copyright (c) 2025
 *
 */

#pragma once

#include <array>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include "bz_mesh.hpp"
// #include "particle_history_mc.hpp"
#include "electron_phonon.hpp"
#include "mesh_tetra.hpp"
#include "vector_bz.hpp"

namespace uepm::fbmc {

enum class particle_type { electron = -1, hole = 1 };

using vector3 = uepm::mesh_bz::vector3;

struct scattering_rate {
    std::array<double, 8> m_phonon_rate;
    double                m_phonon_rate_total      = 0.0;
    double                m_impurity_rate          = 0.0;
    double                m_impact_ionization_rate = 0.0;
};

struct particle_state {
    double      m_time             = 0.0;
    std::size_t m_iter             = 0;
    vector3     m_position         = vector3(0.0, 0.0, 0.0);
    vector3     m_k_vector         = vector3(0.0, 0.0, 0.0);
    vector3     m_velocity         = vector3(0.0, 0.0, 0.0);
    double      m_energy           = 0.0;
    double      m_gamma            = 0.0;
    double      m_free_flight_time = 0.0;
    int         m_band_index       = 0;
};

struct particle_history {
    std::size_t m_index_particle = 0;
    // Total number of steps (including self-scattering)
    std::size_t m_total_nb_steps = 0;
    // Recorded data at each step (not including self-scattering)

    std::vector<double>  m_time_history;
    std::vector<vector3> m_positions;
    std::vector<vector3> m_k_vectors;
    std::vector<vector3> m_velocities;
    std::vector<double>  m_energies;
    std::vector<double>  m_gammas;

    // 8 phonon branches + impact ionization + self-scattering
    std::array<std::size_t, 10> m_scattering_events = {0};
    std::vector<double>         m_band_occupations;

    particle_history() : m_index_particle(0) {}
    particle_history(std::size_t index) : m_index_particle(index) {}

    void reserve(std::size_t n_steps) {
        m_time_history.reserve(n_steps);
        m_positions.reserve(n_steps);
        m_k_vectors.reserve(n_steps);
        m_velocities.reserve(n_steps);
        m_energies.reserve(n_steps);
        m_gammas.reserve(n_steps);
        m_band_occupations.reserve(n_steps);
    }

    void add_particle_state(const particle_state& state) {
        m_time_history.push_back(state.m_time);
        m_positions.push_back(state.m_position);
        m_k_vectors.push_back(state.m_k_vector);
        m_velocities.push_back(state.m_velocity);
        m_energies.push_back(state.m_energy);
        m_gammas.push_back(state.m_gamma);
        m_band_occupations.push_back(static_cast<double>(state.m_band_index));
    }
    std::size_t get_number_of_steps() const { return m_positions.size(); }

    /**
     * @brief Add a scattering event to the history.
     *
     * @param event_index
     */
    void add_event(std::size_t event_index) {
        if (event_index < m_scattering_events.size()) {
            m_scattering_events[event_index]++;
        }
    }
};

class particle {
 protected:
    particle_state                         m_state;
    uepm::mesh_bz::ElectronPhonon*         m_mesh_bz                  = nullptr;
    std::size_t                            m_index                    = 0;
    particle_type                          m_type                     = particle_type::electron;
    uepm::mesh_bz::Tetra*                  m_containing_bz_mesh_tetra = nullptr;
    std::mt19937                           m_random_generator         = std::mt19937(std::random_device{}());
    std::uniform_real_distribution<double> m_random_distribution = std::uniform_real_distribution<double>(0.0, 1.0);
    particle_history                       m_history;

 public:
    particle() = default;
    particle(std::size_t index, particle_type type, uepm::mesh_bz::ElectronPhonon* mesh);
    particle(const particle& other)            = default;
    particle& operator=(const particle& other) = default;
    ~particle()                                = default;

    std::size_t           get_index() const { return m_index; }
    void                  set_index(std::size_t index) { m_index = index; }
    particle_type         get_type() const { return m_type; }
    double                get_signed_charge() const { return static_cast<double>(m_type); }
    const particle_state& state() const { return m_state; }
    particle_state&       state() { return m_state; }
    uepm::mesh_bz::Tetra* get_containing_bz_mesh_tetra() const { return m_containing_bz_mesh_tetra; }
    void                  set_containing_bz_mesh_tetra(uepm::mesh_bz::Tetra* containing_bz_mesh_tetra) {
        m_containing_bz_mesh_tetra = containing_bz_mesh_tetra;
    }
    void set_random_generator(std::mt19937 random_generator) { m_random_generator = random_generator; }
    void reserve_history(std::size_t n_steps) { m_history.reserve(n_steps); }

    std::array<double, 8> interpolate_phonon_scattering_rate_at_location(const vector3& location);
    void                  compute_post_phonon_scattering_state();

    void draw_free_flight_time(double p_gamma);
    void advance_time_by_free_flight();

    void update_k_vector(const vector3& v_electric_field);
    void update_energy();
    void update_group_velocity();

    void update_position();
    void update_position(const vector3& velocity, double dt);

    void select_final_state_after_phonon_scattering(std::size_t idx_phonon_branch);
    void select_final_state_after_impact_ionization(double energy_threshold_eV);

    std::mt19937& get_random_generator() { return m_random_generator; }
    void          draw_random_k_point_at_energy(double energy, std::size_t idx_band) {
        m_mesh_bz->draw_random_k_point_at_energy(energy, idx_band, get_random_generator());
    }
    void update_history() { m_history.add_particle_state(m_state); }
    void add_scattering_event_to_history(std::size_t event_index) { m_history.add_event(event_index); }
    void print_history_summary() const;
    const particle_history& get_history() const { return m_history; }
    void                    reset_history() { m_history = particle_history(m_index); }
    void                    export_history_to_csv(const std::string& filename) const;
    double                  compute_mean_energy() const;
    double                  compute_mean_energy(double start_time_s) const;
    double                  extract_impact_ionization_coeff() const;
    double                  extract_global_average_velocity() const;
    double                  extract_global_average_velocity(double start_time_s) const;
};

}  // namespace uepm::fbmc
