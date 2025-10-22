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
#include "vector.hpp"
// #include "particle_history_mc.hpp"
#include "electron_phonon.hpp"
#include "mesh_tetra.hpp"
#include "vector.hpp"

namespace uepm::fbmc {

enum class particle_type { electron = -1, hole = 1 };

using vector3 = uepm::mesh_bz::vector3;

struct scattering_rate {
    std::array<double, 8> m_phonon_rate;
    double                m_phonon_rate_total      = 0.0;
    double                m_impurity_rate          = 0.0;
    double                m_impact_ionization_rate = 0.0;
};

struct particle_history {
    std::size_t          m_index;
    std::vector<double>  m_time_history;
    std::vector<vector3> m_positions;
    std::vector<vector3> m_k_vectors;
    std::vector<vector3> m_velocities;
    std::vector<double>  m_energies;
    std::vector<double>  m_gammas;

    particle_history() : m_index(0), m_positions(), m_k_vectors(), m_velocities(), m_energies(), m_gammas() {}
    particle_history(std::size_t index) : m_index(index), m_positions(), m_k_vectors(), m_velocities(), m_energies(), m_gammas() {}

    void reserve(std::size_t n_steps) {
        m_time_history.reserve(n_steps);
        m_positions.reserve(n_steps);
        m_k_vectors.reserve(n_steps);
        m_velocities.reserve(n_steps);
        m_energies.reserve(n_steps);
        m_gammas.reserve(n_steps);
    }
    
    void add_time(double time) { m_time_history.push_back(time); }
    void add_position(const vector3& position) { m_positions.push_back(position); }
    void add_k_vector(const vector3& k_vector) { m_k_vectors.push_back(k_vector); }
    void add_velocity(const vector3& velocity) { m_velocities.push_back(velocity); }
    void add_energy(double energy) { m_energies.push_back(energy); }
    void add_gamma(double gamma) { m_gammas.push_back(gamma); }

    void add_particle_state(const double&  my_time,
                            const vector3& my_position,
                            const vector3& my_k_vector,
                            const vector3& my_velocity,
                            double         my_energy,
                            double         my_gamma) {
        add_time(my_time);
        add_position(my_position);
        add_k_vector(my_k_vector);
        add_velocity(my_velocity);
        add_energy(my_energy);
        add_gamma(my_gamma);
    }
    std::size_t get_number_of_steps() const { return m_positions.size(); }
};

class particle {
 protected:
    /**
     * @brief Pointer to the Brillouin Zone Mesh.
     *
     */
    uepm::mesh_bz::ElectronPhonon* m_mesh_bz = nullptr;

    /**
     * @brief Index of the particle in the simulation.
     *
     */
    std::size_t m_index = 0;

    /**
     * @brief Type of the particle (electron or hole).
     *
     */
    particle_type m_type = particle_type::electron;

    /**
     * @brief Time of the particle in the simulation.
     *
     */
    double m_time = 0.0;

    /**
     * @brief Position of the particle in the real space.
     *
     */
    vector3 m_position = {0.0, 0.0, 0.0};

    /**
     * @brief Wave vector of the particle.
     *
     */
    vector3 m_k_vector = {0.0, 0.0, 0.0};

    /**
     * @brief Band index of the particle.
     *
     */
    int m_band_index = 0;

    /**
     * @brief Velocity of the particle.
     *
     */
    vector3 m_velocity = {0.0, 0.0, 0.0};

    /**
     * @brief Energy of the particle.
     *
     */
    double m_energy = 0.0;

    /**
     * @brief Free flight time of the particle.
     *
     */
    double m_current_free_flight_time = 0.0;

    /**
     * @brief Scattering rate gamma of the particle.
     *
     */
    double m_gamma = 0.0;

    /**
     * @brief Pointer to the tetrahedron in which the particle lies.
     *
     */
    uepm::mesh_bz::Tetra* m_containing_bz_mesh_tetra = nullptr;

    /**
     * @brief Random number generator.
     *
     */
    std::mt19937 m_random_generator = std::mt19937(std::random_device{}());

    /**
     * @brief Random number distribution.
     *
     */
    std::uniform_real_distribution<double> m_random_distribution = std::uniform_real_distribution<double>(0.0, 1.0);

    /**
     * @brief History of the particle during the simulation.
     *
     */
    particle_history m_history;

 public:
    particle() = default;
    particle(std::size_t index, particle_type type, uepm::mesh_bz::ElectronPhonon* mesh);
    particle(const particle& other)            = default;
    particle& operator=(const particle& other) = default;
    ~particle()                                = default;

    std::size_t           get_index() const { return m_index; }
    void                  set_index(std::size_t index) { m_index = index; }
    particle_type         get_type() const { return m_type; }
    double                get_charge_sign() const { return static_cast<double>(m_type); }
    double                get_time() const { return m_time; }
    void                  set_time(double time) { m_time = time; }
    int                   get_band_index() const { return m_band_index; }
    void                  set_band_index(int band_index) { m_band_index = band_index; }
    const vector3&        get_position() const { return m_position; }
    void                  set_position(const vector3& position) { m_position = position; }
    const vector3&        get_k_vector() const { return m_k_vector; }
    void                  set_k_vector(const vector3& k_vector) { m_k_vector = k_vector; }
    const vector3&        get_velocity() const { return m_velocity; }
    void                  set_velocity(const vector3& velocity) { m_velocity = velocity; }
    double                get_gamma() const { return m_gamma; }
    void                  set_gamma(double gamma) { m_gamma = gamma; }
    double                get_energy() const { return m_energy; }
    void                  set_energy(double energy) { m_energy = energy; }
    double                get_current_free_flight_time() const { return m_current_free_flight_time; }
    uepm::mesh_bz::Tetra* get_containing_bz_mesh_tetra() const { return m_containing_bz_mesh_tetra; }
    void                  set_containing_bz_mesh_tetra(uepm::mesh_bz::Tetra* containing_bz_mesh_tetra) {
        m_containing_bz_mesh_tetra = containing_bz_mesh_tetra;
    }
    void set_random_generator(std::mt19937 random_generator) { m_random_generator = random_generator; }
    void reserve_history(std::size_t n_steps) { m_history.reserve(n_steps); }

    std::array<double, 8> interpolate_phonon_scattering_rate_at_location(const vector3& location);
    void                  compute_post_phonon_scattering_state();

    void draw_free_flight_time(double p_gamma);
    void update_k_vector(const vector3& v_electric_field);
    void update_energy();
    void update_group_velocity();

    void select_final_state_after_phonon_scattering(std::size_t idx_phonon_branch);

    std::mt19937& get_random_generator() { return m_random_generator; }
    void          draw_random_k_point_at_energy(double energy, std::size_t idx_band) {
        m_mesh_bz->draw_random_k_point_at_energy(energy, idx_band, get_random_generator());
    }
    void update_history() { m_history.add_particle_state(m_time, m_position, m_k_vector, m_velocity, m_energy, m_gamma); }
    const particle_history& get_history() const { return m_history; }
    void                    reset_history() { m_history = particle_history(m_index); }
};

}  // namespace uepm::fbmc