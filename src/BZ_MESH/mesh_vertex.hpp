/**
 * @file mesh_vertex.hpp
 * @author Rémi Helleboid (remi.helleboid@st.com)
 * @brief
 * @version 0.1
 * @date 2022-07-14
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "vector.hpp"

namespace bz_mesh {

class Vertex {
 private:
    /**
     * @brief Index of the vertex within the mesh.
     *
     */
    std::size_t m_index;

    /**
     * @brief Position of the k-vector.
     *
     */
    vector3 m_position;

    /**
     * @brief The energy of the band with index b_idx at this vertex
     * is stored as m_band_energies[b_idx].
     * For example, m_band_energies[3] is the energy of the band 3 at the k-point with position m_position.
     *
     */
    std::vector<double> m_band_energies;

    /**
     * @brief The energy gradient at each band at this vertex.
     * The energy gradient is stored as a vector3 for each band.
     * For example, m_energy_gradient_at_bands[3] is the energy gradient of the band 3 at the k-point with position m_position.
     *
     */
    std::vector<vector3> m_energy_gradient_at_bands = {};

    /**
     * @brief Electron-phonon rates for each band.
     * The rates are stored as an array of 8 values:
     * [ALO, ALA, ATO, ATA, ELO, ELA, ETO, ETA] (ALO: absorption longitudinal optical, ALA: absorption longitudinal acoustic, ...)
     *
     */
    std::vector<std::array<double, 8>> m_electron_phonon_rates = {};

 public:
    /**
     * @brief Default constructor of a new Vertex object.
     *
     */
    Vertex() : m_index{0}, m_position{} {}

    /**
     * @brief Construct a new Vertex object with a given index.
     *
     * @param index
     */
    explicit Vertex(std::size_t index) : m_index(index), m_position{} {}

    /**
     * @brief Construct a new Vertex object with a given index and position.
     *
     * @param index
     * @param postion
     */
    Vertex(std::size_t index, const vector3& position) : m_index(index), m_position{position} {}

    /**
     * @brief Construct a new Vertex object with a given index and position (x, y, z).
     *
     * @param index
     * @param x
     * @param y
     * @param z
     */
    Vertex(std::size_t index, double x, double y, double z) : m_index(index), m_position{x, y, z} {}

    /**
     * @brief Get the index of the vertex.
     *
     * @return std::size_t
     */
    std::size_t get_index() const { return m_index; }

    /**
     * @brief Get the position of the vertex.
     *
     * @return vector3
     */
    const vector3& get_position() const { return m_position; }

    /**
     * @brief Set the position of the vertex.
     *
     * @param position
     */
    void set_position(const vector3& position) { m_position = position; }

    /**
     * @brief Shift the position of the vertex by a given vector.
     *
     * @param shift
     */
    void shift_position(const vector3& shift) { m_position += shift; }

    /**
     * @brief Add an energy at the end of the list of energy. (One energy per band).
     *
     * @param energy
     */
    void add_band_energy_value(double energy) { m_band_energies.push_back(energy); }

    void remove_band_energy(std::size_t band_index) {
        if (band_index >= m_band_energies.size()) {
            throw std::invalid_argument("The energy of band " + std::to_string(band_index) +
                                        " cannot be removed because it does not exists.");
        }
        m_band_energies.erase(m_band_energies.begin() + band_index);
    }

    /**
     * @brief Set the band energy at a given band index.
     *
     * @param index_band
     * @param new_energy
     */
    void set_band_energy(std::size_t index_band, double new_energy) {
        if (index_band >= m_band_energies.size()) {
            throw std::invalid_argument("The energy of band " + std::to_string(index_band) +
                                        " cannot be modify because it does not exist.");
        }
        m_band_energies[index_band] = new_energy;
    }

    /**
     * @brief Swap the energy of two bands.
     *
     * @param old_index
     * @param new_index
     */
    void swap_bands(std::size_t old_index, std::size_t new_index) {
        if (old_index >= m_band_energies.size() || new_index >= m_band_energies.size()) {
            throw std::invalid_argument("The energy of band " + std::to_string(old_index) + " or " + std::to_string(new_index) +
                                        " cannot be swapped because it does not exists.");
        }
        std::swap(m_band_energies[old_index], m_band_energies[new_index]);
    }

    /**
     * @brief Get the number of band energies stored in the Vertex object.
     *
     * @return std::size_t
     */
    std::size_t get_number_bands() const { return m_band_energies.size(); }

    /**
     * @brief Get the energy at band index.
     *
     * @param band_index
     * @return double
     */
    double get_energy_at_band(std::size_t band_index) const { return m_band_energies[band_index]; }

    /**
     * @brief Get the list of band energies.
     *
     * @return const std::vector<double>&
     */
    const std::vector<double>& get_band_energies() const { return m_band_energies; }

    const std::vector<vector3>& get_energy_gradient_at_bands() const { return m_energy_gradient_at_bands; }

    const vector3& get_energy_gradient_at_band(std::size_t band_index) const { return m_energy_gradient_at_bands[band_index]; }
    void           set_energy_gradient_at_bands(const std::vector<vector3>& gradients) { m_energy_gradient_at_bands = gradients; }
    void           set_energy_gradient_at_band(std::size_t band_index, const vector3& gradient) {
        if (band_index >= m_energy_gradient_at_bands.size()) {
            throw std::invalid_argument("The energy gradient of band " + std::to_string(band_index) +
                                        " cannot be modify because it does not exist.");
        }
        m_energy_gradient_at_bands[band_index] = gradient;
    }
    void resize_energy_gradient_at_bands(std::size_t new_size) { m_energy_gradient_at_bands.resize(new_size); }
    void push_back_energy_gradient_at_band(const vector3& gradient) { m_energy_gradient_at_bands.push_back(gradient); }

    /**
     * @brief Add the electron-phonon rates for a given band.
     *
     * @param rates
     */
    void add_electron_phonon_rates(const std::array<double, 8>& rates) { m_electron_phonon_rates.push_back(rates); }

    /**
     * @brief Get the electron-phonon rates for a given band.
     *
     * @param band_index
     * @return const std::array<double, 8>&
     */
    const std::array<double, 8>& get_electron_phonon_rates(std::size_t band_index) const { return m_electron_phonon_rates[band_index]; }

    const std::vector<std::array<double, 8>>& get_electron_phonon_rates_all_bands() const { return m_electron_phonon_rates; }
};

}  // namespace bz_mesh