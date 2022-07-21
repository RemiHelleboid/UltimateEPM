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

 public:
    Vertex() : m_index{0}, m_position{} {}
    explicit Vertex(std::size_t index) : m_index(index), m_position{} {}
    Vertex(std::size_t index, const vector3& postion) : m_index(index), m_position{} {}
    Vertex(std::size_t index, double x, double y, double z) : m_index(index), m_position{x, y, z} {}

    std::size_t    get_index() const { return m_index; }
    const vector3& get_position() const { return m_position; }

    void add_band_energy_value(double energy) { m_band_energies.push_back(energy); }
    void set_band_energy(std::size_t index_band, double new_energy) {
        if (index_band > m_band_energies.size()) {
            throw std::invalid_argument("The energy of valence band " + std::to_string(index_band) +
                                        " cannot be modify because it does not exists.");
        }
        m_band_energies[index_band] = new_energy;
    }
    std::size_t get_number_bands() const { return m_band_energies.size(); }
    double      get_energy_at_band(std::size_t band_index) const { return m_band_energies[band_index]; }
};

}  // namespace bz_mesh