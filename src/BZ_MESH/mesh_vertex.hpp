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
#include <vector>

#include "vector.hpp"

namespace bz_mesh {

class Vertex {
 private:
    std::size_t m_index;

    vector3 m_position;

    /**
     * @brief The energy of the conduction band with index b_idx at this vertex
     * is stored as m_conduction_band_energies[b_idx]
     *
     */
    std::vector<double> m_conduction_band_energies;

    /**
     * @brief The energy of the valance band with index b_idx at this vertex
     * is stored as m_conduction_band_energies[b_idx]
     *
     */
    std::vector<double> m_valance_band_energies;

 public:
    Vertex() : m_index{0}, m_position{} {}
    Vertex(std::size_t index) : m_index(index), m_position{} {}
    Vertex(std::size_t index, vector3 postion) : m_index(index), m_position{} {}
    Vertex(std::size_t index, double x, double y, double z) : m_index(index), m_position{x, y, z} {}

    std::size_t    get_index() const { return m_index; }
    const vector3& get_position() const { return m_position; }

    std::size_t                get_number_valence_bands() const { return m_valance_band_energies.size(); }
    std::size_t                get_number_conduction_bands() const { return m_conduction_band_energies.size(); }
    const std::vector<double>& get_valence_energies() const { return m_valance_band_energies; }
    const std::vector<double>& get_conduction_energies() const { return m_conduction_band_energies; }
    double get_energy_at_conduction_band(std::size_t band_index) const { return m_conduction_band_energies[band_index]; }
    double get_energy_at_valance_band(std::size_t band_index) const { return m_valance_band_energies[band_index]; }

};

}  // namespace bz_mesh