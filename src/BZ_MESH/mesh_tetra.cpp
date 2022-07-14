/**
 * @file mesh_tetra.cpp
 * @author Rémi Helleboid (remi.helleboid@st.com)
 * @brief
 * @version 0.1
 * @date 2022-07-14
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "mesh_tetra.hpp"

#include <algorithm>
#include <array>
#include <vector>

namespace bz_mesh {

/**
 * @brief Construct a new Tetra by passing directly the array of the four pointers to the vertices.
 *
 * @param list_vertices
 */
Tetra::Tetra(const std::array<Vertex*, 4>& list_vertices)
    : m_list_vertices(list_vertices),
      m_nb_conduction_bands{m_list_vertices[0]->get_number_conduction_bands()},
      m_nb_valence_bands{m_list_vertices[0]->get_number_valence_bands()} {
    m_list_edges[0] = compute_edge(1, 0);
    m_list_edges[1] = compute_edge(2, 0);
    m_list_edges[2] = compute_edge(3, 0);
    m_list_edges[3] = compute_edge(2, 1);
    m_list_edges[4] = compute_edge(3, 1);
    m_list_edges[5] = compute_edge(3, 2);
    m_signed_volume = compute_signed_volume();
}

double Tetra::compute_signed_volume() const {
    return (1.0 / 6.0) * scalar_triple_product(m_list_edges[0], m_list_edges[1], m_list_edges[2]);
}

/**
 * @brief Return the values of the energy of the index_band valence band at the 4 vertices of the tetrahedra.
 *
 * @param index_band
 * @return std::vector<double>
 */
std::vector<double> Tetra::get_valence_band_energies_at_vtx(std::size_t index_band) const {
    std::vector<double> list_energies_valence;
    for (auto& p_vtx : m_list_vertices) {
        list_energies_valence.push_back(p_vtx->get_energy_at_valance_band(index_band));
    }
    return list_energies_valence;
}

/**
 * @brief Return the values of the energy of the index_band conduction band at the 4 vertices of the tetrahedra.
 *
 * @param index_band
 * @return std::vector<double>
 */
std::vector<double> Tetra::get_conduction_band_energies_at_vtx(std::size_t index_band) const {
    std::vector<double> list_energies_conduction;
    for (auto& p_vtx : m_list_vertices) {
        list_energies_conduction.push_back(p_vtx->get_energy_at_conduction_band(index_band));
    }
    return list_energies_conduction;
}

/**
 * @brief Compute the minimum and maximum energy for each band.
 * The results are set in the corresponding class members.
 *
 */
void Tetra::compute_and_set_minmax_energies() {
    min_energy_per_valence_band.clear();
    max_energy_per_valence_band.clear();
    min_energy_per_conduction_band.clear();
    max_energy_per_conduction_band.clear();
    for (std::size_t band_index = 0; band_index < m_nb_valence_bands; ++band_index) {
        std::vector<double> energies_at_vertices = get_valence_band_energies_at_vtx(band_index);
        auto                min_max_energy       = std::minmax(energies_at_vertices.begin(), energies_at_vertices.end());
        min_energy_per_valence_band.push_back(*min_max_energy.first);
        min_energy_per_valence_band.push_back(*min_max_energy.second);
    }
    for (std::size_t band_index = 0; band_index < m_nb_conduction_bands; ++band_index) {
        std::vector<double> energies_at_vertices = get_conduction_band_energies_at_vtx(band_index);
        auto                min_max_energy       = std::minmax(energies_at_vertices.begin(), energies_at_vertices.end());
        min_energy_per_conduction_band.push_back(*min_max_energy.first);
        max_energy_per_conduction_band.push_back(*min_max_energy.second);
    }
}

/**
 * @brief Compute the edge vector between two vertices of the tetrahedra.
 * The result is: vtx_1 - vtx_2.
 *
 * @param index_vtx_1
 * @param index_vtx_2
 * @return vector3
 */
vector3 Tetra::compute_edge(std::size_t index_vtx_1, std::size_t index_vtx_2) const {
    if (index_vtx_1 > 3 || index_vtx_2 > 3) {
        throw std::invalid_argument("In Tetra::compute_edge, the index of vertex must be between 0 and 3.");
    }
    return m_list_vertices[index_vtx_1]->get_position() - m_list_vertices[index_vtx_2]->get_position();
}

/**
 * @brief Compute the barycentric coordinate of a given location within the tetrahedra.
 * The returned array of size 4 contains the barycentric coordinates with respect to the vertices in the following order :
 *  0, 1, 2 and 3, respectively.
 *
 * @warning warning message: Do not use this function to check if the location lies in the tetrahedra,
 * The computation relies on the hypothesis that the location do lies in it. Use Tetra::is_location_inside instead.
 *
 * @param location
 * @return std::array<double, 4>
 */
std::array<double, 4> Tetra::compute_barycentric_coordinates(const vector3& location) const {
    const vector3 v_loc1            = location - m_list_vertices[0]->get_position();
    const double  tetra_determinant = 6.0 * m_signed_volume;
    const double  lambda_2          = scalar_triple_product(v_loc1, m_list_edges[1], m_list_edges[2]) / tetra_determinant;
    const double  lambda_3          = scalar_triple_product(v_loc1, m_list_edges[2], m_list_edges[0]) / tetra_determinant;
    const double  lambda_4          = scalar_triple_product(v_loc1, m_list_edges[0], m_list_edges[1]) / tetra_determinant;
    const double  lambda_1          = 1.0 - lambda_2 - lambda_3 - lambda_4;
    return {lambda_1, lambda_2, lambda_3, lambda_4};
}

/**
 * @brief Check if a given location lies inside the tetrahedra.
 *
 * @param location
 * @return true
 * @return false
 */
bool Tetra::is_location_inside(const vector3& location) const {
    const vector3 v_loc1            = location - m_list_vertices[0]->get_position();
    const vector3 v_loc2            = location - m_list_vertices[1]->get_position();
    const double  tetra_determinant = 6.0 * m_signed_volume;
    const double  lambda_1          = scalar_triple_product(v_loc2, m_list_edges[4], m_list_edges[3]) / tetra_determinant;
    const double  lambda_2          = scalar_triple_product(v_loc1, m_list_edges[1], m_list_edges[2]) / tetra_determinant;
    const double  lambda_3          = scalar_triple_product(v_loc1, m_list_edges[2], m_list_edges[0]) / tetra_determinant;
    const double  lambda_4          = scalar_triple_product(v_loc1, m_list_edges[0], m_list_edges[1]) / tetra_determinant;
    return (lambda_1 >= 0 && lambda_2 >= 0 && lambda_3 >= 0 && lambda_4 >= 0);
}

/**
 * @brief Compute the euclidean position from barycentric coordinates.
 *
 * @param barycentric_coordinates
 * @return vector3
 */
vector3 Tetra::compute_euclidean_coordinates(const std::array<double, 4>& barycentric_coordinates) const {
    return (
        barycentric_coordinates[0] * m_list_vertices[0]->get_position() + barycentric_coordinates[1] * m_list_vertices[1]->get_position() +
        barycentric_coordinates[2] * m_list_vertices[2]->get_position() + barycentric_coordinates[3] * m_list_vertices[3]->get_position());
}

}  // namespace bz_mesh