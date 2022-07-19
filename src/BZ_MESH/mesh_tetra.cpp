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

#include "iso_triangle.hpp"

namespace bz_mesh {

std::vector<double> Tetra::ms_case_stats = {0, 0, 0, 0, 0};


/**
 * @brief Construct a new Tetra by passing directly the array of the four pointers to the vertices.
 *
 * @param list_vertices
 */
Tetra::Tetra(std::size_t index, const std::array<Vertex*, 4>& list_vertices)
    : m_index(index),
      m_list_vertices(list_vertices),
      m_nb_bands{m_list_vertices[0]->get_number_bands()} {
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
std::vector<double> Tetra::get_band_energies_at_vertices(std::size_t index_band) const {
    std::vector<double> list_energies_valence;
    for (auto& p_vtx : m_list_vertices) {
        list_energies_valence.push_back(p_vtx->get_energy_at_band(index_band));
    }
    return list_energies_valence;
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

vector3 Tetra::compute_euclidean_coordinates_with_indices(const std::array<double, 4>& barycentric_coordinates,
                                                          const std::array<int, 4>&    indices_vertex) const {
    return (barycentric_coordinates[0] * m_list_vertices[indices_vertex[0]]->get_position() +
            barycentric_coordinates[1] * m_list_vertices[indices_vertex[1]]->get_position() +
            barycentric_coordinates[2] * m_list_vertices[indices_vertex[2]]->get_position() +
            barycentric_coordinates[3] * m_list_vertices[indices_vertex[3]]->get_position());
}

/**
 * @brief Return a list of indices a, b, c, d such as, for the conduction band with index index_band,
 * we have Vtx_a <= Vtx_b <= Vtx_c <= Vtx_d in term of energy.
 *
 * @param index_band
 * @return std::array<int, 4>
 */
std::array<int, 4> Tetra::get_index_vertices_with_sorted_energy_at_band(std::size_t index_band) const {
    std::vector<double> energies_at_vertices = get_band_energies_at_vertices(index_band);
    std::array<int, 4>  sorted_index         = {0, 1, 2, 3};
    if (energies_at_vertices[0] > energies_at_vertices[1]) {
        std::swap(energies_at_vertices[0], energies_at_vertices[1]);
        std::swap(sorted_index[0], sorted_index[1]);
    }
    if (energies_at_vertices[2] > energies_at_vertices[3]) {
        std::swap(energies_at_vertices[2], energies_at_vertices[3]);
        std::swap(sorted_index[2], sorted_index[3]);
    }
    if (energies_at_vertices[0] > energies_at_vertices[2]) {
        std::swap(energies_at_vertices[0], energies_at_vertices[2]);
        std::swap(sorted_index[0], sorted_index[2]);
    }
    if (energies_at_vertices[1] > energies_at_vertices[3]) {
        std::swap(energies_at_vertices[1], energies_at_vertices[3]);
        std::swap(sorted_index[1], sorted_index[3]);
    }
    if (energies_at_vertices[1] > energies_at_vertices[2]) {
        std::swap(energies_at_vertices[1], energies_at_vertices[2]);
        std::swap(sorted_index[1], sorted_index[2]);
    }
    return sorted_index;
}

std::vector<vector3> Tetra::compute_band_iso_energy_surface(double iso_energy, std::size_t band_index) const {
    std::vector<double> energies_at_vertices = get_band_energies_at_vertices(band_index);
    std::array<int, 4>  indices_sort         = get_index_vertices_with_sorted_energy_at_band(band_index);
    double              e_0                  = energies_at_vertices[indices_sort[0]];
    double              e_1                  = energies_at_vertices[indices_sort[1]];
    double              e_2                  = energies_at_vertices[indices_sort[2]];
    double              e_3                  = energies_at_vertices[indices_sort[3]];

    bool check_order = (e_0 <= e_1 && e_1 <= e_2 && e_2 <= e_3);
    std::vector<vector3> list_points_iso_surface{};

    if (e_0 >= iso_energy) {
        // std::cout << "Case 1 " << e_0 << "\n";
        ms_case_stats[0]++;
        return {};
    }
    if (e_3 <= iso_energy) {
        // std::cout << "Case 2 " << e_3 << "\n";
        ms_case_stats[1]++;
        return {};
    }
    if (iso_energy < e_1 && iso_energy >= e_0) {
        ms_case_stats[2]++;
        double  lA_U = (iso_energy - e_0) / (e_1 - e_0);
        vector3 U    = compute_euclidean_coordinates_with_indices({1.0 - lA_U, lA_U, 0.0, 0.0}, indices_sort);
        double  lA_V = (iso_energy - e_0) / (e_2 - e_0);
        vector3 V    = compute_euclidean_coordinates_with_indices({1.0 - lA_V, 0.0, lA_V, 0.0}, indices_sort);
        double  lA_W = (iso_energy - e_0) / (e_3 - e_0);
        vector3 W    = compute_euclidean_coordinates_with_indices({1.0 - lA_W, 0.0, 0.0, lA_W}, indices_sort);
        return {U, V, W};
    }
    if (iso_energy < e_2 && iso_energy >= e_1) {
        ms_case_stats[3]++;
        double  lA_U = (iso_energy - e_0) / (e_2 - e_0);
        vector3 U    = compute_euclidean_coordinates_with_indices({1.0 - lA_U, 0.0, lA_U, 0.0}, indices_sort);
        double  lA_V = (iso_energy - e_0) / (e_3 - e_0);
        vector3 V    = compute_euclidean_coordinates_with_indices({1.0 - lA_V, 0.0, 0.0, lA_U}, indices_sort);
        double  lA_W = (e_2 - iso_energy) / (e_2 - e_1);
        vector3 W    = compute_euclidean_coordinates_with_indices({0.0, lA_W, 1.0 - lA_W, 0.0}, indices_sort);
        double  lA_X = (iso_energy - e_1) / (e_3 - e_1);
        vector3 X    = compute_euclidean_coordinates_with_indices({0.0, 1.0 - lA_X, 0.0, lA_X}, indices_sort);
        return {U, V, W, X};
    }
    if (iso_energy >= e_2) {
        ms_case_stats[4]++;
        double  lC_U = (e_3 - iso_energy) / (e_3 - e_2);
        vector3 U    = compute_euclidean_coordinates_with_indices({0.0, 0.0, lC_U, 1.0 - lC_U}, indices_sort);
        double  lB_V = (e_3 - iso_energy) / (e_3 - e_1);
        vector3 V    = compute_euclidean_coordinates_with_indices({0.0, lB_V, 0.0, 1.0 - lB_V}, indices_sort);
        double  lA_W = (e_3 - iso_energy) / (e_3 - e_0);
        vector3 W    = compute_euclidean_coordinates_with_indices({lA_W, 0.0, 1.0 - lA_W}, indices_sort);
        return {U, V, W};
    } else {
        throw std::runtime_error("ISO SURFACE CASE UNKNOWN IN DOS COMPUTATION... ABORT.");
    }
    return {};
}

double Tetra::compute_tetra_dos_band(double energy, std::size_t band_index) const {
    std::vector<vector3> vertices_iso_surface = compute_band_iso_energy_surface(energy, band_index);
    if (vertices_iso_surface.empty()) {
        return 0.0;
    } else if (vertices_iso_surface.size() == 3) {
        IsoTriangle triangle(vertices_iso_surface[0], vertices_iso_surface[1], vertices_iso_surface[2], energy);
        return fabs(triangle.get_signed_surface());
    } else {
        IsoTriangle triangle1(vertices_iso_surface[0], vertices_iso_surface[1], vertices_iso_surface[3], energy);
        IsoTriangle triangle2(vertices_iso_surface[0], vertices_iso_surface[1], vertices_iso_surface[2], energy);
        return fabs(triangle1.get_signed_surface()) + fabs(triangle2.get_signed_surface());
    }
}

}  // namespace bz_mesh