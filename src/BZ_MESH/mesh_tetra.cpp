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
#include <cmath>
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

vector3 Tetra::compute_barycenter() const {
    return (m_list_vertices[0]->get_position() + m_list_vertices[1]->get_position() + m_list_vertices[2]->get_position() +
            m_list_vertices[3]->get_position()) /
           4.0;
}

void Tetra::compute_gradient_energy_at_bands() {
    m_gradient_energy_per_band.clear();
    std::size_t m_nb_bands = m_list_vertices[0]->get_number_bands();
    m_gradient_energy_per_band.reserve(m_nb_bands);
    for (std::size_t band_index = 0; band_index < m_nb_bands; band_index++) {
        const std::array<double, 4> energies_at_vertices = get_band_energies_at_vertices(band_index);
        const std::array<int, 4>    indices_sort         = get_index_vertices_with_sorted_energy_at_band(band_index);
        const double                e_0                  = energies_at_vertices[indices_sort[0]];
        const double                eps_12               = (e_0 - energies_at_vertices[indices_sort[1]]);
        const double                eps_13               = (e_0 - energies_at_vertices[indices_sort[2]]);
        const double                eps_14               = (e_0 - energies_at_vertices[indices_sort[3]]);
        const double                gradient_energy      = sqrt((eps_12 * eps_12 + eps_13 * eps_13 + eps_14 * eps_14));
        if (gradient_energy == 0) {
            std::cout << "Warning: gradient energy is zero for tetra " << m_index << " at band " << band_index << std::endl;
        }
        m_gradient_energy_per_band.push_back(gradient_energy);
    }
}

/**
 * @brief Compute the minimum and maximum energies of the bands within the tetrahedron.
 *  *
 */
void Tetra::compute_min_max_energies_at_bands() {
    m_nb_bands = m_list_vertices[0]->get_number_bands();
    for (std::size_t idx_band = 0; idx_band < m_nb_bands; ++idx_band) {
        auto energies = get_band_energies_at_vertices(idx_band);
        auto minmax   = std::minmax_element(energies.begin(), energies.end());
        m_min_energy_per_band.push_back(*minmax.first);
        m_max_energy_per_band.push_back(*minmax.second);
    }
}

/**
 * @brief Compute the signed volume of the tetrahedron.
 *
 * @return double
 */
double Tetra::compute_signed_volume() const {
    // std::cout << m_list_edges[0] << std::endl;
    // std::cout << m_list_edges[1] << std::endl;
    return (1.0 / 6.0) * scalar_triple_product(m_list_edges[0], m_list_edges[1], m_list_edges[2]);
}

/**
 * @brief Return the values of the energy of the index_band valence band at the 4 vertices of the tetrahedra.
 *
 * @param index_band
 * @return std::vector<double>
 */
std::array<double, 4> Tetra::get_band_energies_at_vertices(std::size_t index_band) const {
    return {m_list_vertices[0]->get_energy_at_band(index_band),
            m_list_vertices[1]->get_energy_at_band(index_band),
            m_list_vertices[2]->get_energy_at_band(index_band),
            m_list_vertices[3]->get_energy_at_band(index_band)};
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

/**
 * @brief Compute the euclidean position from barycentric coordinates, with a given vertices order,
 * that might be different from the vertices of the tetrahedra.
 *
 * @param barycentric_coordinates
 * @param indices_vertex
 * @return vector3
 */
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
 * This function is written explicitely instead of using std::sort functions, because the sorting is done
 * with the minimum number of operations for a 4 values sorting. Other solution might be tested later.
 *
 * @param index_band
 * @return std::array<int, 4>
 */
std::array<int, 4> Tetra::get_index_vertices_with_sorted_energy_at_band(std::size_t index_band) const {
    std::array<double, 4> energies_at_vertices = get_band_energies_at_vertices(index_band);
    std::array<int, 4>    sorted_index         = {0, 1, 2, 3};
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

/**
 * @brief Compute the iso-energy surface within the tetrahedra for a given energy of a given band.
 * The surface is returned as a list of points (3 when the surface is a triangle, 4 when it is a quadrangle).
 *
 * The case of energy being smaller than the minimum energy of the tetrahedra is not taken into account.
 * Same thing for the case of energy being greater than the maximum energy of the tetrahedra.
 * Those two cases are handle by the caller function. This is done to avoid computing the sorted index which is computationallly intensive.
 * The minimum and maximum energies are stored in the member variables
 * m_min_energy_at_vertices and m_max_energy_at_vertices at the construction of the tetrahedra.
 *
 * This is very important because those 2 trivial cases represent usually more than 95% of the cases.
 *
 *
 * @param iso_energy
 * @param band_index
 * @return std::vector<vector3>
 */
std::vector<vector3> Tetra::compute_band_iso_energy_surface(double iso_energy, std::size_t band_index) const {
    std::array<double, 4> energies_at_vertices = get_band_energies_at_vertices(band_index);
    std::array<int, 4>    indices_sort         = get_index_vertices_with_sorted_energy_at_band(band_index);
    double                e_0                  = energies_at_vertices[indices_sort[0]];
    double                e_1                  = energies_at_vertices[indices_sort[1]];
    double                e_2                  = energies_at_vertices[indices_sort[2]];
    double                e_3                  = energies_at_vertices[indices_sort[3]];

    bool check_order = (e_0 <= e_1 && e_1 <= e_2 && e_2 <= e_3);
    if (!check_order) {
        std::cerr << "Error: the order of the energies is not correct" << std::endl;
        throw std::runtime_error("Error: the order of the energies is not correct");
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
        vector3 V    = compute_euclidean_coordinates_with_indices({1.0 - lA_V, 0.0, 0.0, lA_V}, indices_sort);
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

/**
 * @brief Compute the surface of the tetrahedra for a given energy of a given band.
 * The iso-surface is computed by the function compute_band_iso_energy_surface, and then
 * the area of the surface is computed.
 *
 * If the surface is a triangle, the area is computed by the class IsoTriangle class function "get_signed_area".
 * If the surface is a quadrangle, the area is computed by splitting the quadrangle into 2 triangles and then computing the area of each
 * triangle.
 *
 * A more direct way to compute the area in both cases should be tested for performances improvement.
 *  *
 * @param energy
 * @param band_index
 * @return double
 */
double Tetra::compute_tetra_iso_surface_energy_band(double energy, std::size_t band_index) const {
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

/**
 * @brief Main function to compute the DOS of the energy energy of the band with index band_index within the tetrahedra.
 *
 * One could try to store the value eps_12, eps_13, eps_14 for each band at the construction step, to avoid computing them each time the
 * function is called.
 *
 * @param energy
 * @param band_index
 * @return double
 */
double Tetra::compute_tetra_dos_energy_band(double energy, std::size_t band_index) const {
    if (energy < m_min_energy_per_band[band_index] || energy > m_max_energy_per_band[band_index]) {
        return 0.0;
    }
    const double renormalization = 6.0 * fabs(this->m_signed_volume);
    return renormalization * (1.0 / m_gradient_energy_per_band[band_index]) *
           this->compute_tetra_iso_surface_energy_band(energy, band_index);
}

}  // namespace bz_mesh