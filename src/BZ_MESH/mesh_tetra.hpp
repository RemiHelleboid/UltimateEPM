/**
 * @file mesh_tetra.hpp
 * @author Rémi Helleboid (remi.helleboid@st.com)
 * @brief
 * @version 0.1
 * @date 2022-07-14
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <array>
#include <complex>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "bbox_mesh.hpp"
#include "mesh_vertex.hpp"

namespace bz_mesh {

using array4d = std::array<double, 4>;

class Tetra {
 private:
    /**
     * @brief Element index.
     *
     */
    std::size_t m_index;

    /**
     * @brief The 4 vertices of the tetrahedra are stored as an array of pointers on Vertices object.
     *
     */
    std::array<Vertex*, 4> m_list_vertices{nullptr, nullptr, nullptr, nullptr};

    /**
     * @brief List of the edges vectors of the tetrahedra.
     * Stored as follows: [v01, v02, v03, v12, v13, v23]
     * where for example v13 = m_list_vertices[1] - m_list_vertices[3].
     *
     */
    std::array<vector3, 6> m_list_edges{};

    /**
     * @brief Bounding box of the tetrahedra.
     *
     */
    bbox_mesh m_bbox;

    /**
     * @brief Signed volume of the tetrahedra.
     * The sign depends on the "orientation" of the tetrahedra.
     *
     */
    double m_signed_volume = 0.0;

    /**
     * @brief Number of conduction bands.
     *
     */
    std::size_t m_nb_bands = 0;

    /**
     * @brief The value min_energy_per_band[k] is the minimum energy of the k-th band on the tetrahedra.
     * It is pre-computed and stored for optimization purposes.
     *
     */
    std::vector<double> m_min_energy_per_band;

    /**
     * @brief The value max_energy_per_band[k] is the maximum energy of the k-th band on the tetrahedra.
     * It is pre-computed and stored for optimization purposes.
     *
     */
    std::vector<double> m_max_energy_per_band;

    /**
     * @brief G>radient of the energy at the vertices of the tetrahedra for each band.
     * It is pre-computed and stored for optimization purposes.
     *
     */
    std::vector<double> m_gradient_energy_per_band;

 public:
    /**
     * @brief There is not default constructor for Tetra class.
     *
     */
    Tetra() = delete;

    Tetra(std::size_t index, const std::array<Vertex*, 4>& list_vertices);

    const bbox_mesh& get_bounding_box() const;
    bbox_mesh        compute_bounding_box() const;
    void             compute_min_max_energies_at_bands();

    std::size_t                   get_index() const { return m_index; }
    const std::array<Vertex*, 4>& get_list_vertices() const { return m_list_vertices; }
    std::array<std::size_t, 4>    get_list_indices_vertices() const {
        return {m_list_vertices[0]->get_index(),
                   m_list_vertices[1]->get_index(),
                   m_list_vertices[2]->get_index(),
                   m_list_vertices[3]->get_index()};
    }
    std::array<vector3, 6> get_list_edges() const { return m_list_edges; }
    std::size_t            get_nb_bands() const { return m_nb_bands; }

    std::array<double, 4> get_band_energies_at_vertices(std::size_t index_band) const;

    double  compute_signed_volume() const;
    double  get_signed_volume() const { return m_signed_volume; }
    vector3 compute_edge(std::size_t index_vtx_1, std::size_t index_vtx_2) const;
    void    compute_gradient_energy_at_bands();
    vector3 compute_gradient_at_tetra(const array4d& values_at_vertices) const;

    vector3               compute_barycenter() const;
    bool                  is_location_inside(const vector3& location) const;
    std::array<double, 4> compute_barycentric_coordinates(const vector3& location) const;
    vector3               compute_euclidean_coordinates(const std::array<double, 4>& barycentric_coordinates) const;
    vector3               compute_euclidean_coordinates_with_indices(const std::array<double, 4>& barycentric_coordinates,
                                                                     const std::array<int, 4>&    indices_vertex) const;

    bool                 is_energy_inside_band(double energy, std::size_t index_band) const;
    bool                 does_intersect_band_energy_range(double e_min, double e_max, std::size_t index_band) const;
    std::array<int, 4>   get_index_vertices_with_sorted_energy_at_band(std::size_t index_band) const;
    std::vector<vector3> compute_band_iso_energy_surface(double iso_energy, std::size_t band_index) const;
    double               compute_tetra_iso_surface_energy_band(double energy, std::size_t band_index) const;
    double               compute_tetra_dos_energy_band(double energy, std::size_t band_index) const;

    std::array<double, 8> get_tetra_electron_phonon_rates(int band_index) const;

    double interpolate_scalar_at_position(const std::array<double, 4>& barycentric_coordinates,
                                          const std::vector<double>&   scalar_field) const;

    vector3 interpolate_vector_at_position(const std::array<double, 4>& barycentric_coordinates,
                                           const std::vector<vector3>&  vector_field) const;
    template <typename T>
    T interpolate_at_position(const std::array<double, 4>& barycentric_coordinates, const std::vector<T>& field) const {
        T interpolated_value = T::Zero();
        for (std::size_t idx_vtx = 0; idx_vtx < 4; ++idx_vtx) {
            interpolated_value += barycentric_coordinates[idx_vtx] * field[idx_vtx];
        }
        return interpolated_value;
    }
    std::complex<double> interpolate_at_position(const std::array<double, 4>&             barycentric_coordinates,
                                                 const std::vector<std::complex<double>>& field) const {
        std::complex<double> interpolated_value = 0.0;
        for (std::size_t idx_vtx = 0; idx_vtx < 4; ++idx_vtx) {
            interpolated_value += barycentric_coordinates[idx_vtx] * field[idx_vtx];
        }

        return interpolated_value;
    }
};

}  // namespace bz_mesh