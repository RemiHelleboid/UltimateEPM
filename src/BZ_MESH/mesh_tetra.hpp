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

namespace uepm::mesh_bz {

using array4d = std::array<double, 4>;

struct IsoWeightsResult {
    double                 dos_eV{0.0};    // tranche de DOS tétra à Ef [1/eV]
    std::array<double, 4>  w{0, 0, 0, 0};  // poids par sommet, somme = dos_eV
    std::array<vector3, 3> tri{};          // (optionnel) points de l’iso-triangle
    vector3                centroid{};     // (optionnel) centroïde de l’iso-triangle
};

struct UniformDos {
    bool               valid{false};
    double             E0{0.0};
    double             Emax{0.0};
    double             inv_dx{0.0};  // 1 / dx
    uint32_t           N{0};         // number of knots = nb_steps + 1
    std::vector<float> D;            // try float to halve bandwidth; use double if needed

    inline double sample_or_zero(double E) const noexcept {
        if (!valid || E < E0 || E > Emax) {
            return 0.0;
        }
        const double t   = (E - E0) * inv_dx;
        int          idx = static_cast<int>(t);  // floor
        // clamp to [0, N-2]
        const int last = static_cast<int>(N) - 2;
        if (idx < 0) {
            idx = 0;
        }
        if (idx > last) {
            idx = last;
        }

        const float  d0   = D[static_cast<size_t>(idx)];
        const float  d1   = D[static_cast<size_t>(idx + 1)];
        const double frac = t - static_cast<double>(idx);
        // one FMA; cast once to double to keep precision in math path
        return std::fma(frac, static_cast<double>(d1 - d0), static_cast<double>(d0));
    }
};

class Tetra {
 private:
    /**
     * @brief Element index.
     *
     */
    std::size_t m_index;

    /**
     * @brief True if the tetrahedra lies in the irreducible wedge of the BZ.
     *
     */
    bool m_lies_in_iwedge{false};

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
     * @brief Barycenter of the tetrahedra.
     *
     */
    vector3 m_barycenter;

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
    std::vector<vector3> m_gradient_energy_per_band;

    /**
     * @brief For each band, store the indices of the vertices sorted by increasing energy.
     * For example, if for band k, vertex 2 has the lowest energy, then vertex 0, then vertex 3, then vertex 1,
     * then m_sorted_slots_per_band[k] = {2, 0, 3, 1}.
     *
     */
    std::vector<std::array<int, 4>> m_sorted_slots_per_band;

    /**
     * @brief Precomputed DOS on a uniform energy grid for each band.
     *
     */
    // std::vector<UniformDos> m_dos_per_band;

 public:
    /**
     * @brief There is not default constructor for Tetra class.
     *
     */
    Tetra()                        = default;
    Tetra(const Tetra&)            = default;
    Tetra& operator=(const Tetra&) = default;
    Tetra(Tetra&&)                 = default;
    Tetra& operator=(Tetra&&)      = default;

    Tetra(std::size_t index, const std::array<Vertex*, 4>& list_vertices);

    const bbox_mesh& get_bounding_box() const;
    bbox_mesh        compute_bounding_box() const;
    bool             lies_in_irreducible_wedge() const { return m_lies_in_iwedge; }
    void             set_lies_in_irreducible_wedge(bool value) { m_lies_in_iwedge = value; }

    void    compute_min_max_energies_at_bands();
    double  get_min_energy_at_band(std::size_t band_index) const { return m_min_energy_per_band[band_index]; }
    double  get_max_energy_at_band(std::size_t band_index) const { return m_max_energy_per_band[band_index]; }
    vector3 get_barycenter() const { return m_barycenter; }

    std::size_t                   get_index() const { return m_index; }
    const std::array<Vertex*, 4>& get_list_vertices() const { return m_list_vertices; }
    std::array<std::size_t, 4>    get_list_indices_vertices() const {
        return {m_list_vertices[0]->get_index(),
                   m_list_vertices[1]->get_index(),
                   m_list_vertices[2]->get_index(),
                   m_list_vertices[3]->get_index()};
    }
    const std::array<vector3, 6>& get_list_edges() const { return m_list_edges; }
    std::size_t                   get_nb_bands() const { return m_nb_bands; }

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

    void                      pre_compute_sorted_slots_per_band();
    const std::array<int, 4>& get_index_vertices_with_sorted_energy_at_band(std::size_t index_band) const {
        return m_sorted_slots_per_band[index_band];
    }

    void   precompute_dos_on_energy_grid_per_band(double energy_step, double energy_threshold);
    double interpolate_dos_at_energy_per_band(double energy, std::size_t band_index) const noexcept;

    const vector3& get_gradient_energy_at_band(std::size_t band_index) const { return m_gradient_energy_per_band[band_index]; }

    bool                 is_energy_inside_band(double energy, std::size_t index_band) const;
    bool                 does_intersect_band_energy_range(double e_min, double e_max, std::size_t index_band) const;
    std::vector<vector3> compute_band_iso_energy_surface(double iso_energy, std::size_t band_index) const;
    double               compute_tetra_iso_surface_energy_band(double energy, std::size_t band_index) const;
    double               compute_tetra_dos_energy_band(double energy, std::size_t band_index) const;
    vector3              draw_random_uniform_point_at_energy(double iso_energy, std::size_t band_index, std::mt19937& rng) const;

    std::array<double, 8> get_tetra_electron_phonon_rates(int band_index) const;
    std::array<double, 8> interpolate_phonon_scattering_rate_at_location(const vector3& location, const std::size_t& band_index) const;

    double interpolate_scalar_at_position(const vector3& location, const std::vector<double>& scalar_field) const;
    double interpolate_energy_at_band(const vector3& location, std::size_t band_index) const;

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

    void print_info() const {
        std::cout << "Tetra #" << m_index << ":\n";
        std::cout << "  Vertices indices: ";
        for (const auto& vtx : m_list_vertices) {
            std::cout << vtx->get_index() << std::endl;
            int nb_bands = vtx->get_number_bands();
            for (int i = 0; i < nb_bands; ++i) {
                std::cout << "    Band " << i << ": Energy = " << vtx->get_energy_at_band(i) << " eV\n";
            }
        }
        std::cout << "\n";
        std::cout << "  Signed volume: " << m_signed_volume << "\n";
        std::cout << "  Barycenter: " << m_barycenter << "\n";
    }
};

}  // namespace uepm::mesh_bz