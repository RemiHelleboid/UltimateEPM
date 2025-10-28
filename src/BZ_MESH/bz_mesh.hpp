/**
 * @file bz_mesh.hpp
 * @author Rémi Helleboid (remi.helleboid@st.com)
 * @brief
 * @version 0.1
 * @date 2022-07-14
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <memory>
#include <random>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "BandStructure.h"
#include "Material.h"
#include "export_octree_vtu.hpp"
#include "mesh_tetra.hpp"
#include "mesh_vertex.hpp"
#include "octree_bz.hpp"
#include "vector.hpp"

namespace uepm::mesh_bz {

using MapStringToDoubles = std::map<std::string, std::vector<double>>;
using MapStringToVectors = std::map<std::string, std::vector<vector3>>;

enum class MeshParticleType { valence, conduction };

struct BandInfo {
    MeshParticleType type{MeshParticleType::conduction};
    std::size_t      local_index{0};
};

struct BandRange {
    int         global_start_index{-1};
    std::size_t count{0};
};

struct TetraOrderedEnergyMin {
    std::size_t              m_band_idx;
    std::vector<std::size_t> m_ordered_tetra_indices;
    std::vector<double>      m_ordered_energies;
    double                   m_max_energy_spread{0.0};
    double                   max_energy = 1e100;

    auto get_inidices_in_energy_range(double E_min, double E_max) const {
        auto        it_low   = std::lower_bound(m_ordered_energies.begin(), m_ordered_energies.end(), E_min);
        auto        it_high  = std::upper_bound(m_ordered_energies.begin(), m_ordered_energies.end(), E_max);
        std::size_t idx_low  = std::distance(m_ordered_energies.begin(), it_low);
        std::size_t idx_high = std::distance(m_ordered_energies.begin(), it_high);
        return std::make_pair(idx_low, idx_high);
    }

    std::span<const std::size_t> indices_in_min_energy_range(double E_min, double E_max) const noexcept {
        // binary-search on energies
        const auto        it_low  = std::lower_bound(m_ordered_energies.begin(), m_ordered_energies.end(), E_min);
        const auto        it_high = std::upper_bound(m_ordered_energies.begin(), m_ordered_energies.end(), E_max);
        const std::size_t lo      = static_cast<std::size_t>(it_low - m_ordered_energies.begin());
        const std::size_t hi      = static_cast<std::size_t>(it_high - m_ordered_energies.begin());

        // return a view over the corresponding slice of tetra indices
        const auto all = std::span<const std::size_t>(m_ordered_tetra_indices);
        return all.subspan(lo, hi - lo);
    }

    std::span<const std::size_t> candidate_indices(double E_min, double E_max) const noexcept {
        return indices_in_min_energy_range(E_min - m_max_energy_spread, E_max);
    }
};

class MeshBZ {
 protected:
    std::string m_filename_mesh;

    int m_nb_threads_mesh_ops = 1;

    uepm::pseudopotential::Material m_material;

    vector3 m_center{0.0, 0.0, 0.0};

    std::vector<std::size_t>              m_node_tags;
    std::vector<Vertex>                   m_list_vertices;
    std::vector<Tetra>                    m_list_tetrahedra;
    std::vector<std::vector<std::size_t>> m_vertex_to_tetrahedra;

    double m_reduce_bz_factor = 1.0;

    /**
     * @brief List of indices of vertices that lie in the irreducible wedge of the BZ.
     *
     */
    std::vector<std::size_t> m_list_vtx_in_iwedge;

    /**
     * @brief Mapping from k-points in the irreducible Brillouin zone to the full Brillouin zone.
     * m_kstar_ibz_to_bz[i] is the list of indices in the full BZ equivalent to the i-th k-point in the IBZ.
     *
     */
    std::vector<std::vector<std::size_t>> m_kstar_ibz_to_bz;

    /**
     * @brief For each band, list of tetrahedra ordered by increasing minimum energy in the tetrahedron.
     *
     */
    std::vector<TetraOrderedEnergyMin> m_tetra_ordered_energy_min;

    /**
     * @brief Octree search tree for fast spatial queries.
     *
     */
    std::unique_ptr<Octree_mesh> m_search_tree;

    std::size_t           m_nb_bands_total = 0;
    std::vector<BandInfo> m_band_info;
    BandRange             m_valence_bands;
    BandRange             m_conduction_bands;

    std::vector<double> m_min_band;
    std::vector<double> m_max_band;

    double m_max_energy_global = 1e100;

    double m_total_volume    = 0.0;
    double m_spin_degeneracy = 2.0;

    std::vector<vector3> m_Gshifts;

    // Wigner–Seitz / BZ folding helpers
    Eigen::Matrix3d m_recip_B      = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d m_recip_Bi     = Eigen::Matrix3d::Identity();
    double          m_si2red       = 1.0;
    double          m_bz_halfwidth = 1.0;

 public:
    // ---------- ctors/dtor ----------
    MeshBZ() = default;
    explicit MeshBZ(const uepm::pseudopotential::Material& material) : m_material(material) {}
    MeshBZ(const MeshBZ&)                = default;
    MeshBZ& operator=(const MeshBZ&)     = default;
    MeshBZ(MeshBZ&&) noexcept            = default;
    MeshBZ& operator=(MeshBZ&&) noexcept = default;
    // ~MeshBZ();  // out-of-line (needed for unique_ptr<Octree_mesh> with fwd-decl)

    // ---------- light getters / basics ----------
    int  get_number_threads_mesh_ops() const noexcept { return m_nb_threads_mesh_ops; }
    void set_number_threads_mesh_ops(int nb_threads) noexcept { m_nb_threads_mesh_ops = nb_threads; }

    double get_max_energy_global() const noexcept { return m_max_energy_global; }
    void   set_max_energy_global(double max_energy) noexcept { m_max_energy_global = max_energy; }

    const vector3& get_vertex_position(std::size_t idx_vtx) const { return m_list_vertices[idx_vtx].get_position(); }
    const vector3& get_center() const noexcept { return m_center; }
    void           shift_bz_center(const vector3& shift);

    double get_reduce_bz_factor() const noexcept { return m_reduce_bz_factor; }
    void   set_reduce_bz_factor(double factor) noexcept { m_reduce_bz_factor = factor; }
    double si_to_reduced_scale() const noexcept;

    std::size_t get_number_vertices() const noexcept { return m_list_vertices.size(); }
    std::size_t get_number_elements() const noexcept { return m_list_tetrahedra.size(); }
    double      get_volume() const noexcept { return m_total_volume; }

    const std::vector<Vertex>&      get_list_vertices() const noexcept { return m_list_vertices; }
    const std::vector<Tetra>&       get_list_tetrahedra() const noexcept { return m_list_tetrahedra; }
    const std::vector<std::size_t>& get_list_vertex_indices_in_irreducible_wedge() const noexcept { return m_list_vtx_in_iwedge; }
    std::size_t                     get_number_vertices_in_irreducible_wedge() const noexcept { return m_list_vtx_in_iwedge.size(); }
    const std::vector<std::vector<std::size_t>>& get_kstar_ibz_to_bz() const noexcept { return m_kstar_ibz_to_bz; }
    std::size_t                                  get_multiplicity_of_kpoint_in_ibz(std::size_t ibz_kpoint_index) const noexcept {
        return m_kstar_ibz_to_bz[ibz_kpoint_index].size();
    }

    // Band infos
    std::size_t get_number_bands_total() const noexcept { return m_nb_bands_total; }
    std::size_t get_number_valence_bands() const noexcept { return m_valence_bands.count; }
    std::size_t get_number_conduction_bands() const noexcept { return m_conduction_bands.count; }
    std::size_t get_number_bands(MeshParticleType type) const noexcept {
        return (type == MeshParticleType::valence) ? m_valence_bands.count : m_conduction_bands.count;
    }
    std::vector<std::size_t> get_band_indices(MeshParticleType type) const;
    std::pair<int, int>      get_start_end_valence_band_idx() const {
        return {m_valence_bands.global_start_index, m_valence_bands.global_start_index + m_valence_bands.count};
    }
    std::pair<int, int> get_start_end_conduction_band_idx() const {
        return {m_conduction_bands.global_start_index, m_conduction_bands.global_start_index + m_conduction_bands.count};
    }
    std::size_t get_local_band_index(int global_band_index) const {
        if (global_band_index < 0 || global_band_index >= static_cast<int>(m_nb_bands_total)) {
            throw std::out_of_range("Global band index out of range.");
        }
        return m_band_info[global_band_index].local_index;
    }
    std::size_t get_global_band_index(std::size_t local_band_index, MeshParticleType type) const {
        if (type == MeshParticleType::valence) {
            if (local_band_index >= m_valence_bands.count) {
                throw std::out_of_range("Local valence band index out of range.");
            }
            return m_valence_bands.global_start_index + local_band_index;
        } else {
            if (local_band_index >= m_conduction_bands.count) {
                throw std::out_of_range("Local conduction band index out of range.");
            }
            return m_conduction_bands.global_start_index + local_band_index;
        }
    }

    void print_band_info() const;

    // ---------- geometry / search ----------
    bbox_mesh           compute_bounding_box() const;
    void                build_search_tree();
    std::vector<Tetra*> get_list_p_tetra() {
        std::vector<Tetra*> ptrs;
        std::transform(m_list_tetrahedra.begin(), m_list_tetrahedra.end(), std::back_inserter(ptrs), [](Tetra& t) { return &t; });
        return ptrs;
    }
    Tetra* find_tetra_at_location(const vector3& location) const;

    // ---------- export / I/O (heavy headers live in .cpp) ----------
    void export_k_points_to_file(const std::string& filename) const;
    void export_to_vtk(const std::string&        filename,
                       const MapStringToDoubles& point_scalars = {},
                       const MapStringToVectors& point_vectors = {},
                       const MapStringToDoubles& cell_scalars  = {},
                       const MapStringToVectors& cell_vectors  = {}) const;
    void export_energies_and_gradients_to_vtk(const std::string& filename) const;
    void export_octree_to_vtu(const std::string& filename) const;
    void export_selected_bands_to_gmsh(const std::string& out_filename,
                                       std::size_t        nb_valence_to_export,
                                       std::size_t        nb_conduction_to_export,
                                       bool               highest_valence_as_band0,
                                       const std::string& model_name_or_msh_path) const;

    // ---------- reading ----------
    void read_mesh_geometry_from_msh_file(const std::string& filename, bool normalize_by_fourier_factor = true);
    void read_mesh_bands_from_msh_file(const std::string& filename,
                                       int                nb_conduction_bands        = -1,
                                       int                nb_valence_bands           = -1,
                                       bool               auto_shift_conduction_band = false,
                                       bool               set_positive_valence_band  = false);
    void add_new_band_energies_to_vertices(const std::vector<double>& energies_at_vertices);
    void keep_only_bands(std::size_t nb_valence_bands, std::size_t nb_conduction_bands);
    void load_kstar_ibz_to_bz(const std::string& filename = "kstar_ibz_to_bz.txt");

    // ---------- Band structure ----------
    void compute_band_structure_over_mesh(uepm::pseudopotential::BandStructure& band_structure, bool use_iwedge = true);
    void distribute_energies_from_iw_wedge_to_full_bz();

    // ---------- analysis / precompute ----------
    void apply_scissor(double scissor_value);
    void compute_min_max_energies_at_tetras();
    void compute_energy_gradient_at_tetras();
    void auto_shift_conduction_band_energies();
    void auto_set_positive_valence_band_energies();
    void set_bands_in_right_order();
    void recompute_min_max_energies();
    void recompute_tetra_ordered_energies(double max_energy = 1e100);

    /**
     * @brief Get the list of tetrahedron indices ordered by increasing minimum energy at given band,
     * limited to tetrahedra with minimum energy below m_max_energy_global.
     *
     * @param band_index
     * @return const std::vector<std::size_t>&
     */
    const std::vector<std::size_t>& get_ordered_tetra_indices_at_band(std::size_t band_index) const {
        return m_tetra_ordered_energy_min[band_index].m_ordered_tetra_indices;
    }

    void precompute_dos_tetra(double energy_step = 0.01, double energy_max = 100.0);
    void set_energy_gradient_at_vertices_by_averaging_tetras();

    void recompute_energies_data_and_sync(bool   recompute_min_max = true,
                                          bool   recompute_grad    = true,
                                          bool   recompute_dos     = true,
                                          double dos_energy_step   = 0.01,
                                          double dos_energy_max    = 100.0);
    // ---------- queries ----------
    vector3 interpolate_energy_gradient_at_location(const vector3& location, const std::size_t& idx_band) const;

    void    precompute_G_shifts();
    bool    is_inside_mesh_geometry(const vector3& k) const;
    vector3 retrieve_k_inside_mesh_geometry(const vector3& k) const;

    void init_reciprocal_basis(const Eigen::Vector3d& b1_SI,
                               const Eigen::Vector3d& b2_SI,
                               const Eigen::Vector3d& b3_SI,
                               double                 halfwidth_reduced,
                               double                 si_to_reduced);

    vector3                         fold_ws_bcc(const vector3& k_SI) const noexcept;
    bool                            inside_ws_bcc(const vector3& k_SI) const noexcept;
    bool                            is_irreducible_wedge(const vector3& k_SI) const noexcept;
    std::size_t                     get_index_irreducible_wedge(const vector3& k_SI) const;
    const std::vector<std::size_t>& get_all_equivalent_indices_in_bz(const vector3& k_SI) const noexcept {
        return m_kstar_ibz_to_bz[get_index_irreducible_wedge(k_SI)];
    }
    vector3                         get_k_at_index(std::size_t index) const { return m_list_vertices[index].get_position(); }
    std::size_t                     get_nearest_k_index(const vector3& k) const;
    const std::vector<std::size_t>& get_tetrahedra_of_vertex(std::size_t vi) const { return m_vertex_to_tetrahedra[vi]; }

    std::pair<double, double> get_min_max_energy_at_band(const int& band_index) const {
        return {m_min_band[band_index], m_max_band[band_index]};
    }

    // ---------- metrics / DOS ----------
    double compute_mesh_volume() const;
    double compute_iso_surface(double iso_energy, int band_index) const;
    double compute_dos_at_energy_and_band(double iso_energy, int band_index, bool use_interp = false, bool use_iw = false) const;

    std::size_t draw_random_tetrahedron_index_with_dos_probability(double energy, std::size_t idx_band, std::mt19937& rng) const;

    vector3 draw_random_k_point_at_energy(double energy, std::size_t idx_band, std::mt19937& rng) const;

    std::vector<std::vector<double>> compute_dos_band_at_band(int         band_index,
                                                              double      min_energy,
                                                              double      max_energy,
                                                              std::size_t nb_points,
                                                              bool        use_interp = false,
                                                              bool        use_iw     = false) const;

    std::vector<std::vector<double>> compute_dos_band_at_band_auto(int         band_index,
                                                                   std::size_t nb_points,
                                                                   bool        use_interp = false,
                                                                   bool        use_iw     = false) const;
};

}  // namespace uepm::mesh_bz
