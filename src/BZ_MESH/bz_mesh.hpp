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
#include <string>
#include <utility>
#include <vector>

#include "BandStructure.h"
#include "band_catalog.hpp"
#include "bz_domain.hpp"
#include "bz_dos.hpp"
#include "epm_material.hpp"
#include "export_octree_vtu.hpp"
#include "mesh_tetra.hpp"
#include "mesh_vertex.hpp"
#include "octree_bz.hpp"
#include "reciprocal_space.hpp"
#include "tetra_energy_index.hpp"
#include "vector_bz.hpp"

namespace uepm::mesh_bz {

using MapStringToDoubles = std::map<std::string, std::vector<double>>;
using MapStringToVectors = std::map<std::string, std::vector<vector3>>;

class MeshBZ {
 protected:
    void compact_geometry_to_positive_octant();
    void build_positive_octant_kstar();

    std::string m_filename_mesh;

    int m_nb_threads_mesh_ops = 1;

    uepm::pseudopotential::epm_material m_material;

    vector3 m_center{0.0, 0.0, 0.0};

    std::vector<std::size_t>              m_node_tags;
    std::vector<std::size_t>              m_source_vertex_indices;
    std::vector<std::size_t>              m_local_index_from_source_vertex;
    std::vector<Vertex>                   m_list_vertices;
    std::vector<Tetra>                    m_list_tetrahedra;
    std::vector<std::vector<std::size_t>> m_vertex_to_tetrahedra;

    // Corrects small discrepancies between the integrated full-BZ mesh volume
    // and the analytical reciprocal-cell volume.
    double m_bz_volume_correction = 1.0;

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
    TetraEnergyIndex m_tetra_energy_index;

    /**
     * @brief Octree search tree for fast spatial queries.
     *
     */
    std::unique_ptr<Octree_mesh> m_search_tree;

    BandCatalog m_bands;

    double m_max_energy_global = 1e100;

    double m_total_volume    = 0.0;
    double m_spin_degeneracy = 2.0;

    ReciprocalSpace m_reciprocal_space;
    BZDomainMode    m_domain_mode = BZDomainMode::full;

 public:
    // ---------- ctors/dtor ----------
    MeshBZ() = default;
    explicit MeshBZ(const uepm::pseudopotential::epm_material& material) : m_material(material) {}
    MeshBZ(const MeshBZ&)                = delete;
    MeshBZ& operator=(const MeshBZ&)     = delete;
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

    double       get_bz_volume_correction() const noexcept { return m_bz_volume_correction; }
    double       get_spin_degeneracy() const noexcept { return m_spin_degeneracy; }
    void         set_domain_mode(BZDomainMode mode) noexcept { m_domain_mode = mode; }
    BZDomainMode domain_mode() const noexcept { return m_domain_mode; }
    bool         stores_positive_octant() const noexcept { return m_domain_mode == BZDomainMode::positive_octant; }
    double       stored_domain_multiplicity() const noexcept { return stores_positive_octant() ? 8.0 : 1.0; }
    CanonicalK   canonicalize_physical_k(const vector3& physical_k) const noexcept {
        return canonicalize_k(physical_k, m_domain_mode);
    }
    vector3 representative_vector_to_physical(const vector3&            representative,
                                              const std::array<int, 3>& signs) const noexcept {
        return apply_sign_image(representative, signs);
    }
    const auto& physical_sign_images() const noexcept { return positive_octant_images; }
    std::size_t source_vertex_index(std::size_t local_index) const { return m_source_vertex_indices.at(local_index); }
    std::size_t local_vertex_index_from_source(std::size_t source_index) const;
    void        set_bz_volume_correction(double factor) noexcept { m_bz_volume_correction = factor; }
    // Compatibility aliases for existing callers.
    double  get_reduce_bz_factor() const noexcept { return get_bz_volume_correction(); }
    void    set_reduce_bz_factor(double factor) noexcept { set_bz_volume_correction(factor); }
    double  si_to_reduced_scale() const noexcept;
    vector3 si_to_reduced_k(const vector3& k_si) const noexcept;
    vector3 reduced_to_si_k(const vector3& k_reduced) const noexcept;

    std::size_t get_number_vertices() const noexcept { return m_list_vertices.size(); }
    std::size_t get_number_elements() const noexcept { return m_list_tetrahedra.size(); }
    double      get_volume() const noexcept { return m_total_volume; }

    const std::vector<Vertex>&      get_list_vertices() const noexcept { return m_list_vertices; }
    const std::vector<Tetra>&       get_list_tetrahedra() const noexcept { return m_list_tetrahedra; }
    const BandCatalog&              band_catalog() const noexcept { return m_bands; }
    const ReciprocalSpace&          reciprocal_space() const noexcept { return m_reciprocal_space; }
    const TetraEnergyIndex&         tetra_energy_index() const noexcept { return m_tetra_energy_index; }
    const std::vector<std::size_t>& get_list_vertex_indices_in_irreducible_wedge() const noexcept {
        return m_list_vtx_in_iwedge;
    }
    std::size_t get_number_vertices_in_irreducible_wedge() const noexcept { return m_list_vtx_in_iwedge.size(); }
    const std::vector<std::vector<std::size_t>>& get_kstar_ibz_to_bz() const noexcept { return m_kstar_ibz_to_bz; }
    std::size_t get_multiplicity_of_kpoint_in_ibz(std::size_t ibz_kpoint_index) const noexcept {
        return m_kstar_ibz_to_bz[ibz_kpoint_index].size();
    }

    // Band infos
    std::size_t get_number_bands_total() const noexcept { return m_bands.total(); }
    std::size_t get_number_valence_bands() const noexcept { return m_bands.range(MeshParticleType::valence).count; }
    std::size_t get_number_conduction_bands() const noexcept {
        return m_bands.range(MeshParticleType::conduction).count;
    }
    std::size_t get_number_bands(MeshParticleType type) const noexcept { return m_bands.range(type).count; }
    std::vector<std::size_t> get_band_indices(MeshParticleType type) const;
    std::pair<int, int>      get_start_end_valence_band_idx() const {
        const auto& range = m_bands.range(MeshParticleType::valence);
        return {range.global_start_index, range.global_start_index + range.count};
    }
    std::pair<int, int> get_start_end_conduction_band_idx() const {
        const auto& range = m_bands.range(MeshParticleType::conduction);
        return {range.global_start_index, range.global_start_index + range.count};
    }
    std::size_t get_local_band_index(int global_band_index) const { return m_bands.local_index(global_band_index); }
    std::size_t get_global_band_index(std::size_t local_band_index, MeshParticleType type) const {
        return m_bands.global_index(local_band_index, type);
    }

    void print_band_info() const;

    // ---------- geometry / search ----------
    bbox_mesh           compute_bounding_box() const;
    void                build_search_tree();
    std::vector<Tetra*> get_list_p_tetra() {
        std::vector<Tetra*> ptrs;
        std::transform(m_list_tetrahedra.begin(), m_list_tetrahedra.end(), std::back_inserter(ptrs), [](Tetra& t) {
            return &t;
        });
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
                                       const std::string& model_name_or_msh_path,
                                       bool               write_gradients) const;

    // ---------- reading ----------
    /**
     * @brief Read a mesh while storing all internal k coordinates in SI (1/m).
     *
     * @param input_coordinates_are_reduced Set true when the file uses units of
     *        2*pi/a, false when it already contains SI coordinates.
     */
    void read_mesh_geometry_from_msh_file(const std::string& filename, bool input_coordinates_are_reduced = true);
    void read_mesh_bands_from_msh_file(const std::string& filename,
                                       int                nb_conduction_bands        = -1,
                                       int                nb_valence_bands           = -1,
                                       bool               auto_shift_conduction_band = false,
                                       bool               set_positive_valence_band  = false);
    void append_band(MeshParticleType            type,
                     const std::vector<double>&  energies_at_vertices,
                     const std::vector<vector3>& gradients_at_vertices);
    void add_new_band_energies_to_vertices(const std::vector<double>& energies_at_vertices);
    void add_new_gradient_band_energies_to_vertices(const std::vector<double>& gradients_at_vertices);
    void keep_only_bands(std::size_t nb_valence_bands, std::size_t nb_conduction_bands);
    void load_kstar_ibz_to_bz(const std::string& kstarFilePath = "");

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
        return m_tetra_energy_index.at(band_index).ordered_tetra_indices;
    }
    std::span<const std::size_t> get_candidate_tetra_indices_at_band(std::size_t band_index,
                                                                     double      minimum_energy,
                                                                     double      maximum_energy) const noexcept {
        return m_tetra_energy_index.at(band_index).candidate_indices(minimum_energy, maximum_energy);
    }
    const BandTetraEnergyIndex& get_tetra_energy_index_at_band(std::size_t band_index) const {
        return m_tetra_energy_index.at(band_index);
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
    vector3     get_k_at_index(std::size_t index) const { return m_list_vertices[index].get_position(); }
    std::size_t get_nearest_k_index(const vector3& k) const;
    const std::vector<std::size_t>& get_tetrahedra_of_vertex(std::size_t vi) const {
        return m_vertex_to_tetrahedra[vi];
    }

    std::pair<double, double> get_min_max_energy_at_band(const int& band_index) const {
        return m_bands.extrema(static_cast<std::size_t>(band_index));
    }

    // ---------- metrics / DOS ----------
    double compute_mesh_volume() const;
    double compute_iso_surface(double iso_energy, int band_index) const;
    double compute_dos_at_energy_and_band(double iso_energy,
                                          int    band_index,
                                          bool   use_interp = false,
                                          bool   use_iw     = false) const;

    std::size_t draw_random_tetrahedron_index_with_dos_probability(double        energy,
                                                                   std::size_t   idx_band,
                                                                   std::mt19937& rng) const;

    vector3 draw_random_k_point_at_energy(double energy, std::size_t idx_band, std::mt19937& rng) const;
    std::pair<vector3, std::size_t> draw_random_k_point_at_energy(double energy, std::mt19937& rng) const;

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
