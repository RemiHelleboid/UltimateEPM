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

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <filesystem>
#include <regex>
#include <sstream>
#include <unordered_map>

#include "Material.h"
#include "iso_triangle.hpp"
#include "mesh_tetra.hpp"
#include "mesh_vertex.hpp"
#include "octree_bz.hpp"
#include "vector.hpp"

namespace bz_mesh {


using EigenIntSparseMatrix = Eigen::SparseMatrix<int, Eigen::RowMajor>;
using TripletInt           = Eigen::Triplet<int>;

enum class BandType { valence, conduction };

class MeshBZ {
 protected:
    /**
     * @brief The material of the Brillouin zone.
     *
     */
    EmpiricalPseudopotential::Material m_material;

    /**
     * @brief G vector on which the BZ is centered.
     *
     */
    vector3 m_center = {0.0, 0.0, 0.0};

    /**
     * @brief The tags of the nodes in the mesh.
     * These tags are used to identify the vertices in the GMSH files when load is done.
     */
    std::vector<std::size_t> m_node_tags;

    /**
     * @brief List of the vertices of the BZ mesh. Each vertices represent a vector k within the Brillouin Zone.
     *
     */
    std::vector<Vertex> m_list_vertices;

    /**
     * @brief List of the elements (tetrahedra) of the mesh.
     *
     */
    std::vector<Tetra> m_list_tetrahedra;

    /**
     * @brief Stores for each vertex index the list of tetrahedra indices that contains this vertex.
     * For example, m_vertex_to_tetrahedra[5] = {0, 2, 3} means that the vertex with index 5 is part of the tetrahedra with indices 0, 2
     * and 3.
     *
     * This is used to speed up the search of the tetrahedra containing a given vertex.
     */
    std::vector<std::vector<std::size_t>> m_vertex_to_tetrahedra;

    /**
     * @brief Factor applied to each tetra/vertex when we use only a part of the 1st BZ.
     * For example m_reduced_facotr = 48 when we use the irreducible wedge.
     *
     */
    double m_reduce_bz_factor = 1.0;

    /**
     * @brief Octree used to search for the tetrahedra that are overlapping with a given point.
     *
     */
    std::unique_ptr<Octree_mesh> m_search_tree;

    /**
     * @brief The indexes of the valence bands within the Vertex list of energies.
     * For example, if indices_valence_bands = [0, 1, 2] it means that
     * the m_band_energies[0], m_band_energies[1], m_band_energies[2] correspond to valence bands.
     *
     */
    std::vector<int> m_indices_valence_bands{};

    /**
     * @brief The indexes of the conduction bands within the Vertex list of energies.
     * For example, if indices_conduction_bands = [3, 4, 5] it means that
     * the values m_band_energies[3], m_band_energies[4], m_band_energies[5] correspond to conduction bands.
     *
     */
    std::vector<int> m_indices_conduction_bands{};

    /**
     * @brief The value m_min_band[i] is the minimum energy of the band with index i.
     *
     */
    std::vector<double> m_min_band{};

    /**
     * @brief The value m_max_band[i] is the maximum energy of the band with index i.
     *
     */
    std::vector<double> m_max_band{};

    /**
     * @brief Total volume of the BZ mesh in m^3.
     *
     */
    double m_total_volume = 0.0;

    /**
     * @brief Spin degeneracy factor (2 when spin-orbit coupling is not considered).
     *
     */
    double m_spin_degeneracy = 2.0;

    /**
     * @brief Possible G vectors to fold back k vectors within the first BZ.
     *
     */
    std::vector<vector3> m_Gshifts;

    // --- O(1) WS-BZ folding (internal state) ---
    /**
     * @brief Reciprocal primitive basis (columns) in SI [1/m] and its inverse.
     */
    Eigen::Matrix3d m_recip_B  = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d m_recip_Bi = Eigen::Matrix3d::Identity();

    /**
     * @brief Scale from SI k (1/m) to the reduced frame used by plane tests.
     * Keep this consistent with is_inside_mesh_geometry().
     */
    double m_si2red = 1.0;

    /**
     * @brief Half-width in the reduced frame (e.g. 0.5 for [-0.5,0.5]).
     */
    double m_bz_halfwidth = 0.5;



 public:
    MeshBZ() = default;
    MeshBZ(const EmpiricalPseudopotential::Material& material) : m_material(material) {};
    MeshBZ(const MeshBZ&) = default;

    void export_k_points_to_file(const std::string& filename) const;

    vector3 get_vertex_position(std::size_t idx_vtx) const { return m_list_vertices[idx_vtx].get_position(); }

    vector3 get_center() const { return m_center; }
    void    shift_bz_center(const vector3& shift);

    double        get_reduce_bz_factor() const { return m_reduce_bz_factor; }
    void          set_reduce_bz_factor(double factor) { m_reduce_bz_factor = factor; }
    inline double si_to_reduced_scale() const;

    bbox_mesh           compute_bounding_box() const;
    void                build_search_tree();
    std::vector<Tetra*> get_list_p_tetra() {
        std::vector<Tetra*> tetra_pointers;
        tetra_pointers.reserve(m_list_tetrahedra.size());
        std::transform(m_list_tetrahedra.begin(), m_list_tetrahedra.end(), std::back_inserter(tetra_pointers), [](Tetra& tetra) {
            return &tetra;
        });
        return tetra_pointers;
    }
    Tetra* find_tetra_at_location(const vector3& location) const;

    void read_mesh_geometry_from_msh_file(const std::string& filename, bool normalize_by_fourier_factor = true);
    void read_mesh_bands_from_msh_file(const std::string& filename, int nb_bands_to_load = -1);
    void read_mesh_bands_from_multi_band_files(const std::string& dir_bands, int nb_bands_to_load = 100);
    void add_new_band_energies_to_vertices(const std::vector<double>& energies_at_vertices);
    void keep_only_bands(const int nb_valence_bands, const int nb_conduction_bands);
    void compute_min_max_energies_at_tetras();
    void compute_energy_gradient_at_tetras();
    void auto_shift_conduction_band_energies();
    void auto_set_positive_valence_band_energies();
    void set_bands_in_right_order();
    void recompute_min_max_energies();
    void precompute_dos_tetra(double energy_step = 0.01, double energy_threshold = 100.0);

    std::size_t      get_number_vertices() const { return m_list_vertices.size(); }
    std::size_t      get_number_elements() const { return m_list_tetrahedra.size(); }
    double           get_volume() const { return m_total_volume; }
    std::vector<int> get_indices_valence_bands() const { return m_indices_valence_bands; }
    std::vector<int> get_indices_conduction_bands() const { return m_indices_conduction_bands; }
    int              get_nb_bands() const { return m_indices_valence_bands.size() + m_indices_conduction_bands.size(); }

    vector3 interpolate_energy_gradient_at_location(const vector3& location, const std::size_t& idx_band) const;

    void                     precompute_G_shifts();
    bool                     is_inside_mesh_geometry(const vector3& k) const;
    vector3                  retrieve_k_inside_mesh_geometry(const vector3& k) const;
    void                     init_reciprocal_basis(const Eigen::Vector3d& b1_SI,
                                                   const Eigen::Vector3d& b2_SI,
                                                   const Eigen::Vector3d& b3_SI,
                                                   double                 halfwidth_reduced,
                                                   double                 si_to_reduced);
    vector3                  fold_ws_bcc(const vector3& k_SI) const noexcept;
    bool                     inside_ws_bcc(const vector3& k_SI) const noexcept;
    bool                     is_irreducible_wedge(const vector3& k_SI) const noexcept;
    std::size_t              get_index_irreducible_wedge(const vector3& k_SI) const noexcept;
    std::vector<std::size_t> get_all_equivalent_indices_in_bz(const vector3& k_SI) const noexcept;

    vector3                         get_k_at_index(std::size_t index) const { return m_list_vertices[index].get_position(); }
    std::size_t                     get_nearest_k_index(const vector3& k) const;
    const std::vector<std::size_t>& get_tetrahedra_of_vertex(std::size_t vertex_index) const {
        return m_vertex_to_tetrahedra[vertex_index];
    }

    std::size_t               get_number_bands() const { return m_min_band.size(); }
    std::pair<double, double> get_min_max_energy_at_band(const int& band_index) const {
        return std::make_pair(m_min_band[band_index], m_max_band[band_index]);
    }

    const std::vector<Vertex>& get_list_vertices() const { return m_list_vertices; }
    const std::vector<Tetra>&  get_list_elements() const { return m_list_tetrahedra; }

    double      compute_mesh_volume() const;
    double      compute_iso_surface(double iso_energy, int band_index) const;
    double      compute_dos_at_energy_and_band(double iso_energy, int band_index, bool use_interp = false) const;
    std::size_t draw_random_tetrahedron_index_with_dos_probability(double        energy,
                                                                   std::size_t   idx_band,
                                                                   std::mt19937& random_generator) const;

    vector3 draw_random_k_point_at_energy(double energy, std::size_t idx_band, std::mt19937& random_generator) const;

    std::vector<std::vector<double>> compute_dos_band_at_band(int         band_index,
                                                              double      min_energy,
                                                              double      max_energy,
                                                              int         num_threads,
                                                              std::size_t nb_points,
                                                              bool        use_interp = false) const;

    std::vector<std::vector<double>> compute_dos_band_at_band_auto(int         band_index,
                                                                   std::size_t nb_points,
                                                                   int         num_threads,
                                                                   bool        use_interp = false) const;


};

}  // namespace bz_mesh