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

    double m_total_volume = 0.0;

 public:
    MeshBZ() = default;
    MeshBZ(const EmpiricalPseudopotential::Material& material) : m_material(material) {};
    MeshBZ(const MeshBZ&) = default;

    vector3 get_vertex_position(std::size_t idx_vtx) const { return m_list_vertices[idx_vtx].get_position(); }

    vector3 get_center() const { return m_center; }
    void    shift_bz_center(const vector3& shift);

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
    void read_mesh_bands_from_msh_file(const std::string& filename);
    void read_mesh_bands_from_multi_band_files(const std::string& dir_bands);
    void add_new_band_energies_to_vertices(const std::vector<double>& energies_at_vertices);
    void compute_min_max_energies_at_tetras();
    void compute_energy_gradient_at_tetras();

    
    std::size_t get_number_vertices() const { return m_list_vertices.size(); }
    std::size_t get_number_elements() const { return m_list_tetrahedra.size(); }
    double      get_volume() const { return m_total_volume; }

    bool    is_inside_mesh_geometry(const vector3& k) const;
    bool    is_inside_mesh_geometry(const Vector3D<double>& k) const;
    vector3 retrieve_k_inside_mesh_geometry(const vector3& k) const;

    vector3     get_k_at_index(std::size_t index) const { return m_list_vertices[index].get_position(); }
    std::size_t get_nearest_k_index(const Vector3D<double>& k) const;
    std::size_t get_nearest_k_index(const vector3& k) const;

    std::size_t               get_number_bands() const { return m_min_band.size(); }
    std::pair<double, double> get_min_max_energy_at_band(const int& band_index) const {
        return std::make_pair(m_min_band[band_index], m_max_band[band_index]);
    }

    const std::vector<Vertex>& get_list_vertices() const { return m_list_vertices; }
    const std::vector<Tetra>&  get_list_elements() const { return m_list_tetrahedra; }

    double compute_mesh_volume() const;
    double compute_iso_surface(double iso_energy, int band_index) const;
    double compute_dos_at_energy_and_band(double iso_energy, int band_index) const;
    double compute_dos_like_integral(double iso_energy, const std::vector<double>& f_values) const;
    inline double si_to_reduced_scale() const;


    std::vector<std::vector<double>> compute_dos_band_at_band(int         band_index,
                                                              double      min_energy,
                                                              double      max_energy,
                                                              int         num_threads,
                                                              std::size_t nb_points) const;
    std::vector<std::vector<double>> compute_dos_band_at_band_auto(int band_index, std::size_t nb_points, int num_threads) const;

    void export_k_points_to_file(const std::string& filename) const;
};

}  // namespace bz_mesh