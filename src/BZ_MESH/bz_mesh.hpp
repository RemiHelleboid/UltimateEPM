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

#include "iso_triangle.hpp"
#include "mesh_tetra.hpp"
#include "mesh_vertex.hpp"

namespace bz_mesh {

enum class BandType { valence, conduction };

class MeshBZ {
 private:
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

 public:
    MeshBZ() = default;

    void read_mesh_geometry_from_msh_file(const std::string& filename, double lattice_constant);
    void read_mesh_bands_from_msh_file(const std::string& filename);
    void add_new_band_energies_to_vertices(const std::vector<double>& energies_at_vertices);

    std::size_t get_number_vertices() const { return m_list_vertices.size(); }
    std::size_t get_number_elements() const { return m_list_tetrahedra.size(); }
    std::size_t get_number_bands() const { return m_min_band.size(); }

    const std::vector<Vertex>& get_list_vertices() const { return m_list_vertices; }
    const std::vector<Tetra>&  get_list_elements() const { return m_list_tetrahedra; }

    double compute_mesh_volume() const;
    double compute_iso_surface(double iso_energy, int band_index) const;
    double compute_dos_at_energy_and_band(double iso_energy, int band_index) const;

    std::vector<std::vector<double>> compute_dos_band_at_band(int         band_index,
                                                              double      min_energy,
                                                              double      max_energy,
                                                              std::size_t nb_points) const;
    std::vector<std::vector<double>> compute_dos_band_at_band_auto(int band_index, std::size_t nb_points) const;
};

}  // namespace bz_mesh