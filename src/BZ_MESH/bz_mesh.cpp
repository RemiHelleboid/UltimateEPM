/**
 * @file bz_mesh.cpp
 * @author Rémi Helleboid (remi.helleboid@st.com)
 * @brief
 * @version 0.1
 * @date 2022-07-14
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "bz_mesh.hpp"

#include <chrono>

#include "gmsh.h"
#include "omp.h"

#pragma omp declare reduction(merge : std::vector <double> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))


namespace bz_mesh {


/**
 * @brief Read the geometry of the mesh from the .msh file: the vertices and the elements are added to
 * the m_list_vertices and m_list_elements lists.
 *
 * @param filename
 */
void MeshBZ::read_mesh_geometry_from_msh_file(const std::string& filename) {
    std::cout << "Opening file " << filename << std::endl;
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 1);
    gmsh::open(filename);
    std::vector<std::size_t> nodeTags;
    std::vector<double>      nodeCoords;
    std::vector<double>      nodeParams;
    std::cout << "Reading vertices ..." << std::endl;
    gmsh::model::mesh::reclassifyNodes();
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, -1, -1, false, false);
    std::size_t size_nodes_tags        = nodeTags.size();
    std::size_t size_nodes_coordinates = nodeCoords.size();

    if (size_nodes_coordinates != 3 * size_nodes_tags) {
        throw std::runtime_error("Number of coordinates is not 3 times the number of vertices. Abort.");
    }

    m_list_vertices.reserve(size_nodes_tags);
    for (std::size_t index_vertex = 0; index_vertex < size_nodes_tags; ++index_vertex) {
        m_list_vertices.push_back(
            Vertex(index_vertex, nodeCoords[3 * index_vertex], nodeCoords[3 * index_vertex + 1], nodeCoords[3 * index_vertex + 2]));
    }
    std::cout << "Number of k-points vertices: " << m_list_vertices.size() << std::endl;

    std::cout << "Reading elements ..." << std::endl;
    // Get the mesh elements for the entity (dim, tag):
    const int                             dim = 3;
    const int                             tag = -1;
    std::vector<int>                      elemTypes;
    std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);
    if (elemTags.empty()) {
        std::cout << "ElementTags is zero when the mesh was imported... Abort.\n";
        throw std::runtime_error("ElementTags is zero when the mesh was imported... Abort.");
    }
    std::size_t number_elements = elemTags[0].size();

    std::cout << "Number of elements: " << number_elements << std::endl;
    std::cout << "Number of elements: " << elemNodeTags[0].size() << std::endl;

    if (elemNodeTags[0].size() != 4 * number_elements) {
        throw std::runtime_error("Number of elements vertices index is not 4 x NumberOfElements. Abort.");
    }

    m_list_tetrahedra.reserve(number_elements);
    for (std::size_t index_element = 0; index_element < number_elements; ++index_element) {
        const std::array<Vertex*, 4> array_element_vertices = {&m_list_vertices[elemNodeTags[0][4 * index_element] - 1],
                                                               &m_list_vertices[elemNodeTags[0][4 * index_element + 1] - 1],
                                                               &m_list_vertices[elemNodeTags[0][4 * index_element + 2] - 1],
                                                               &m_list_vertices[elemNodeTags[0][4 * index_element + 3] - 1]};
        Tetra                        new_tetra(index_element, array_element_vertices);
        m_list_tetrahedra.push_back(new_tetra);
    }
    std::cout << "Number of elements added: " << m_list_tetrahedra.size() << std::endl;

    gmsh::finalize();
}

/**
 * @brief Read the energy values for each band at every k-points (vertices) of the mesh.
 *
 * @param filename
 */
void MeshBZ::read_mesh_bands_from_msh_file(const std::string& filename) {
    std::cout << "Opening file " << filename << std::endl;
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 1);
    gmsh::open(filename);
    std::cout << "Read gmsh views (band energy values) ..." << std::endl;
    std::vector<int> viewTags;
    gmsh::view::getTags(viewTags);
    std::cout << "Number of view (bands) found: " << viewTags.size() << std::endl;
    int count_band = 0;
    for (auto&& tag : viewTags) {
        const int   index_view  = gmsh::view::getIndex(tag);
        std::string name_object = "View[" + std::to_string(index_view) + "].Name";
        std::string name_view;
        try {
            gmsh::option::getString(name_object, name_view);
        } catch (const std::exception& e) {
            std::cerr << e.what() << '\n';
        }

        std::string              type;
        std::vector<std::size_t> tags;
        double                   time;
        int                      numComp;
        std::vector<double>      data_view;
        gmsh::view::getHomogeneousModelData(tag, 0, type, tags, data_view, time, numComp);
        bool is_valence = data_view[0] <= 0.0;
        if (is_valence) {
            m_indices_valence_bands.push_back(count_band);
        } else {
            m_indices_conduction_bands.push_back(count_band);
        }
        count_band++;
        add_new_band_energies_to_vertices(data_view);
        auto minmax_band = std::minmax_element(data_view.begin(), data_view.end());
        m_min_band.push_back(*(minmax_band.first));
        m_max_band.push_back(*(minmax_band.second));
    }
    gmsh::finalize();
}

void MeshBZ::add_new_band_energies_to_vertices(const std::vector<double>& energies_at_vertices) {
    if (energies_at_vertices.size() != m_list_vertices.size()) {
        throw std::invalid_argument("The number of energy values does not match the number of vertices. Abort.");
    }
    for (std::size_t index_vtx = 0; index_vtx < m_list_vertices.size(); ++index_vtx) {
        m_list_vertices[index_vtx].add_band_energy_value(energies_at_vertices[index_vtx]);
    }
}

double MeshBZ::compute_mesh_volume() const {
    double total_volume = 0.0;
    for (auto&& tetra : m_list_tetrahedra) {
        total_volume += fabs(tetra.get_signed_volume());
    }
    return total_volume;
}

double MeshBZ::compute_iso_surface(double iso_energy, int band_index) const {
    auto   start     = std::chrono::high_resolution_clock::now();
    double total_dos = 0.0;

    for (auto&& tetra : m_list_tetrahedra) {
        total_dos += tetra.compute_tetra_dos_band(iso_energy, band_index);
    }
    auto end              = std::chrono::high_resolution_clock::now();
    auto total_time_count = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    return total_dos;
}

std::vector<std::vector<double>> MeshBZ::compute_dos_band_at_band(int         band_index,
                                                                  double      min_energy,
                                                                  double      max_energy,
                                                                  std::size_t nb_points) const {
    auto   start       = std::chrono::high_resolution_clock::now();
    double energy_step = (max_energy - min_energy) / (nb_points - 1);

    std::vector<double> list_energies{};
    std::vector<double> list_dos{};
#pragma omp parallel for schedule(dynamic) num_threads(32) reduction(merge : list_energies) reduction(merge : list_dos)
    for (std::size_t index_energy = 0; index_energy < nb_points; ++index_energy) {
        double energy = min_energy + index_energy * energy_step;
        double dos    = compute_iso_surface(energy, band_index);
        int    tid    = omp_get_thread_num();
#pragma omp critical
        list_energies.push_back(energy);
        list_dos.push_back(dos);
        //         std::cout << "\rComputing density of state at energy " << index_energy << "/" << nb_points << std::flush;
    }
    auto end              = std::chrono::high_resolution_clock::now();
    auto total_time_count = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nDOS for 1 band computed in  " << total_time_count / 1000.0 << "s" << std::endl;
    return {list_energies, list_dos};
}
std::vector<std::vector<double>> MeshBZ::compute_dos_band_at_band_auto(int band_index, std::size_t nb_points) const {
    const double margin_energy = 0.5;
    double min_energy = m_min_band[band_index] - margin_energy;
    double max_energy = m_max_band[band_index] + margin_energy;
    auto   start       = std::chrono::high_resolution_clock::now();
    double energy_step = (max_energy - min_energy) / (nb_points - 1);

    std::vector<double> list_energies{};
    std::vector<double> list_dos{};
#pragma omp parallel for schedule(dynamic) num_threads(32) reduction(merge : list_energies) reduction(merge : list_dos)
    for (std::size_t index_energy = 0; index_energy < nb_points; ++index_energy) {
        double energy = min_energy + index_energy * energy_step;
        double dos    = compute_iso_surface(energy, band_index);
        int    tid    = omp_get_thread_num();
#pragma omp critical
        list_energies.push_back(energy);
        list_dos.push_back(dos);
        //         std::cout << "\rComputing density of state at energy " << index_energy << "/" << nb_points << std::flush;
    }
    auto end              = std::chrono::high_resolution_clock::now();
    auto total_time_count = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nDOS for 1 band computed in  " << total_time_count / 1000.0 << "s" << std::endl;
    return {list_energies, list_dos};
}

}  // namespace bz_mesh