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

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <vector>

#include "gmsh.h"
#include "omp.h"
#include "rapidcsv.h"

#pragma omp declare reduction(merge : std::vector<double> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

namespace bz_mesh {

void MeshBZ::shift_bz_center(const vector3& center) {
    m_center = center;
    for (auto&& vtx : m_list_vertices) {
        vtx.shift_position(center);
    }
}

inline double MeshBZ::si_to_reduced_scale() const {
    // reduced k = (a / (2π)) * k_SI
    return m_material.get_lattice_constant_meter() / (2.0 * M_PI);
}

/**
 * @brief Read the geometry of the mesh from the .msh file: the vertices and the elements are added to
 * the m_list_vertices and m_list_elements lists.
 * All the points coordinates are re-normalized by the lattice constant passed as argument.
 *
 * @param filename
 * @param lattice_constant
 */
void MeshBZ::read_mesh_geometry_from_msh_file(const std::string& filename, bool normalize_by_fourier_factor) {
    std::cout << "Opening file " << filename << std::endl;
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 1000);
    gmsh::open(filename);
    std::vector<std::size_t> nodeTags;
    std::vector<double>      nodeCoords;
    std::vector<double>      nodeParams;
    std::cout << "Reading vertices ..." << std::endl;
    gmsh::model::mesh::reclassifyNodes();
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, -1, -1, false, false);
    std::size_t size_nodes_tags        = nodeTags.size();
    std::size_t size_nodes_coordinates = nodeCoords.size();
    m_node_tags                        = nodeTags;
    std::cout << "Number of nodes: " << size_nodes_tags << std::endl;

    if (size_nodes_coordinates != 3 * size_nodes_tags) {
        throw std::runtime_error("Number of coordinates is not 3 times the number of vertices. Abort.");
    }

    m_list_vertices.reserve(size_nodes_tags);
    double lattice_constant = m_material.get_lattice_constant_meter();
    std::cout << "Lattice const: " << lattice_constant << std::endl;
    std::cout << "V: " << std::pow(2.0 * M_PI, 3) / std::pow(lattice_constant, 3.0) << std::endl;
    const double fourier_factor = 2.0 * M_PI / lattice_constant;
    // const double fourier_factor       = 1;
    double normalization_factor = normalize_by_fourier_factor ? fourier_factor : 1.0;
    for (std::size_t index_vertex = 0; index_vertex < size_nodes_tags; ++index_vertex) {
        m_list_vertices.push_back(Vertex(index_vertex,
                                         normalization_factor * nodeCoords[3 * index_vertex],
                                         normalization_factor * nodeCoords[3 * index_vertex + 1],
                                         normalization_factor * nodeCoords[3 * index_vertex + 2]));
    }
    std::cout << "Number of k-points vertices: " << m_list_vertices.size() << std::endl;

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

    gmsh::finalize();
    m_total_volume = compute_mesh_volume();
    std::cout << "Total mesh volume: " << m_total_volume << std::endl;
    // Compute the reduced BZ volume
    double       Vcell      = std::pow(m_material.get_lattice_constant_meter(), 3) / 4.0;
    const double VBZ_theory = std::pow(2.0 * M_PI, 3) / Vcell;  // m^-3

    m_reduce_bz_factor = m_total_volume / VBZ_theory;
    std::cout << "Reduce BZ factor: " << m_reduce_bz_factor << std::endl;

    precompute_G_shifts();
    // vector3 b1 = {-1.0, 1.0, 1.0};
    // vector3 b2 = {1.0, -1.0, 1.0};
    // vector3 b3 = {1.0, 1.0, -1.0};
    Eigen::Vector3d  b1_SI                = {-1.0, 1.0, 1.0};
    Eigen::Vector3d  b2_SI                = {1.0, -1.0, 1.0};
    Eigen::Vector3d  b3_SI                = {1.0, 1.0, -1.0};
    constexpr double halfwidth_reduced    = 1.0;
    const double     ssi_to_reduced_scale = si_to_reduced_scale();
    init_reciprocal_basis(b1_SI, b2_SI, b3_SI, halfwidth_reduced, ssi_to_reduced_scale);
}

bbox_mesh MeshBZ::compute_bounding_box() const {
    double x_min = std::numeric_limits<double>::max();
    double y_min = std::numeric_limits<double>::max();
    double z_min = std::numeric_limits<double>::max();
    double x_max = std::numeric_limits<double>::min();
    double y_max = std::numeric_limits<double>::min();
    double z_max = std::numeric_limits<double>::min();
    for (auto&& vtx : m_list_vertices) {
        const vector3& position = vtx.get_position();
        x_min                   = std::min(x_min, position.x());
        y_min                   = std::min(y_min, position.y());
        z_min                   = std::min(z_min, position.z());
        x_max                   = std::max(x_max, position.x());
        y_max                   = std::max(y_max, position.y());
        z_max                   = std::max(z_max, position.z());
    }
    vector3 min_corner(x_min, y_min, z_min);
    vector3 max_corner(x_max, y_max, z_max);
    return bbox_mesh(min_corner, max_corner);
}

void MeshBZ::build_search_tree() {
    bbox_mesh mesh_bbox = compute_bounding_box();
    std::cout << "Mesh bounding box: " << mesh_bbox << std::endl;
    mesh_bbox.dilate(1.05);
    // bbox_mesh.translate({0.0, 0.0, 0.0});
    auto start    = std::chrono::high_resolution_clock::now();
    m_search_tree = std::make_unique<Octree_mesh>(get_list_p_tetra(), mesh_bbox);
    auto end      = std::chrono::high_resolution_clock::now();
    auto total    = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Octree built in " << total / 1000.0 << "s" << std::endl;
}

Tetra* MeshBZ::find_tetra_at_location(const vector3& location) const { return m_search_tree->find_tetra_at_location(location); }

/**
 * @brief Get the nearest k index object.
 * Brute force search of the nearest k-point index. :(
 *
 * @param k
 * @return std::size_t
 */
std::size_t MeshBZ::get_nearest_k_index(const Vector3D<double>& k) const {
    vector3     K(k.X, k.Y, k.Z);
    std::size_t index_nearest_k = 0;
    double      min_distance    = std::numeric_limits<double>::max();
    for (std::size_t index_k = 0; index_k < m_list_vertices.size(); ++index_k) {
        double distance = (K - m_list_vertices[index_k].get_position()).norm();
        if (distance < min_distance) {
            min_distance    = distance;
            index_nearest_k = index_k;
        }
    }
    return index_nearest_k;
}

std::size_t MeshBZ::get_nearest_k_index(const vector3& k) const {
    std::size_t index_nearest_k = 0;
    double      min_distance    = std::numeric_limits<double>::max();
    for (std::size_t index_k = 0; index_k < m_list_vertices.size(); ++index_k) {
        double distance = (k - m_list_vertices[index_k].get_position()).norm();
        if (distance < min_distance) {
            min_distance    = distance;
            index_nearest_k = index_k;
        }
    }
    return index_nearest_k;
}

/**
 * @brief Read the energy values for each band at every k-points (vertices) of the mesh.
 *
 * @param filename
 */
void MeshBZ::read_mesh_bands_from_msh_file(const std::string& filename, int nb_bands_to_load) {
    std::cout << "Opening file " << filename << std::endl;
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 0);
    gmsh::open(filename);
    // std::cout << "Read gmsh views (band energy values) ..." << std::endl;
    std::vector<int> viewTags;
    gmsh::view::getTags(viewTags);
    // std::cout << "Number of view (bands) found: " << viewTags.size() << std::endl;
    int count_band = 0;

    for (auto&& tag : viewTags) {
        if (nb_bands_to_load != -1 && count_band >= nb_bands_to_load) break;

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
        bool is_valence = data_view[0] <= 0.1;
        if (is_valence) {
            m_indices_valence_bands.push_back(count_band);
        } else {
            m_indices_conduction_bands.push_back(count_band);
        }
        count_band++;
        auto minmax_band = std::minmax_element(data_view.begin(), data_view.end());
        m_min_band.push_back(*(minmax_band.first));
        m_max_band.push_back(*(minmax_band.second));
        add_new_band_energies_to_vertices(data_view);
    }
    gmsh::finalize();
    // PRINT INDICES
    std::cout << "Number of bands loaded: " << m_min_band.size() << std::endl;
    std::cout << "Number of valence bands: " << m_indices_valence_bands.size() << std::endl;
    std::cout << "Number of conduction bands: " << m_indices_conduction_bands.size() << std::endl;
    for (auto&& idx : m_indices_valence_bands) {
        std::cout << "Valence band index: " << idx << std::endl;
    }
    for (auto&& idx : m_indices_conduction_bands) {
        std::cout << "Conduction band index: " << idx << std::endl;
    }

    auto_shift_conduction_band_energies();
    auto_set_positive_valence_band_energies();
    compute_min_max_energies_at_tetras();
    compute_energy_gradient_at_tetras();
    for (auto&& tetra : m_list_tetrahedra) {
        tetra.pre_compute_sorted_slots_per_band();
    }

    // Print band info
    std::cout << "Number of bands loaded: " << m_min_band.size() << std::endl;
    std::cout << "Number of valence bands: " << m_indices_valence_bands.size() << std::endl;
    std::cout << "Number of conduction bands: " << m_indices_conduction_bands.size() << std::endl;
    for (std::size_t i = 0; i < m_min_band.size(); ++i) {
        std::cout << "Band " << i << ": min = " << m_min_band[i] << " eV, max = " << m_max_band[i] << " eV";
        if (std::find(m_indices_valence_bands.begin(), m_indices_valence_bands.end(), i) != m_indices_valence_bands.end()) {
            std::cout << " (valence band)";
        } else if (std::find(m_indices_conduction_bands.begin(), m_indices_conduction_bands.end(), i) != m_indices_conduction_bands.end()) {
            std::cout << " (conduction band)";
        }
        std::cout << std::endl;
    }
}

void MeshBZ::precompute_dos_tetra(double energy_step) {
    std::cout << "Precomputing DOS per tetrahedra with energy step = " << energy_step << " eV ..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
#pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < m_list_tetrahedra.size(); ++i) {
        m_list_tetrahedra[i].precompute_dos_on_energy_grid_per_band(energy_step);
    }
    
}


void MeshBZ::recompute_min_max_energies() {
    m_min_band.clear();
    m_max_band.clear();
    int nb_bands = m_list_vertices[0].get_number_bands();
    m_min_band.resize(nb_bands, std::numeric_limits<double>::max());
    m_max_band.resize(nb_bands, std::numeric_limits<double>::lowest());
    for (auto&& vtx : m_list_vertices) {
        const auto& energies = vtx.get_band_energies();
        for (int i = 0; i < nb_bands; ++i) {
            if (energies[i] < m_min_band[i]) m_min_band[i] = energies[i];
            if (energies[i] > m_max_band[i]) m_max_band[i] = energies[i];
        }
    }
    // Print band info
    std::cout << "Number of bands: " << m_min_band.size() << std::endl;
    for (std::size_t i = 0; i < m_min_band.size(); ++i) {
        std::cout << "Band " << i << ": min = " << m_min_band[i] << " eV, max = " << m_max_band[i] << " eV";
        if (std::find(m_indices_valence_bands.begin(), m_indices_valence_bands.end(), i) != m_indices_valence_bands.end()) {
            std::cout << " (valence band)";
        } else if (std::find(m_indices_conduction_bands.begin(), m_indices_conduction_bands.end(), i) != m_indices_conduction_bands.end()) {
            std::cout << " (conduction band)";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

/**
 * @brief Keep only a subset of bands (valence and conduction) in the mesh. Remove higher (relative) bands.
 * @param nb_valence_bands Number of valence bands to keep
 * @param nb_conduction_bands Number of conduction bands to keep
 */
void MeshBZ::keep_only_bands(const int nb_valence_bands, const int nb_conduction_bands) {
    int nb_valence_to_remove    = static_cast<int>(m_indices_valence_bands.size()) - nb_valence_bands;
    int nb_conduction_to_remove = static_cast<int>(m_indices_conduction_bands.size()) - nb_conduction_bands;
    if (nb_valence_to_remove < 0 || nb_conduction_to_remove < 0) {
        throw std::runtime_error("Cannot keep more bands than available.");
    }
    std::cout << "Removing " << nb_valence_to_remove << " valence bands and " << nb_conduction_to_remove << " conduction bands."
              << std::endl;
    // Remove higher valence bands
    for (int i = 0; i < nb_valence_to_remove; ++i) {
        int band_to_remove = m_indices_valence_bands.back();
        m_indices_valence_bands.pop_back();
        for (auto&& vtx : m_list_vertices) {
            vtx.remove_band_energy(band_to_remove);
        }
    }
    // Remove higher conduction bands
    for (int i = 0; i < nb_conduction_to_remove; ++i) {
        int band_to_remove = m_indices_conduction_bands.back();
        m_indices_conduction_bands.pop_back();
        for (auto&& vtx : m_list_vertices) {
            vtx.remove_band_energy(band_to_remove);
        }
    }
    // Recompute min/max band energies
    recompute_min_max_energies();
    std::cout << "After removing bands:" << std::endl;
    std::cout << "Number of valence bands: " << m_indices_valence_bands.size() << std::endl;
    std::cout << "Number of conduction bands: " << m_indices_conduction_bands.size() << std::endl;
}

void MeshBZ::add_new_band_energies_to_vertices(const std::vector<double>& energies_at_vertices) {
    if (energies_at_vertices.size() != m_list_vertices.size()) {
        throw std::invalid_argument("The number of energy values does not match the number of vertices. Abort.");
    }
    for (std::size_t index_vtx = 0; index_vtx < m_list_vertices.size(); ++index_vtx) {
        m_list_vertices[index_vtx].add_band_energy_value(energies_at_vertices[index_vtx]);
    }
}

void MeshBZ::compute_min_max_energies_at_tetras() {
    for (auto&& tetra : m_list_tetrahedra) {
        tetra.compute_min_max_energies_at_bands();
    }
}

void MeshBZ::auto_set_positive_valence_band_energies() {
    for (int vband : m_indices_valence_bands) {
        for (auto&& vtx : m_list_vertices) {
            double energy = vtx.get_energy_at_band(vband);
            vtx.set_band_energy(vband, std::fabs(energy));
        }
    }
}

void MeshBZ::auto_shift_conduction_band_energies() {
    if (m_indices_conduction_bands.empty()) return;

    // Find the highest valence band maximum
    double max_valence = std::numeric_limits<double>::lowest();
    for (int vband : m_indices_valence_bands) {
        max_valence = std::max(max_valence, m_max_band[vband]);
    }
    if (m_indices_valence_bands.empty()) {
        max_valence = 0.0;
    }
    std::cout << "Max valence band energy: " << max_valence << " eV\n";

    // Find the lowest conduction band minimum
    double min_conduction = std::numeric_limits<double>::max();
    for (int cband : m_indices_conduction_bands) {
        min_conduction = std::min(min_conduction, m_min_band[cband]);
    }
    std::cout << "Min conduction band energy: " << min_conduction << " eV\n";

    double band_gap = min_conduction - max_valence;
    std::cout << "Computed band gap: " << band_gap << " eV\n";

    if (band_gap < 0.1) {
        std::cout << "Warning: computed band gap is very small or negative (" << band_gap << " eV). No shift applied.\n";
        return;
    }

    // Shift all conduction bands so that the conduction band minimum aligns with the valence band maximum
    double shift_amount = -min_conduction + max_valence;
    std::cout << "Shifting conduction bands by " << shift_amount << " eV to align CBM with VBM.\n";

    // void set_band_energy(std::size_t index_band, double new_energy) {

    for (auto&& vtx : m_list_vertices) {
        for (int cband : m_indices_conduction_bands) {
            double old_energy = vtx.get_energy_at_band(cband);
            double new_energy = old_energy + shift_amount;
            vtx.set_band_energy(cband, new_energy);
        }
    }

    // Recompute per-band min/max
    m_min_band.clear();
    m_max_band.clear();
    for (std::size_t b = 0; b < m_list_vertices[0].get_number_bands(); ++b) {
        double bmin = std::numeric_limits<double>::max();
        double bmax = std::numeric_limits<double>::lowest();
        for (auto&& vtx : m_list_vertices) {
            double e = vtx.get_energy_at_band(b);
            bmin     = std::min(bmin, e);
            bmax     = std::max(bmax, e);
        }
        m_min_band.push_back(bmin);
        m_max_band.push_back(bmax);
    }

    std::cout << "Post-shift band extrema:\n";
    for (std::size_t b = 0; b < m_min_band.size(); ++b) {
        std::cout << " Band " << b << ": min = " << m_min_band[b] << " eV, max = " << m_max_band[b] << " eV\n";
    }
    std::cout << std::endl;
}

void MeshBZ::set_bands_in_right_order() {
    if (m_indices_valence_bands.empty() || m_indices_conduction_bands.empty()) {
        std::cout << "Warning: valence or conduction band indices are empty. Cannot reorder bands.\n";
        return;
    }

    // DEBUG PROVISOIRE
    for (auto&& vtx : m_list_vertices) {
        vtx.swap_bands(m_indices_valence_bands[0], m_indices_valence_bands.back());
        vtx.swap_bands(m_indices_valence_bands[1], m_indices_valence_bands[m_indices_valence_bands.size() - 2]);
        std::cout << "Swapped bands for vertex " << vtx.get_index() << ": band " << m_indices_valence_bands[0] << " with band "
                  << m_indices_valence_bands.back() << ", band " << m_indices_valence_bands[1] << " with band "
                  << m_indices_valence_bands[m_indices_valence_bands.size() - 2] << std::endl;
    }
    recompute_min_max_energies();
}

void MeshBZ::compute_energy_gradient_at_tetras() {
    for (auto&& tetra : m_list_tetrahedra) {
        tetra.compute_gradient_energy_at_bands();
    }
}

double MeshBZ::compute_mesh_volume() const {
    double total_volume = 0.0;
    for (auto&& tetra : m_list_tetrahedra) {
        total_volume += std::fabs(tetra.get_signed_volume());
        // std::cout << "Tetra " << tetra.get_index() << " volume: " << tetra.get_signed_volume() << std::endl;
    }
    // total_volume *= (1.0 / pow(2.0 * M_PI, 3.0));
    return total_volume;
}

double MeshBZ::compute_iso_surface(double iso_energy, int band_index) const {
    double total_dos = 0.0;
    for (auto&& tetra : m_list_tetrahedra) {
        total_dos += tetra.compute_tetra_iso_surface_energy_band(iso_energy, band_index);
    }

    return total_dos;
}

double MeshBZ::compute_dos_at_energy_and_band(double iso_energy, int band_index, bool use_interp) const {
    double total_dos = 0.0;
    for (auto&& tetra : m_list_tetrahedra) {
        if (use_interp) {
            total_dos += tetra.interpolate_dos_at_energy_per_band(iso_energy, band_index);
        } else {
            total_dos += tetra.compute_tetra_dos_energy_band(iso_energy, band_index);
        }
    }
    total_dos *= get_reduce_bz_factor();
    total_dos *= m_spin_degeneracy;
    return total_dos;
}

std::vector<std::vector<double>> MeshBZ::compute_dos_band_at_band(int         band_index,
                                                                  double      min_energy,
                                                                  double      max_energy,
                                                                  int         num_threads,
                                                                  std::size_t nb_points,
                                                                  bool        use_interp) const {
    auto   start       = std::chrono::high_resolution_clock::now();
    double energy_step = (max_energy - min_energy) / (nb_points - 1);

    std::vector<double> list_energies{};
    std::vector<double> list_dos{};
#pragma omp parallel for schedule(dynamic) num_threads(num_threads) reduction(merge : list_energies) reduction(merge : list_dos)
    for (std::size_t index_energy = 0; index_energy < nb_points; ++index_energy) {
        double energy = min_energy + index_energy * energy_step;
        double dos    = compute_dos_at_energy_and_band(energy, band_index);
        list_energies.push_back(energy);
        list_dos.push_back(dos);
        //         std::cout << "\rComputing density of state at energy " << index_energy << "/" << nb_points << std::flush;
    }
    auto end              = std::chrono::high_resolution_clock::now();
    auto total_time_count = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nDOS for 1 band computed in  " << total_time_count / 1000.0 << "s" << std::endl;
    return {list_energies, list_dos};
}

std::vector<std::vector<double>> MeshBZ::compute_dos_band_at_band_auto(int         band_index,
                                                                       std::size_t nb_points,
                                                                       int         num_threads,
                                                                       bool        use_interp) const {
    if (band_index < 0 || band_index >= static_cast<int>(m_min_band.size())) {
        throw std::out_of_range("Band index out of range in compute_dos_band_at_band_auto.");
    }

    const double margin_energy = 0.1;
    double       min_energy    = m_min_band[band_index] - margin_energy;
    double       max_energy    = m_max_band[band_index] + margin_energy;
    auto         start         = std::chrono::high_resolution_clock::now();
    double       energy_step   = (max_energy - min_energy) / (nb_points - 1);

    std::vector<double> list_energies{};
    std::vector<double> list_dos{};
#pragma omp parallel for schedule(dynamic) num_threads(num_threads) reduction(merge : list_energies) reduction(merge : list_dos)
    for (std::size_t index_energy = 0; index_energy < nb_points; ++index_energy) {
        double energy = min_energy + index_energy * energy_step;
        double dos    = compute_dos_at_energy_and_band(energy, band_index, use_interp);
        list_energies.push_back(energy);
        list_dos.push_back(dos);
    }
    auto end              = std::chrono::high_resolution_clock::now();
    auto total_time_count = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nDOS for 1 band computed in  " << total_time_count / 1000.0 << "s" << std::endl;
    return {list_energies, list_dos};
}

/**
 * @brief Draw a random tetrahedron index in the mesh, on a iso-energy surface.
 *
 * @param energy
 * @param idx_band
 * @param random_generator
 * @return std::size_t
 */
std::size_t MeshBZ::draw_random_tetrahedron_index_with_dos_probability(double        energy,
                                                                       std::size_t   idx_band,
                                                                       std::mt19937& random_generator) const {
    std::vector<double> list_dos;
    list_dos.reserve(m_list_tetrahedra.size());
    for (auto&& tetra : m_list_tetrahedra) {
        double dos = tetra.compute_tetra_dos_energy_band(energy, idx_band);
        list_dos.push_back(dos);
    }
    std::discrete_distribution<std::size_t> distribution(list_dos.begin(), list_dos.end());
    return distribution(random_generator);
}

/**
 * @brief Draw a random k-vector in the mesh, on a iso-energy surface.
 * The k-vector is drawn with a probability locally proportional to the DOS.
 *
 * @param energy
 * @param idx_band
 * @param random_generator
 * @return vector3
 */
vector3 MeshBZ::draw_random_k_point_at_energy(double energy, std::size_t idx_band, std::mt19937& random_generator) const {
    if (energy < m_min_band[idx_band] || energy > m_max_band[idx_band]) {
        throw std::runtime_error("Energy is out of range");
    }
    const std::size_t index_tetra = draw_random_tetrahedron_index_with_dos_probability(energy, idx_band, random_generator);
    return m_list_tetrahedra[index_tetra].draw_random_uniform_point_at_energy(energy, idx_band, random_generator);
}

void MeshBZ::export_k_points_to_file(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::invalid_argument("Could not open file " + filename + " for writing.");
    }
    for (auto&& k_point : m_list_vertices) {
        file << k_point.get_position().x() << "," << k_point.get_position().y() << "," << k_point.get_position().z() << std::endl;
    }
    file.close();
}


// helper: half-width in reduced units.
static inline double bz_halfwidth_reduced() { return 1.0; }  // <- your case

bool MeshBZ::is_inside_mesh_geometry(const vector3& k) const {
    const double     s   = si_to_reduced_scale();   // converts SI (1/m) -> reduced (unitless)
    const double     hw  = bz_halfwidth_reduced();  // 1.0 for [-1,1], 0.5 for [-0.5,0.5]
    constexpr double eps = 1e-12;

    const double kx = k.x() * s;
    const double ky = k.y() * s;
    const double kz = k.z() * s;

    const double ax = std::fabs(kx);
    const double ay = std::fabs(ky);
    const double az = std::fabs(kz);

    const bool cond1 = (ax <= hw + eps) && (ay <= hw + eps) && (az <= hw + eps);
    const bool cond2 = (ax + ay + az) <= (1.5 * hw + eps);

    return cond1 && cond2;
}

void MeshBZ::precompute_G_shifts() {
    vector3      b1 = {-1.0, 1.0, 1.0};
    vector3      b2 = {1.0, -1.0, 1.0};
    vector3      b3 = {1.0, 1.0, -1.0};
    const double s  = si_to_reduced_scale();
    b1 /= s;
    b2 /= s;
    b3 /= s;

    const int maxShell = 5;  // adjust as needed

    // shell 0
    m_Gshifts.push_back({0, 0, 0});

    // by increasing |n1|+|n2|+|n3|
    for (int L1 = 1; L1 <= maxShell; ++L1) {
        for (int n1 = -L1; n1 <= L1; ++n1) {
            for (int n2 = -L1; n2 <= L1; ++n2) {
                for (int n3 = -L1; n3 <= L1; ++n3) {
                    if (n1 == 0 && n2 == 0 && n3 == 0) continue;
                    if (std::abs(n1) + std::abs(n2) + std::abs(n3) == L1) {
                        m_Gshifts.push_back(n1 * b1 + n2 * b2 + n3 * b3);
                    }
                }
            }
        }
    }
}

// vector3 MeshBZ::fold_ws_bcc_brut(const vector3d& k_SI) const noexcept {
//     if (inside_ws_bcc(k_SI)) return k_SI;
//     for (const auto& G : m_Gshifts) {
//         const vector3d k_try = k_SI - G;
//         if (inside_ws_bcc(k_try)) return k_try;
//     }
//     throw std::runtime_error("Could not fold k into WS-BZ");
// }

vector3 MeshBZ::retrieve_k_inside_mesh_geometry(const vector3& k) const {
    for (const auto& G : m_Gshifts) {
        const vector3 kG = k + G;
        if (is_inside_mesh_geometry(kG)) {
            return kG;
        }
    }
    std::cout << "No k-point inside the mesh geometry found for k: " << k << std::endl;
    throw std::runtime_error("No k-point inside the mesh geometry found.");
}

// ---------------------------------------------------------
// O(1) WS-BZ folding: init (call once after you know the lattice)
// ---------------------------------------------------------
void MeshBZ::init_reciprocal_basis(const Eigen::Vector3d& b1_SI,
                                   const Eigen::Vector3d& b2_SI,
                                   const Eigen::Vector3d& b3_SI,
                                   double                 halfwidth_reduced,
                                   double                 si_to_reduced) {
    m_recip_B.col(0) = b1_SI;
    m_recip_B.col(1) = b2_SI;
    m_recip_B.col(2) = b3_SI;
    m_recip_Bi       = m_recip_B.inverse();

    m_bz_halfwidth = halfwidth_reduced;  // e.g. 0.5 for [-0.5,0.5]
    m_si2red       = si_to_reduced;      // must match your is_inside_mesh_geometry() scale
}

// ---------------------------------------------------------
// O(1) WS-BZ folding: pure function
// ---------------------------------------------------------
vector3 MeshBZ::fold_ws_bcc(const vector3& k_SI) const noexcept {
    // Convert to Eigen for the tiny linear algebra steps
    Eigen::Vector3d ke(k_SI.x(), k_SI.y(), k_SI.z());

    // A) Nearest-lattice wrap (Babai rounding) into the primitive parallelepiped
    Eigen::Vector3d r  = m_recip_Bi * ke;    // reduced coords in basis of reciprocal vectors
    Eigen::Array3d  n  = r.array().round();  // nearest reciprocal-lattice node
    Eigen::Vector3d k0 = ke - m_recip_B * n.matrix();

    // B) Apply WS (truncated-octahedron) plane tests in your reduced frame
    Eigen::Vector3d kr = m_si2red * k0;

    const double hw      = m_bz_halfwidth;
    const double sum_lim = 1.5 * hw;
    const double eps     = 1e-12 * std::max(1.0, hw);

    // A few corrections always suffice after Babai; 3 passes is plenty.
    for (int it = 0; it < 3; ++it) {
        const double ax = std::abs(kr.x()), ay = std::abs(kr.y()), az = std::abs(kr.z());
        bool         moved = false;

        // 6 square faces: |x|,|y|,|z| ≤ hw
        if (ax > hw + eps) {
            kr.x() -= (kr.x() > 0 ? 1.0 : -1.0) * 2.0 * hw;
            moved = true;
        }
        if (ay > hw + eps) {
            kr.y() -= (kr.y() > 0 ? 1.0 : -1.0) * 2.0 * hw;
            moved = true;
        }
        if (az > hw + eps) {
            kr.z() -= (kr.z() > 0 ? 1.0 : -1.0) * 2.0 * hw;
            moved = true;
        }

        // 8 hex faces: |x| + |y| + |z| ≤ 1.5*hw
        if ((ax + ay + az) > (sum_lim + eps)) {
            const int    idx    = (ax >= ay && ax >= az) ? 0 : (ay >= az ? 1 : 2);
            const double excess = (ax + ay + az) - sum_lim;
            if (idx == 0)
                kr.x() -= (kr.x() > 0 ? 1.0 : -1.0) * excess;
            else if (idx == 1)
                kr.y() -= (kr.y() > 0 ? 1.0 : -1.0) * excess;
            else
                kr.z() -= (kr.z() > 0 ? 1.0 : -1.0) * excess;
            moved = true;
        }

        if (!moved) break;  // inside WS
    }

    // Back to SI and convert to your vector3
    Eigen::Vector3d kf = kr / m_si2red;
    return vector3{kf.x(), kf.y(), kf.z()};
}

// ---------------------------------------------------------
// Predicate (branch-lean), consistent with your existing checker
// ---------------------------------------------------------
bool MeshBZ::inside_ws_bcc(const vector3& k_SI) const noexcept {
    Eigen::Vector3d kr = m_si2red * Eigen::Vector3d(k_SI.x(), k_SI.y(), k_SI.z());
    const double    hw = m_bz_halfwidth, eps = 1e-12 * std::max(1.0, hw);
    const double    ax = std::abs(kr.x()), ay = std::abs(kr.y()), az = std::abs(kr.z());
    // use bitwise & so all comparisons are evaluated (branchless)
    return (ax <= hw + eps) & (ay <= hw + eps) & (az <= hw + eps) & ((ax + ay + az) <= (1.5 * hw + eps));
}

// bool is_in_first_BZ(const Vector3D<double>& k, bool one_eighth = false) {
//     bool cond_1      = fabs(k.X) <= 1.0 && fabs(k.Y) <= 1.0 && fabs(k.Z) <= 1.0;
//     bool cond_2      = fabs(k.X) + fabs(k.Y) + fabs(k.Z) <= 3.0 / 2.0;
//     bool cond_eighth = (k.X >= 0.0 && k.Y >= 0.0 && k.Z >= 0.0);
//     return cond_1 && cond_2 && (one_eighth ? cond_eighth : true);
// }

// def IsInIrreducibleWedge(k):
//     return (k[2] >= 0.0 and k[2] <= k[1] and k[1] <= k[0] and k[0] <= 1.0) and \
//         (np.sum(k) <= 3.0/2.0)


bool MeshBZ::is_irreducible_wedge(const vector3& k_SI) const noexcept {
    const double x = k_SI.x() * si_to_reduced_scale();
    const double y = k_SI.y() * si_to_reduced_scale();
    const double z = k_SI.z() * si_to_reduced_scale();
    // std::cout << "Checking irreducible wedge for k = (" << x << ", " << y << ", " << z << ")" << std::endl;
    constexpr double eps = 1e-12;
    return (x >= -eps) && (y >= -eps) && (z >= -eps) && (x <= y + eps) && (y <= z + eps) && (z <= 1.0 + eps) && ((x + y + z) <= 1.5 + eps);
}

/**
 * @brief Read phonon scattering rates from a file.
 *
 * @param path
 */
void MeshBZ::read_phonon_scattering_rates_from_file(const std::filesystem::path& path) {
    std::cout << "Reading phonon scattering rates (CSV) from file " << path.string() << " ..." << std::endl;

    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open file " + path.string());
    }

    m_list_phonon_scattering_rates.clear();
    m_list_phonon_scattering_rates.resize(m_list_vertices.size());

    std::string line;
    std::size_t line_no = 0;

    for (std::size_t idx_vtx = 0; idx_vtx < m_list_vertices.size(); ++idx_vtx) {
        const auto&       vertex   = m_list_vertices[idx_vtx];
        const std::size_t nbands   = vertex.get_number_bands();
        auto&             per_band = m_list_phonon_scattering_rates[idx_vtx];
        per_band.resize(nbands);

        for (std::size_t idx_band = 0; idx_band < nbands; ++idx_band) {
            do {
                if (!std::getline(in, line)) {
                    std::cout << "Line no: " << line_no << " " << line << std::endl;
                    throw std::runtime_error("Unexpected EOF at vertex " + std::to_string(idx_vtx) + ", band " + std::to_string(idx_band));
                }
                ++line_no;
            } while (line.empty() || line[0] == '#' || line[0] == ';');

            for (char& c : line)
                if (c == ',') c = ' ';
            std::istringstream iss(line);

            std::size_t band_idx_file{};
            double      energy_file{};
            Rate8       rates{};

            if (!(iss >> band_idx_file >> energy_file >> rates[static_cast<std::size_t>(PhMode::ALO)] >>
                  rates[static_cast<std::size_t>(PhMode::ALA)] >> rates[static_cast<std::size_t>(PhMode::ATO)] >>
                  rates[static_cast<std::size_t>(PhMode::ATA)] >> rates[static_cast<std::size_t>(PhMode::ELO)] >>
                  rates[static_cast<std::size_t>(PhMode::ELA)] >> rates[static_cast<std::size_t>(PhMode::ETO)] >>
                  rates[static_cast<std::size_t>(PhMode::ETA)])) {
                throw std::runtime_error("Malformed CSV line " + std::to_string(line_no));
            }
            if (band_idx_file == nbands) {
                // More bands in the file than in the mesh: ignore extra bands
                break;
            }

            per_band[idx_band] = rates;
        }
    }
    in.close();

    std::cout << "Finished reading phonon scattering rates for " << m_list_vertices.size() << " vertices." << std::endl;
}

/**
 * @brief Interpolate the phonon scattering rate at a given location for a given band.
 *
 * @param location
 * @param idx_band
 * @return Rate8
 */
Rate8 MeshBZ::interpolate_phonon_scattering_rate_at_location(const vector3& location, const std::size_t& idx_band) const {
    // Find the tetrahedron containing the location
    const Tetra* tetra = find_tetra_at_location(location);
    if (!tetra) {
        throw std::runtime_error("Location is not inside any tetrahedron");
    }

    // Get the vertex indices of the tetrahedron
    const auto& vertex_indices = tetra->get_index_vertices_with_sorted_energy_at_band(idx_band);

    // Interpolate the scattering rates at the vertices
    Rate8 rates;
    for (std::size_t i = 0; i < vertex_indices.size(); ++i) {
        const auto& vertex = m_list_vertices[vertex_indices[i]];
        for (std::size_t idx_mode = 0; idx_mode < rates.size(); ++idx_mode) {
            rates[idx_mode] += m_list_phonon_scattering_rates[vertex.get_index()][idx_band][idx_mode];
        }
    }

    // Average the rates
    for (auto& rate : rates) {
        rate /= static_cast<double>(vertex_indices.size());
    }
    return rates;
}

// Helper: sum 8 phonon-mode rates (ALO, ALA, ATO, ATA, ELO, ELA, ETO, ETA)
inline double MeshBZ::sum_modes(const Rate8& r) const noexcept {
    double s = 0.0;
    for (std::size_t m = 0; m < kModeCount; ++m)
        s += r[m];
    return s;
}

double MeshBZ::compute_P_Gamma() const {
    double pgamma_max = 0.0;
    for (const auto& perVertex : m_list_phonon_scattering_rates) {
        for (const auto& rate8 : perVertex) {
            const double tot = sum_modes(rate8);
            if (tot > pgamma_max) {
                pgamma_max = tot;
            }
        }
    }
    return pgamma_max;
}

}  // namespace bz_mesh