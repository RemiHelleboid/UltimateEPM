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

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "gmsh.h"
#include "integrals.hpp"
#include "octree_bz.hpp"
#include "omp.h"
#include "physical_constants.hpp"
#include "rapidcsv.h"
#include "vector.hpp"

#pragma omp declare reduction(merge : std::vector<double> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

namespace uepm::mesh_bz {

void MeshBZ::shift_bz_center(const vector3& center) {
    m_center = center;
    for (auto&& vtx : m_list_vertices) {
        vtx.shift_position(center);
    }
}

double MeshBZ::si_to_reduced_scale() const noexcept {
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
    m_vertex_to_tetrahedra.resize(m_list_vertices.size());
    for (std::size_t index_element = 0; index_element < number_elements; ++index_element) {
        const std::array<Vertex*, 4> array_element_vertices = {&m_list_vertices[elemNodeTags[0][4 * index_element] - 1],
                                                               &m_list_vertices[elemNodeTags[0][4 * index_element + 1] - 1],
                                                               &m_list_vertices[elemNodeTags[0][4 * index_element + 2] - 1],
                                                               &m_list_vertices[elemNodeTags[0][4 * index_element + 3] - 1]};
        Tetra                        new_tetra(index_element, array_element_vertices);
        m_list_tetrahedra.push_back(new_tetra);
        for (std::size_t i = 0; i < 4; ++i) {
            m_vertex_to_tetrahedra[elemNodeTags[0][4 * index_element + i] - 1].push_back(index_element);
        }
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
    Eigen::Vector3d  b1_SI                = {-1.0, 1.0, 1.0};
    Eigen::Vector3d  b2_SI                = {1.0, -1.0, 1.0};
    Eigen::Vector3d  b3_SI                = {1.0, 1.0, -1.0};
    constexpr double halfwidth_reduced    = 1.0;
    const double     ssi_to_reduced_scale = si_to_reduced_scale();
    init_reciprocal_basis(b1_SI, b2_SI, b3_SI, halfwidth_reduced, ssi_to_reduced_scale);

    // build_search_tree();
}

bbox_mesh MeshBZ::compute_bounding_box() const {
    double x_min = std::numeric_limits<double>::max();
    double y_min = std::numeric_limits<double>::max();
    double z_min = std::numeric_limits<double>::max();
    double x_max = std::numeric_limits<double>::lowest();
    double y_max = std::numeric_limits<double>::lowest();
    double z_max = std::numeric_limits<double>::lowest();
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
    mesh_bbox.dilate(1.10);
    // Translate the bounding box by a quite random small vector to avoid degenerate cases...
    const double  translation_unit_factor = 1.876473876e-6;
    const vector3 translation             = translation_unit_factor * (mesh_bbox.max_corner() - mesh_bbox.min_corner());
    std::cout << "Translate bbox by " << translation << std::endl;
    mesh_bbox.translate(translation);
    auto start    = std::chrono::high_resolution_clock::now();
    m_search_tree = std::make_unique<Octree_mesh>(get_list_p_tetra(), mesh_bbox);
    auto end      = std::chrono::high_resolution_clock::now();
    auto total    = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Octree built in " << total / 1000.0 << "s" << std::endl;
}

Tetra* MeshBZ::find_tetra_at_location(const vector3& location) const { return m_search_tree->find_tetra_at_location(location); }

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
void MeshBZ::read_mesh_bands_from_msh_file(const std::string& filename, int nb_conduction_bands, int nb_valence_bands) {
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
            BandInfo band_info = {MeshParticleType::valence, m_valence_bands.count};
            m_band_info.push_back(band_info);
            m_valence_bands.global_start_index = (m_valence_bands.count == 0) ? count_band : m_valence_bands.global_start_index;
            m_valence_bands.count++;
        } else {
            BandInfo band_info = {MeshParticleType::conduction, m_conduction_bands.count};
            m_band_info.push_back(band_info);
            m_conduction_bands.global_start_index = (m_conduction_bands.count == 0) ? count_band : m_conduction_bands.global_start_index;
            m_conduction_bands.count++;
        }
        count_band++;
        auto minmax_band = std::minmax_element(data_view.begin(), data_view.end());
        m_min_band.push_back(*(minmax_band.first));
        m_max_band.push_back(*(minmax_band.second));
        add_new_band_energies_to_vertices(data_view);
    }
    m_nb_bands_total = count_band;
    gmsh::finalize();
    // PRINT INFO
    std::cout << "Number of bands loaded: " << get_number_bands_total() << std::endl;
    print_band_info();
    set_bands_in_right_order();
    if (nb_valence_bands > 0 || nb_conduction_bands > 0) {
        keep_only_bands((nb_valence_bands > 0) ? nb_valence_bands : m_valence_bands.count,
                        (nb_conduction_bands > 0) ? nb_conduction_bands : m_conduction_bands.count);
    }
    auto_shift_conduction_band_energies();
    // auto_set_positive_valence_band_energies();
    compute_min_max_energies_at_tetras();
    compute_energy_gradient_at_tetras();
    set_energy_gradient_at_vertices_by_averaging_tetras();

    for (auto&& tetra : m_list_tetrahedra) {
        tetra.pre_compute_sorted_slots_per_band();
    }
    print_band_info();
    std::cout << "Done." << std::endl;
}

void MeshBZ::print_band_info() const {
    std::cout << "Band info:\n";
    std::cout << "Index\tType\tLocalIdx\tGlobalIdx\tMinE(eV)\tMaxE(eV)\n";
    for (std::size_t i = 0; i < m_band_info.size(); ++i) {
        const auto& band     = m_band_info[i];
        std::string type_str = (band.type == MeshParticleType::valence) ? "Valence   " : "Conduction";
        std::cout << i << "\t" << type_str << "\t" << band.local_index << "\t\t" << i << "\t\t" << std::fixed << std::setprecision(4)
                  << m_min_band[i] << "\t" << m_max_band[i] << "\n";
    }
}

void MeshBZ::precompute_dos_tetra(double energy_step, double energy_max) {
    std::cout << "Precomputing DOS per tetrahedra with energy step = " << energy_step << " eV ..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
#pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < m_list_tetrahedra.size(); ++i) {
        m_list_tetrahedra[i].precompute_dos_on_energy_grid_per_band(energy_step, energy_max);
    }
}

void MeshBZ::set_energy_gradient_at_vertices_by_averaging_tetras() {
    std::cout << "Setting energy gradient at vertices by volume-weighted averaging ..." << std::endl;

    constexpr double eps = 1e-18;

    for (std::size_t i = 0; i < m_list_vertices.size(); ++i) {
        const std::size_t nb_bands = m_list_vertices[i].get_number_bands();

        for (std::size_t b = 0; b < nb_bands; ++b) {
            vector3 accum(0.0, 0.0, 0.0);
            double  wsum = 0.0;

            // Loop over incident tetrahedra
            for (auto t_idx : m_vertex_to_tetrahedra[i]) {
                const auto&  T  = m_list_tetrahedra[t_idx];
                const double VT = std::abs(T.get_signed_volume());
                if (VT <= eps) {
                    continue;
                }
                const vector3 gT = T.get_gradient_energy_at_band(b);  // P1: constant per tet
                accum += VT * gT;
                wsum += VT;
            }

            const vector3 g_i = (wsum > 0.0) ? (accum / wsum) : vector3(0.0, 0.0, 0.0);
            m_list_vertices[i].push_back_energy_gradient_at_band(g_i);
        }
    }

    std::cout << "Done." << std::endl;
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
            if (energies[i] < m_min_band[i]) {
                m_min_band[i] = energies[i];
            }
            if (energies[i] > m_max_band[i]) {
                m_max_band[i] = energies[i];
            }
        }
    }
}

std::vector<std::size_t> MeshBZ::get_band_indices(MeshParticleType type) const {
    std::vector<std::size_t> indices;
    if (type == MeshParticleType::valence) {
        indices.resize(m_valence_bands.count);
        for (std::size_t i = 0; i < m_valence_bands.count; ++i) {
            indices[i] = m_valence_bands.global_start_index + i;
        }
    } else {
        indices.resize(m_conduction_bands.count);
        for (std::size_t i = 0; i < m_conduction_bands.count; ++i) {
            indices[i] = m_conduction_bands.global_start_index + i;
        }
    }
    return indices;
}

/**
 * @brief Keep only a subset of bands.
 * Must be called after reading the band energies from the .msh file and after shifting/setting valence
 * absolute energies and so on.
 * @param required_nb_bands Number of bands to keep (from the top for valence bands, from the bottom for conduction bands).
 */
void MeshBZ::keep_only_bands(std::size_t nb_valence_bands, std::size_t nb_conduction_bands) {
    if (nb_valence_bands + nb_conduction_bands >= get_number_bands_total()) {
        throw std::runtime_error("Cannot keep more bands than available.");
    }
    // Valence 
    std::vector<std::size_t> valence_indices = get_band_indices(MeshParticleType::valence);
    if (nb_valence_bands < valence_indices.size()) {
        int nb_valence_to_remove = valence_indices.size() - nb_valence_bands;
        std::cout << "Removing " << nb_valence_to_remove << " valence bands (keeping " << nb_valence_bands << ")\n";
        for (int i = 0; i < nb_valence_to_remove; ++i) {
            std::size_t band_index_to_remove = valence_indices[valence_indices.size() - 1 - i];
            for (auto&& vtx : m_list_vertices) {
                vtx.remove_band_energy(band_index_to_remove);
            }
            m_band_info.erase(m_band_info.begin() + band_index_to_remove);
            m_valence_bands.count--;
            m_conduction_bands.global_start_index--;
        }
    }
    // Conduction
    std::vector<std::size_t> conduction_indices = get_band_indices(MeshParticleType::conduction);
    if (nb_conduction_bands < conduction_indices.size()) {
        int nb_conduction_to_remove = conduction_indices.size() - nb_conduction_bands;
        std::cout << "Removing " << nb_conduction_to_remove << " conduction bands (keeping " << nb_conduction_bands << ")\n";
        for (int i = 0; i < nb_conduction_to_remove; ++i) {
            std::size_t band_index_to_remove = conduction_indices[conduction_indices.size() - 1 - i];
            for (auto&& vtx : m_list_vertices) {
                vtx.remove_band_energy(band_index_to_remove);
            }
            m_band_info.erase(m_band_info.begin() + band_index_to_remove);
            m_conduction_bands.count--;
        }
    }
    m_nb_bands_total = nb_valence_bands + nb_conduction_bands;
    // Recompute min/max band energies
    recompute_min_max_energies();
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
    std::vector<std::size_t> valence_band_indices = get_band_indices(MeshParticleType::valence);
    for (auto&& vtx : m_list_vertices) {
        for (auto idx_band : valence_band_indices) {
            double energy = vtx.get_energy_at_band(idx_band);
            if (energy < 0.0) {
                vtx.set_band_energy(idx_band, std::fabs(energy));
            }
        }
    }
}

void MeshBZ::auto_shift_conduction_band_energies() {

    constexpr double band_reference_eV = 0.0;  // Set the conduction band minimum to 0 eV
    // Find the lowest conduction band minimum
    std::vector<std::size_t> conduction_band_indices = get_band_indices(MeshParticleType::conduction);
    if (conduction_band_indices.empty()) {
        std::cout << "No conduction bands found. Skipping auto-shift of conduction band energies.\n";
        return;
    }
    double min_conduction = std::numeric_limits<double>::max();
    for (auto idx_band : conduction_band_indices) {
        if (m_min_band[idx_band] < min_conduction) {
            min_conduction = m_min_band[idx_band];
        }
    }
    std::cout << "Min conduction band energy: " << min_conduction << " eV\n";

    double band_gap = min_conduction - band_reference_eV;
    std::cout << "Computed band gap: " << band_gap << " eV\n";

    if (band_gap < 0.1) {
        std::cout << "Warning: computed band gap is very small or negative (" << band_gap << " eV). No shift applied.\n";
        return;
    }

    // Shift all conduction bands so that the conduction band minimum aligns with the valence band maximum
    double shift_amount = -band_gap;
    std::cout << "Shifting conduction bands by " << shift_amount << " eV to align CBM with VBM.\n";

    for (auto&& vtx : m_list_vertices) {
        for (auto idx_band : conduction_band_indices) {
            double old_energy = vtx.get_energy_at_band(idx_band);
            double new_energy = old_energy + shift_amount;
            vtx.set_band_energy(idx_band, new_energy);
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

    // Reverse the order of valence bands only.
    std::vector<std::size_t> valence_indices = get_band_indices(MeshParticleType::valence);
    // Check if already in the right order
    bool already_in_order = true;
    for (std::size_t i = 1; i < valence_indices.size(); ++i) {
        if (valence_indices[i] < valence_indices[i - 1]) {
            already_in_order = false;
            break;
        }
    }
    if (already_in_order) {
        std::cout << "Valence bands already in the right order. No change applied.\n";
        return;
    }
    // void reverse_energies_order(int first_idx, int last_idx) {
    int global_idx_first = valence_indices.front();
    int global_idx_last  = valence_indices.back();
    for (auto&& vtx : m_list_vertices) {
        vtx.reverse_energies_order(global_idx_first, global_idx_last);
    }

    recompute_min_max_energies();
}

void MeshBZ::compute_energy_gradient_at_tetras() {
    for (auto&& tetra : m_list_tetrahedra) {
        tetra.compute_gradient_energy_at_bands();
    }
}

vector3 MeshBZ::interpolate_energy_gradient_at_location(const vector3& location, const std::size_t& idx_band) const {
    Tetra* tetra = find_tetra_at_location(location);
    if (tetra == nullptr) {
        throw std::runtime_error("Location is outside the mesh. Cannot interpolate energy gradient.");
    }
    return tetra->get_gradient_energy_at_band(idx_band);
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

    std::vector<double> list_energies(nb_points);
    std::vector<double> list_dos(nb_points);
#pragma omp parallel for schedule(dynamic) num_threads(num_threads)
    for (std::size_t index_energy = 0; index_energy < nb_points; ++index_energy) {
        double energy = min_energy + index_energy * energy_step;
        double dos    = compute_dos_at_energy_and_band(energy, band_index, use_interp);
        list_energies[index_energy] = energy;
        list_dos[index_energy] = dos;
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
    std::cout << "Computing DOS weights for tetrahedra at energy " << energy << " eV ..." << std::endl;
    for (auto&& tetra : m_list_tetrahedra) {
        std::cout << "\rTetra " << tetra.get_index() << "/" << m_list_tetrahedra.size() << std::flush;
        double dos = tetra.compute_tetra_dos_energy_band(energy, idx_band);
        list_dos.push_back(dos);
    }
    std::cout << "\nDrawing tetrahedron with DOS weights ..." << std::endl;
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
    std::cout << "Selected tetrahedron index: " << index_tetra << std::endl;
    if (index_tetra >= m_list_tetrahedra.size()) {
        throw std::runtime_error("Selected tetrahedron index is out of range");
    }
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
                    if (n1 == 0 && n2 == 0 && n3 == 0) {
                        continue;
                    }
                    if (std::abs(n1) + std::abs(n2) + std::abs(n3) == L1) {
                        m_Gshifts.push_back(n1 * b1 + n2 * b2 + n3 * b3);
                    }
                }
            }
        }
    }
}

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

// --- O(1) WS-BZ folding (public API) ---
/**
 * @brief Initialize the reciprocal basis and reduced-frame parameters for O(1) folding.
 * @param b1_SI First reciprocal primitive vector (SI, 1/m)
 * @param b2_SI Second reciprocal primitive vector (SI, 1/m)
 * @param b3_SI Third reciprocal primitive vector (SI, 1/m)
 * @param halfwidth_reduced Half-width in reduced coords (0.5 for [-0.5,0.5])
 * @param si_to_reduced Scale factor from SI (1/m) to reduced coords used by plane tests
 */
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

/**
 * @brief Fold k into the Wigner–Seitz BZ (truncated octahedron of bcc). Pure, no side effects.
 * @param k_SI Input wavevector (SI, 1/m)
 * @return Folded wavevector inside the first BZ (SI, 1/m)
 */
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
            if (idx == 0) {
                kr.x() -= (kr.x() > 0 ? 1.0 : -1.0) * excess;
            } else if (idx == 1) {
                kr.y() -= (kr.y() > 0 ? 1.0 : -1.0) * excess;
            } else {
                kr.z() -= (kr.z() > 0 ? 1.0 : -1.0) * excess;
            }
            moved = true;
        }

        if (!moved) {
            break;  // inside WS
        }
    }

    // Back to SI and convert to your vector3
    Eigen::Vector3d kf = kr / m_si2red;
    return vector3{kf.x(), kf.y(), kf.z()};
}

/**
 * @brief Check if k is inside the WS BZ using the same plane tests as is_inside_mesh_geometry().
 * @param k_SI Input wavevector (SI, 1/m)
 * @return true if inside, false otherwise
 */
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
 * @brief Find the index of the vertex in the irreducible wedge that represents k_SI.
 * In this version, we first first fold k_SI into the first IW, then search for the closest vertex.
 *
 *
 * @param k_SI
 * @return std::size_t
 */
std::size_t MeshBZ::get_index_irreducible_wedge(const vector3& k_SI) const {
    vector3 k_folded = {std::fabs(k_SI.x()), std::fabs(k_SI.y()), std::fabs(k_SI.z())};
    // Test the 6 permutations of (|kx|, |ky|, |kz|)
    std::array<vector3, 6> permutations = {vector3{k_folded.x(), k_folded.y(), k_folded.z()},
                                           vector3{k_folded.x(), k_folded.z(), k_folded.y()},
                                           vector3{k_folded.y(), k_folded.x(), k_folded.z()},
                                           vector3{k_folded.y(), k_folded.z(), k_folded.x()},
                                           vector3{k_folded.z(), k_folded.x(), k_folded.y()},
                                           vector3{k_folded.z(), k_folded.y(), k_folded.x()}};
    bool                   found        = false;
    for (auto&& perm : permutations) {
        if (is_irreducible_wedge(perm)) {
            k_folded = perm;
            found    = true;
            break;
        }
    }
    if (!found) {
        std::cout << "Could not fold k = " << k_SI << " into the irreducible wedge." << std::endl;
        throw std::runtime_error("Could not fold k into the irreducible wedge.");
    }
    // Now search for the closest vertex in the IW
    double       min_dist = std::numeric_limits<double>::max();
    std::size_t  idx_min  = 0;
    const double s        = si_to_reduced_scale();
    for (const auto& vtx : m_list_vertices) {
        if (is_irreducible_wedge(vtx.get_position())) {
            double dist = (vtx.get_position() - k_folded).norm_squared();
            if (dist < min_dist) {
                min_dist = dist;
                idx_min  = vtx.get_index();
            }
        }
    }
    return idx_min;
}

static inline void bz_write_vtk_scalars(std::ofstream&             out,
                                        const std::string&         name,
                                        const std::vector<double>& vals,
                                        const char*                loc_keyword,  // "POINT_DATA" or "CELL_DATA" already emitted
                                        std::size_t                expected_count) {
    if (vals.size() != expected_count) {
        throw std::runtime_error("VTK export: scalar field '" + name + "' has size " + std::to_string(vals.size()) + ", expected " +
                                 std::to_string(expected_count));
    }
    out << "SCALARS " << name << " double 1\n";
    out << "LOOKUP_TABLE default\n";
    out << std::setprecision(17);
    for (double v : vals) {
        out << v << "\n";
    }
}

static inline void bz_write_vtk_vectors(std::ofstream&              out,
                                        const std::string&          name,
                                        const std::vector<vector3>& vecs,
                                        const char*                 loc_keyword,
                                        std::size_t                 expected_count) {
    if (vecs.size() != expected_count) {
        throw std::runtime_error("VTK export: vector field '" + name + "' has size " + std::to_string(vecs.size()) + ", expected " +
                                 std::to_string(expected_count));
    }
    out << "VECTORS " << name << " double\n";
    out << std::setprecision(17);
    for (const auto& v : vecs) {
        out << v.x() << " " << v.y() << " " << v.z() << "\n";
    }
}

void MeshBZ::export_to_vtk(const std::string&        filename,
                           const MapStringToDoubles& point_scalars,
                           const MapStringToVectors& point_vectors,
                           const MapStringToDoubles& cell_scalars,
                           const MapStringToVectors& cell_vectors) const {
    const std::size_t n_points = m_list_vertices.size();
    const std::size_t n_cells  = m_list_tetrahedra.size();

    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Cannot open VTK file '" + filename + "' for writing.");
    }

    // --- VTK legacy header (ASCII, UnstructuredGrid) ---
    out << "# vtk DataFile Version 4.2\n";
    out << "BZ mesh export\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // --- POINTS ---
    out << "POINTS " << n_points << " double\n";
    out << std::setprecision(17);
    for (const auto& v : m_list_vertices) {
        const auto& p = v.get_position();
        out << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    // --- CELLS (tetra) ---
    // For legacy VTK: total_indices = n_cells * (1 + 4) because each line: "4 i0 i1 i2 i3"
    const std::size_t total_idx = n_cells * (1 + 4);
    out << "CELLS " << n_cells << " " << total_idx << "\n";
    for (const auto& t : m_list_tetrahedra) {
        const auto ids = t.get_list_indices_vertices();  // assumes 0-based vertex ids
        out << "4 " << ids[0] << " " << ids[1] << " " << ids[2] << " " << ids[3] << "\n";
    }

    // --- CELL_TYPES (tetra = 10) ---
    out << "CELL_TYPES " << n_cells << "\n";
    for (std::size_t i = 0; i < n_cells; ++i) {
        out << "10\n";
    }

    // --- Optional per-vertex data ---
    if (!point_scalars.empty() || !point_vectors.empty()) {
        out << "POINT_DATA " << n_points << "\n";
        for (const auto& [name, vals] : point_scalars) {
            bz_write_vtk_scalars(out, name, vals, "POINT_DATA", n_points);
        }
        for (const auto& [name, vecs] : point_vectors) {
            bz_write_vtk_vectors(out, name, vecs, "POINT_DATA", n_points);
        }
    }

    // --- Optional per-cell data ---
    if (!cell_scalars.empty() || !cell_vectors.empty()) {
        out << "CELL_DATA " << n_cells << "\n";
        for (const auto& [name, vals] : cell_scalars) {
            bz_write_vtk_scalars(out, name, vals, "CELL_DATA", n_cells);
        }
        for (const auto& [name, vecs] : cell_vectors) {
            bz_write_vtk_vectors(out, name, vecs, "CELL_DATA", n_cells);
        }
    }
    out.close();
}

void MeshBZ::export_energies_and_gradients_to_vtk(const std::string& filename) const {
    std::map<std::string, std::vector<double>>  point_scalars;
    std::map<std::string, std::vector<vector3>> point_vectors;

    // Per-vertex band energies
    for (std::size_t b = 0; b < m_list_vertices[0].get_number_bands(); ++b) {
        std::string         name = "band_energy_" + std::to_string(b);
        std::vector<double> vals;
        vals.reserve(m_list_vertices.size());
        for (const auto& v : m_list_vertices) {
            vals.push_back(v.get_energy_at_band(b));
        }
        point_scalars[name] = vals;
    }

    // Per-vertex band energy gradients
    for (std::size_t b = 0; b < m_list_vertices[0].get_energy_gradient_at_bands().size(); ++b) {
        std::string          name = "band_grad_" + std::to_string(b);
        std::vector<vector3> vecs;
        vecs.reserve(m_list_vertices.size());
        for (const auto& v : m_list_vertices) {
            vecs.push_back(v.get_energy_gradient_at_band(b));
        }
        point_vectors[name] = vecs;
    }
    export_to_vtk(filename, point_scalars, point_vectors, {}, {});
}

void MeshBZ::export_octree_to_vtu(const std::string& filename) const {
    if (!m_search_tree) {
        throw std::runtime_error("Octree not initialized. Cannot export.");
    }
    bool export_leaves_only = true;
    write_octree_as_vtu(*m_search_tree, filename, export_leaves_only);
}

}  // namespace uepm::mesh_bz