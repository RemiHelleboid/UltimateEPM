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
    m_node_tags = nodeTags;
     std::cout << "Number of nodes: " << size_nodes_tags << std::endl;

    if (size_nodes_coordinates != 3 * size_nodes_tags) {
        throw std::runtime_error("Number of coordinates is not 3 times the number of vertices. Abort.");
    }

    m_list_vertices.reserve(size_nodes_tags);
    double lattice_constant = m_material.get_lattice_constant_meter();
    std::cout << "Lattice const: " << lattice_constant << std::endl;
    std::cout << "V: " << std::pow(2.0 * M_PI, 3) / std::pow(lattice_constant, 3.0) << std::endl;
    const double fourier_factor       = 2.0 * M_PI / lattice_constant;
    // const double fourier_factor       = 1;
    double       normalization_factor = normalize_by_fourier_factor ? fourier_factor : 1.0;
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
void MeshBZ::read_mesh_bands_from_msh_file(const std::string& filename) {
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
    compute_min_max_energies_at_tetras();
    compute_energy_gradient_at_tetras();
}

void MeshBZ::read_mesh_bands_from_multi_band_files(const std::string& dir_bands, int nb_bands_to_load){

    if (m_list_vertices.empty()) throw std::runtime_error("Geometry must be loaded before bands (m_list_vertices is empty).");

    if (m_node_tags.empty()) throw std::runtime_error("m_node_tags is empty. Store node tags when reading geometry to preserve mapping.");

    // 1) Collect band files: band_<idx>.msh
    std::vector<std::pair<int, std::filesystem::path>> band_files;
    std::regex                            re(R"(band_(\d+)\.msh$)");
    for (auto& p : std::filesystem::directory_iterator(dir_bands)) {
        if (!p.is_regular_file()) continue;
        const std::string name = p.path().filename().string();
        std::smatch  m;
        if (std::regex_search(name, m, re)) {
            int idx = std::stoi(m[1].str());
            band_files.emplace_back(idx, p.path());
        }
    }
    if (band_files.empty()) throw std::runtime_error("No band_*.msh files found in directory: " + dir_bands);

    std::sort(band_files.begin(), band_files.end(), [](auto& a, auto& b) { return a.first < b.first; });

    // 2) Build tag->vertex index map once (reference from geometry session)
    std::unordered_map<std::size_t, std::size_t> tag2idx;
    tag2idx.reserve(m_node_tags.size());
    for (std::size_t i = 0; i < m_node_tags.size(); ++i)
        tag2idx[m_node_tags[i]] = i;

    // 3) Clear existing band metadata (if reloading)
    m_indices_valence_bands.clear();
    m_indices_conduction_bands.clear();
    m_min_band.clear();
    m_max_band.clear();

    // 4) Loop over band files and ingest view values
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 0);

    int count_band = 0;
    for (auto& [bidx, path] : band_files) {
        if (nb_bands_to_load > 0 && count_band >= nb_bands_to_load) {
            std::cout << "Loaded " << count_band << " bands, stopping as requested." << std::endl;
            break;
        }
        std::cout << "Opening band file: " << path << std::endl;

        gmsh::clear();
        gmsh::open(path.string());

        // Expect exactly one view
        std::vector<int> vtags;
        gmsh::view::getTags(vtags);
        if (vtags.empty()) throw std::runtime_error("No views found in " + path.string());
        if (vtags.size() > 1) std::cerr << "[warn] " << path << " contains " << vtags.size() << " views; using the first one.\n";

        const int vtag = vtags.front();

        std::string              type;
        std::vector<std::size_t> tags;
        double                   time;
        int                      numComp;
        std::vector<double>      data_view;

        gmsh::view::getHomogeneousModelData(vtag, 0, type, tags, data_view, time, numComp);

        if (type != "NodeData") throw std::runtime_error("Unexpected view type in " + path.string() + ": " + type);
        if (numComp != 1) throw std::runtime_error("Non-scalar NodeData in " + path.string());
        if (tags.size() != data_view.size()) throw std::runtime_error("tags.size() != data_view.size() in " + path.string());

        if (tags.size() != m_list_vertices.size()) {
            std::ostringstream oss;
            oss << "Node count mismatch in " << path << ": file has " << tags.size() << " nodes, mesh has " << m_list_vertices.size();
            throw std::runtime_error(oss.str());
        }

        // Map values to our vertex order via node tag
        std::vector<double> energies_in_vertex_order(m_list_vertices.size(), 0.0);
        for (std::size_t i = 0; i < tags.size(); ++i) {
            auto it = tag2idx.find(tags[i]);
            if (it == tag2idx.end()) {
                std::ostringstream oss;
                oss << "Node tag " << tags[i] << " from " << path << " not present in reference tag map.";
                throw std::runtime_error(oss.str());
            }
            energies_in_vertex_order[it->second] = data_view[i];
            // HARD CODED SHIFT FOR CONDUCTION BANDS - DEBUG
            double band_gap_si = 1.05;
            if (energies_in_vertex_order[it->second] > 0.5) {
                energies_in_vertex_order[it->second] -= band_gap_si;
            }
        }

        

        // Attach this band to vertices
        add_new_band_energies_to_vertices(energies_in_vertex_order);

        // Classify + min/max for this band
        auto minmax_band = std::minmax_element(energies_in_vertex_order.begin(), energies_in_vertex_order.end());
        m_min_band.push_back(*minmax_band.first);
        m_max_band.push_back(*minmax_band.second);

        // Adjust criterion if your energy zero differs
        const bool is_valence = (*minmax_band.second) <= 0.0;
        if (is_valence)
            m_indices_valence_bands.push_back(count_band);
        else
            m_indices_conduction_bands.push_back(count_band);

        ++count_band;
    }

    gmsh::finalize();

    // 5) Recompute per-tetra helpers
    compute_min_max_energies_at_tetras();
    compute_energy_gradient_at_tetras();

    std::cout << "Loaded " << count_band << " bands from " << band_files.size() << " files in directory: " << dir_bands << std::endl;
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

double MeshBZ::compute_dos_at_energy_and_band(double iso_energy, int band_index) const {
    double total_dos = 0.0;
    for (auto&& tetra : m_list_tetrahedra) {
        total_dos += tetra.compute_tetra_dos_energy_band(iso_energy, band_index);
    }
    return total_dos;
}

std::vector<std::vector<double>> MeshBZ::compute_dos_band_at_band(int         band_index,
                                                                  double      min_energy,
                                                                  double      max_energy,
                                                                  int         num_threads,
                                                                  std::size_t nb_points) const {
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

std::vector<std::vector<double>> MeshBZ::compute_dos_band_at_band_auto(int band_index, std::size_t nb_points, int num_threads) const {
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
        double dos    = compute_dos_at_energy_and_band(energy, band_index);
        list_energies.push_back(energy);
        list_dos.push_back(dos);
    }
    auto end              = std::chrono::high_resolution_clock::now();
    auto total_time_count = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nDOS for 1 band computed in  " << total_time_count / 1000.0 << "s" << std::endl;
    return {list_energies, list_dos};
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

inline double MeshBZ::si_to_reduced_scale() const {
    // reduced k = (a / (2π)) * k_SI
    return m_material.get_lattice_constant_meter() / (2.0 * M_PI);
}

// helper: half-width of your cube in *reduced* coords.
// If your mesh spans [-1,1], half-width = 1.0 (not 0.5).
static inline double bz_halfwidth_reduced() { return 1.0; }  // <- your case

bool MeshBZ::is_inside_mesh_geometry(const vector3& k) const {
    const double s  = si_to_reduced_scale();   // converts SI (1/m) -> reduced (unitless)
    const double hw = bz_halfwidth_reduced();  // 1.0 for [-1,1], 0.5 for [-0.5,0.5]

    const double kx = k.x() * s, ky = k.y() * s, kz = k.z() * s;  // reduced coords
    const double ax = std::abs(kx), ay = std::abs(ky), az = std::abs(kz);

    // cube bound
    const bool cond1 = (ax <= hw) && (ay <= hw) && (az <= hw);

    // truncation planes (fcc BZ): |kx|+|ky|+|kz| <= (3/2)*hw
    const bool cond2 = (ax + ay + az) <= (1.5 * hw);

    return cond1 && cond2;
}

bool MeshBZ::is_inside_mesh_geometry(const Vector3D<double>& k) const {
    const double s  = si_to_reduced_scale();
    const double hw = bz_halfwidth_reduced();

    const double kx = k.X * s, ky = k.Y * s, kz = k.Z * s;
    const double ax = std::abs(kx), ay = std::abs(ky), az = std::abs(kz);

    const bool cond1 = (ax <= hw) && (ay <= hw) && (az <= hw);
    const bool cond2 = (ax + ay + az) <= (1.5 * hw);

    return cond1 && cond2;
}

vector3 MeshBZ::retrieve_k_inside_mesh_geometry(const vector3& k) const {
    vector3 b1 = {-1.0, 1.0, 1.0};
    vector3 b2 = {1.0, -1.0, 1.0};
    vector3 b3 = {1.0, 1.0, -1.0};
    b1 /= si_to_reduced_scale();
    b2 /= si_to_reduced_scale();
    b3 /= si_to_reduced_scale();

    std::vector<int> list_n_k = {0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5};
    for (auto&& n_k_x : list_n_k) {
        for (auto&& n_k_y : list_n_k) {
            for (auto&& n_k_z : list_n_k) {
                vector3 k_plus_G = k + n_k_x * b1 + n_k_y * b2 + n_k_z * b3;
                if (is_inside_mesh_geometry(k_plus_G)) {
                    // std::cout << "k: " << k << std::endl;
                    // std::cout << " G: " << n_k_x * b1 + n_k_y * b2 + n_k_z * b3 << std::endl;
                    // std::cout << "k + G: " << k_plus_G * si_to_reduced_scale() << std::endl;
                    return k_plus_G;
                }
            }
        }
    }
    std::cout << "No k-point inside the mesh geometry found for k: " << k  << std::endl;
    throw std::runtime_error("No k-point inside the mesh geometry found.");
}

}  // namespace bz_mesh