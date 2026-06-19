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

#include <fmt/chrono.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "gmsh.h"
#include "gmsh_guard.hpp"
#include "integrals.hpp"
#include "numerical_helper.hpp"
#include "octree_bz.hpp"
#include "omp.h"
#include "physical_constants.hpp"
#include "rapidcsv.h"
#include "string_utils.hpp"
#include "vector_bz.hpp"

#pragma omp declare reduction(merge : std::vector<double> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

namespace uepm::mesh_bz {

namespace {

constexpr std::size_t invalid_vertex_index = std::numeric_limits<std::size_t>::max();

struct ReducedCoordinateKey {
    std::int64_t x;
    std::int64_t y;
    std::int64_t z;

    bool operator==(const ReducedCoordinateKey&) const = default;
};

struct ReducedCoordinateKeyHash {
    std::size_t operator()(const ReducedCoordinateKey& key) const noexcept {
        std::size_t seed = std::hash<std::int64_t>{}(key.x);
        seed ^= std::hash<std::int64_t>{}(key.y) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<std::int64_t>{}(key.z) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

}  // namespace

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

vector3 MeshBZ::si_to_reduced_k(const vector3& k_si) const noexcept { return k_si * si_to_reduced_scale(); }

vector3 MeshBZ::reduced_to_si_k(const vector3& k_reduced) const noexcept { return k_reduced / si_to_reduced_scale(); }

/**
 * @brief Read the geometry of the mesh from the .msh file: the vertices and the elements are added to
 * the m_list_vertices and m_list_elements lists.
 * All the points coordinates are re-normalized by the lattice constant passed as argument.
 *
 * @param filename
 * @param lattice_constant
 */
void MeshBZ::read_mesh_geometry_from_msh_file(const std::string& filename, bool input_coordinates_are_reduced) {
    m_filename_mesh = filename;
    m_search_tree.reset();
    m_list_vtx_in_iwedge.clear();
    m_kstar_ibz_to_bz.clear();
    std::cout << "Opening file " << filename << std::endl;
    GmshSession gmsh_session;
    gmsh::open(filename);
    std::vector<std::size_t> nodeTags;
    std::vector<double>      nodeCoords;
    std::vector<double>      nodeParams;
    std::cout << "Reading vertices ..." << std::endl;
    // gmsh::model::mesh::reclassifyNodes();
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, -1, -1, false, false);
    std::size_t size_nodes_tags        = nodeTags.size();
    std::size_t size_nodes_coordinates = nodeCoords.size();
    m_node_tags                        = nodeTags;
    m_source_vertex_indices.resize(size_nodes_tags);
    std::iota(m_source_vertex_indices.begin(), m_source_vertex_indices.end(), 0);
    m_local_index_from_source_vertex = m_source_vertex_indices;

    std::cerr << "sizeof(vector3)=" << sizeof(vector3) << " sizeof(bbox_mesh)=" << sizeof(bbox_mesh)
              << " sizeof(UniformDos)=" << sizeof(UniformDos) << " sizeof(Tetra)=" << sizeof(Tetra) << "\n";
    std::cout << "Number of nodes: " << size_nodes_tags << std::endl;

    if (size_nodes_coordinates != 3 * size_nodes_tags) {
        throw std::runtime_error("Number of coordinates is not 3 times the number of vertices. Abort.");
    }

    m_list_vertices.resize(size_nodes_tags);
    double lattice_constant = m_material.get_lattice_constant_meter();
    std::cout << "Lattice const: " << lattice_constant << std::endl;
    std::cout << "V: " << std::pow(2.0 * M_PI, 3) / std::pow(lattice_constant, 3.0) << std::endl;
    const double fourier_factor = 2.0 * M_PI / lattice_constant;
    // const double fourier_factor       = 1;
    const double normalization_factor = input_coordinates_are_reduced ? fourier_factor : 1.0;
    std::cout << "Nb of threads for mesh ops: " << m_nb_threads_mesh_ops << std::endl;
#pragma omp parallel for schedule(static) num_threads(m_nb_threads_mesh_ops)
    for (std::size_t index_vertex = 0; index_vertex < size_nodes_tags; ++index_vertex) {
        m_list_vertices[index_vertex] = Vertex(index_vertex,
                                               normalization_factor * nodeCoords[3 * index_vertex],
                                               normalization_factor * nodeCoords[3 * index_vertex + 1],
                                               normalization_factor * nodeCoords[3 * index_vertex + 2]);
    }
    std::cout << "Number of k-points vertices: " << m_list_vertices.size() << std::endl;

    std::unordered_map<std::size_t, std::size_t> local_index_from_node_tag;
    local_index_from_node_tag.reserve(nodeTags.size());
    for (std::size_t local_index = 0; local_index < nodeTags.size(); ++local_index) {
        const auto [unused, inserted] = local_index_from_node_tag.emplace(nodeTags[local_index], local_index);
        if (!inserted) {
            throw std::runtime_error("Duplicate Gmsh node tag " + std::to_string(nodeTags[local_index]) + ".");
        }
    }
    const auto local_index_from_tag = [&](std::size_t node_tag) {
        const auto it = local_index_from_node_tag.find(node_tag);
        if (it == local_index_from_node_tag.end()) {
            throw std::runtime_error("A tetrahedron references unknown Gmsh node tag " + std::to_string(node_tag) +
                                     ".");
        }
        return it->second;
    };

    // Get the mesh elements for the entity (dim, tag):
    std::cout << "Reading elements ..." << std::endl;
    const int                             dim = 3;
    const int                             tag = -1;
    std::vector<int>                      elemTypes;
    std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);
    if (elemTags.empty()) {
        std::cout << "ElementTags is zero when the mesh was imported... Abort.\n";
        throw std::runtime_error("ElementTags is zero when the mesh was imported... Abort.");
    }
    std::size_t number_elements = 0;
    for (std::size_t block_index = 0; block_index < elemTypes.size(); ++block_index) {
        if (elemTypes[block_index] != 4) {
            continue;
        }
        if (elemNodeTags[block_index].size() != 4 * elemTags[block_index].size()) {
            throw std::runtime_error("Invalid connectivity size in a linear tetrahedron element block.");
        }
        number_elements += elemTags[block_index].size();
    }
    if (number_elements == 0) {
        throw std::runtime_error("The mesh contains no linear tetrahedron elements.");
    }
    std::cout << "Number of elements: " << number_elements << " ... reserving memory." << std::endl;
    std::cout << "Size of tetrahedra: " << sizeof(Tetra) << " bytes." << std::endl;

    m_list_tetrahedra.clear();
    m_list_tetrahedra.reserve(number_elements);
    m_vertex_to_tetrahedra.clear();
    m_vertex_to_tetrahedra.resize(m_list_vertices.size());

    std::size_t      degenerate_count = 0;
    constexpr double vol_threshold    = 1e-12;

    std::size_t source_element_index = 0;
    for (std::size_t block_index = 0; block_index < elemTypes.size(); ++block_index) {
        if (elemTypes[block_index] != 4) {
            continue;
        }
        const auto& connectivity = elemNodeTags[block_index];
        for (std::size_t index_element = 0; index_element < elemTags[block_index].size();
             ++index_element, ++source_element_index) {
            std::array<std::size_t, 4> vertex_indices{};
            std::array<Vertex*, 4>     array_element_vertices{};
            for (std::size_t local_vertex = 0; local_vertex < 4; ++local_vertex) {
                vertex_indices[local_vertex] = local_index_from_tag(connectivity[4 * index_element + local_vertex]);
                array_element_vertices[local_vertex] = &m_list_vertices[vertex_indices[local_vertex]];
            }

            Tetra        new_tetra(0, array_element_vertices);
            const double signed_vol = new_tetra.get_signed_volume();

            if (std::abs(signed_vol) < vol_threshold) {
                std::cerr << "Warning: Tetrahedron " << source_element_index
                          << " has a very small volume (|6*V| = " << std::abs(signed_vol)
                          << "). This may lead to numerical instability in barycentric coordinate computation."
                          << std::endl;
                ++degenerate_count;
                continue;
            }

            const bool in_iwedge = is_irreducible_wedge(new_tetra.get_barycenter());
            new_tetra.set_lies_in_irreducible_wedge(in_iwedge);

            const std::size_t compact_index = m_list_tetrahedra.size();
            new_tetra.set_index(compact_index);

            m_list_tetrahedra.push_back(std::move(new_tetra));

            for (const std::size_t vertex_index : vertex_indices) {
                m_vertex_to_tetrahedra[vertex_index].push_back(compact_index);
            }
        }
    }

    m_list_tetrahedra.shrink_to_fit();

    if (stores_positive_octant()) {
        compact_geometry_to_positive_octant();
    }

    fmt::print("Number of degenerate tetrahedra: {} / {} ({:.2f}%)\n",
               degenerate_count,
               number_elements,
               100.0 * static_cast<double>(degenerate_count) / static_cast<double>(number_elements));

    m_total_volume = compute_mesh_volume();
    std::cout << "Total mesh volume: " << m_total_volume << std::endl;
    // Compute the reduced BZ volume
    double       Vcell      = std::pow(m_material.get_lattice_constant_meter(), 3) / 4.0;
    const double VBZ_theory = std::pow(2.0 * M_PI, 3) / Vcell;  // m^-3

    if (!(m_total_volume > 0.0)) {
        throw std::runtime_error("The tetrahedral mesh has zero total volume.");
    }
    m_bz_volume_correction = VBZ_theory / m_total_volume;
    std::cout << "BZ volume correction: " << m_bz_volume_correction << std::endl;

    precompute_G_shifts();
    const double     reduced_to_si        = 1.0 / si_to_reduced_scale();
    Eigen::Vector3d  b1_SI                = reduced_to_si * Eigen::Vector3d{-1.0, 1.0, 1.0};
    Eigen::Vector3d  b2_SI                = reduced_to_si * Eigen::Vector3d{1.0, -1.0, 1.0};
    Eigen::Vector3d  b3_SI                = reduced_to_si * Eigen::Vector3d{1.0, 1.0, -1.0};
    constexpr double halfwidth_reduced    = 1.0;
    const double     ssi_to_reduced_scale = si_to_reduced_scale();
    init_reciprocal_basis(b1_SI, b2_SI, b3_SI, halfwidth_reduced, ssi_to_reduced_scale);

    if (stores_positive_octant()) {
        build_positive_octant_kstar();
    } else {
        load_kstar_ibz_to_bz();
    }
}

void MeshBZ::compact_geometry_to_positive_octant() {
    constexpr double tolerance = 1e-12;

    std::vector<std::array<std::size_t, 4>> retained_source_tetrahedra;
    retained_source_tetrahedra.reserve(m_list_tetrahedra.size() / 8 + 1);
    std::vector<bool> source_vertex_is_used(m_list_vertices.size(), false);

    for (const Tetra& tetra : m_list_tetrahedra) {
        const auto source_indices = tetra.get_list_indices_vertices();
        bool       inside_octant  = true;
        for (std::size_t source_index : source_indices) {
            const vector3& position = m_list_vertices[source_index].get_position();
            inside_octant = inside_octant && position.x() >= -tolerance && position.y() >= -tolerance &&
                            position.z() >= -tolerance;
        }
        if (!inside_octant) {
            continue;
        }
        retained_source_tetrahedra.push_back(source_indices);
        for (std::size_t source_index : source_indices) {
            source_vertex_is_used[source_index] = true;
        }
    }

    std::vector<std::size_t> local_from_source(m_list_vertices.size(), invalid_vertex_index);
    std::vector<Vertex>      compact_vertices;
    std::vector<std::size_t> compact_node_tags;
    std::vector<std::size_t> compact_source_indices;
    compact_vertices.reserve(std::count(source_vertex_is_used.begin(), source_vertex_is_used.end(), true));
    compact_node_tags.reserve(compact_vertices.capacity());
    compact_source_indices.reserve(compact_vertices.capacity());

    for (std::size_t source_index = 0; source_index < m_list_vertices.size(); ++source_index) {
        if (!source_vertex_is_used[source_index]) {
            continue;
        }
        const std::size_t local_index = compact_vertices.size();
        local_from_source[source_index] = local_index;
        compact_vertices.emplace_back(local_index, m_list_vertices[source_index].get_position());
        compact_node_tags.push_back(m_node_tags[source_index]);
        compact_source_indices.push_back(source_index);
    }

    std::vector<Tetra> compact_tetrahedra;
    compact_tetrahedra.reserve(retained_source_tetrahedra.size());
    std::vector<std::vector<std::size_t>> compact_vertex_to_tetrahedra(compact_vertices.size());
    for (const auto& source_indices : retained_source_tetrahedra) {
        std::array<Vertex*, 4> vertices{};
        for (std::size_t slot = 0; slot < vertices.size(); ++slot) {
            vertices[slot] = &compact_vertices.at(local_from_source.at(source_indices[slot]));
        }
        const std::size_t tetra_index = compact_tetrahedra.size();
        compact_tetrahedra.emplace_back(tetra_index, vertices);
        compact_tetrahedra.back().set_lies_in_irreducible_wedge(
            is_irreducible_wedge(compact_tetrahedra.back().get_barycenter()));
        for (Vertex* vertex : vertices) {
            compact_vertex_to_tetrahedra[vertex->get_index()].push_back(tetra_index);
        }
    }

    if (compact_vertices.empty() || compact_tetrahedra.empty()) {
        throw std::runtime_error("Positive-octant reduction produced an empty mesh");
    }

    m_list_tetrahedra.clear();
    m_list_vertices.clear();
    m_list_vertices            = std::move(compact_vertices);
    m_list_tetrahedra          = std::move(compact_tetrahedra);
    m_vertex_to_tetrahedra     = std::move(compact_vertex_to_tetrahedra);
    m_node_tags                = std::move(compact_node_tags);
    m_source_vertex_indices    = std::move(compact_source_indices);
    m_local_index_from_source_vertex = std::move(local_from_source);

    fmt::print("Compacted BZ storage to positive octant: {} vertices, {} tetrahedra\n",
               m_list_vertices.size(),
               m_list_tetrahedra.size());
}

void MeshBZ::build_positive_octant_kstar() {
    if (!stores_positive_octant()) {
        throw std::logic_error("Positive-octant symmetry orbits require positive-octant storage");
    }

    constexpr double coordinate_tolerance_reduced = 1e-10;
    const auto coordinate_key = [&](const vector3& position) {
        const vector3 reduced = si_to_reduced_k(position);
        return ReducedCoordinateKey{
            std::llround(reduced.x() / coordinate_tolerance_reduced),
            std::llround(reduced.y() / coordinate_tolerance_reduced),
            std::llround(reduced.z() / coordinate_tolerance_reduced),
        };
    };

    std::unordered_map<ReducedCoordinateKey, std::size_t, ReducedCoordinateKeyHash> vertex_by_position;
    vertex_by_position.reserve(m_list_vertices.size());
    for (const Vertex& vertex : m_list_vertices) {
        const auto [it, inserted] = vertex_by_position.emplace(coordinate_key(vertex.get_position()), vertex.get_index());
        if (!inserted) {
            throw std::runtime_error("Positive-octant mesh contains duplicate vertices within symmetry tolerance");
        }
    }

    m_list_vtx_in_iwedge.clear();
    m_kstar_ibz_to_bz.assign(m_list_vertices.size(), {});
    for (Vertex& vertex : m_list_vertices) {
        vertex.set_lies_in_irreducible_wedge(false);
    }

    for (const Vertex& vertex : m_list_vertices) {
        const vector3 representative = fold_positive_octant_to_irreducible_wedge(vertex.get_position());
        const auto    it             = vertex_by_position.find(coordinate_key(representative));
        if (it == vertex_by_position.end()) {
            const vector3 representative_reduced = si_to_reduced_k(representative);
            double        nearest_distance       = std::numeric_limits<double>::infinity();
            vector3       nearest_reduced;
            for (const Vertex& candidate : m_list_vertices) {
                const vector3 candidate_reduced = si_to_reduced_k(candidate.get_position());
                const double  distance          = (candidate_reduced - representative_reduced).norm();
                if (distance < nearest_distance) {
                    nearest_distance = distance;
                    nearest_reduced  = candidate_reduced;
                }
            }
            fmt::print(stderr,
                       "[warn] positive-octant mesh is not closed under coordinate permutations: representative "
                       "({:.12g}, {:.12g}, {:.12g}) is missing (nearest vertex ({:.12g}, {:.12g}, {:.12g}), "
                       "reduced distance {:.3e}). IW reduction is unavailable; generate a symmetry-exact octant "
                       "mesh with bz_meshing --bz-domain octant.\n",
                       representative_reduced.x(),
                       representative_reduced.y(),
                       representative_reduced.z(),
                       nearest_reduced.x(),
                       nearest_reduced.y(),
                       nearest_reduced.z(),
                       nearest_distance);
            m_list_vtx_in_iwedge.clear();
            m_kstar_ibz_to_bz.clear();
            return;
        }
        m_kstar_ibz_to_bz[it->second].push_back(vertex.get_index());
    }

    std::size_t mapped_vertices = 0;
    for (std::size_t index = 0; index < m_kstar_ibz_to_bz.size(); ++index) {
        auto& orbit = m_kstar_ibz_to_bz[index];
        if (orbit.empty()) {
            continue;
        }
        if (!is_irreducible_wedge(m_list_vertices[index].get_position())) {
            throw std::runtime_error("Positive-octant symmetry orbit representative is outside the irreducible wedge");
        }
        std::sort(orbit.begin(), orbit.end());
        orbit.erase(std::unique(orbit.begin(), orbit.end()), orbit.end());
        if (orbit.size() > 6) {
            throw std::runtime_error("Positive-octant symmetry orbit has more than six coordinate permutations");
        }
        m_list_vertices[index].set_lies_in_irreducible_wedge(true);
        m_list_vtx_in_iwedge.push_back(index);
        mapped_vertices += orbit.size();
    }

    if (mapped_vertices != m_list_vertices.size()) {
        throw std::runtime_error("Positive-octant symmetry orbits do not cover every stored vertex");
    }

    fmt::print("Built positive-octant IW mapping: {} representatives cover {} vertices\n",
               m_list_vtx_in_iwedge.size(),
               mapped_vertices);
}

std::size_t MeshBZ::local_vertex_index_from_source(std::size_t source_index) const {
    if (source_index >= m_local_index_from_source_vertex.size()) {
        return invalid_vertex_index;
    }
    return m_local_index_from_source_vertex[source_index];
}

/**
 * @brief Load the mapping from k-points in the irreducible Brillouin zone to the full Brillouin zone.
 *
 * @param filename
 */
void MeshBZ::load_kstar_ibz_to_bz(const std::string& kstarFilePath) {
    m_kstar_ibz_to_bz.clear();
    m_list_vtx_in_iwedge.clear();
    m_kstar_ibz_to_bz.reserve(m_list_vertices.size());

    const std::filesystem::path mesh_path(m_filename_mesh);
    const std::filesystem::path kstar_path     = uepm::utils::find_kstar_file_for_mesh(mesh_path, kstarFilePath);
    const std::string           kstar_filename = kstar_path.string();
    std::cout << "Loading kstar_ibz_to_bz from file : " << kstar_filename << std::endl;
    std::ifstream in(kstar_filename);
    if (!in) {
        throw std::runtime_error("load_kstar_file: can't open " + kstar_filename);
    }

    std::string line;
    std::size_t max_iw = 0;

    struct Row {
        std::size_t              iw;
        std::vector<std::size_t> ids;
    };

    std::vector<Row> tmp;
    std::size_t      count_vtx = 0;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream iss(line);
        std::size_t        iw, m;
        if (!(iss >> iw >> m)) {
            continue;
        }
        std::vector<std::size_t> ids(m);
        for (std::size_t k = 0; k < m; ++k) {
            iss >> ids[k];
        }
        max_iw = std::max(max_iw, iw);
        tmp.push_back({iw, std::move(ids)});
    }
    std::cout << "load_kstar_file: read " << tmp.size() << " rows, max iw = " << max_iw << std::endl;
    std::cout << "Nb of vertices in the mesh: " << m_list_vertices.size() << std::endl;

    m_kstar_ibz_to_bz.assign(max_iw + 1, {});
    for (auto& r : tmp) {
        std::size_t idx_iw        = r.iw;
        m_kstar_ibz_to_bz[idx_iw] = std::move(r.ids);
        const Vertex& vtx         = m_list_vertices[idx_iw];
        if (!is_irreducible_wedge(vtx.get_position())) {
            std::cout << "load_kstar_file: vertex " << idx_iw << " at " << vtx.get_position()
                      << " does not lie in the irreducible wedge, but it should" << std::endl;
            throw std::runtime_error("load_kstar_file: vertex " + std::to_string(idx_iw) +
                                     " does not lie in the irreducible wedge, but it should");
        }
        m_list_vtx_in_iwedge.push_back(idx_iw);
        count_vtx += m_kstar_ibz_to_bz[idx_iw].size();
    }

    if (count_vtx != m_list_vertices.size()) {
        throw std::runtime_error("load_kstar_file: number of k-points in IBZ (" + std::to_string(count_vtx) +
                                 ") does not match number of vertices in mesh (" +
                                 std::to_string(m_list_vertices.size()) + ")");
    }

    for (std::size_t i = 0; i < m_kstar_ibz_to_bz.size(); ++i) {
        for (auto& id : m_kstar_ibz_to_bz[i]) {
            if (id >= m_list_vertices.size()) {
                throw std::runtime_error("load_kstar_file: index " + std::to_string(id) + " out of range");
            }
        }
    }
    for (auto idx : m_list_vtx_in_iwedge) {
        m_list_vertices[idx].set_lies_in_irreducible_wedge(true);
    }
}

bbox_mesh MeshBZ::compute_bounding_box() const {
    if (m_list_vertices.empty()) {
        throw std::logic_error("Cannot compute the bounding box of an empty BZ mesh.");
    }
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
    fmt::print("Building octree search tree for {} tetrahedra ...\n", m_list_tetrahedra.size());
    bbox_mesh mesh_bbox = compute_bounding_box();
    fmt::print("Mesh bounding box: x=[{}, {}], y=[{}, {}], z=[{}, {}]\n",
               mesh_bbox.get_x_min(),
               mesh_bbox.get_x_max(),
               mesh_bbox.get_y_min(),
               mesh_bbox.get_y_max(),
               mesh_bbox.get_z_min(),
               mesh_bbox.get_z_max());
    // const double dilatation_factor = 1.10;
    // mesh_bbox.dilate(dilatation_factor);
    // fmt::print("Dilated x{} mesh bounding box: x=[{}, {}], y=[{}, {}], z=[{}, {}]\n",
    //            dilatation_factor,
    //            mesh_bbox.get_x_min(),
    //            mesh_bbox.get_x_max(),
    //            mesh_bbox.get_y_min(),
    //            mesh_bbox.get_y_max(),
    //            mesh_bbox.get_z_min(),
    //            mesh_bbox.get_z_max());
    auto start    = std::chrono::high_resolution_clock::now();
    m_search_tree = std::make_unique<Octree_mesh>(get_list_p_tetra(), mesh_bbox);
    auto end      = std::chrono::high_resolution_clock::now();
    auto total    = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    fmt::print("Octree built in {}s\n", total / 1000.0);
}

Tetra* MeshBZ::find_tetra_at_location(const vector3& location) const {
    if (!m_search_tree) {
        throw std::logic_error("The BZ search tree has not been built.");
    }
    return m_search_tree->find_tetra_at_location(canonicalize_physical_k(location).representative);
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
void MeshBZ::read_mesh_bands_from_msh_file(const std::string& filename,
                                           int                nb_conduction_bands,
                                           int                nb_valence_bands,
                                           bool               auto_shift_conduction_band,
                                           bool               set_positive_valence_band) {
    m_filename_mesh = filename;
    std::cout << "Opening file " << filename << std::endl;
    GmshSession gmsh_session;
    gmsh::open(filename);
    for (auto& vertex : m_list_vertices) {
        vertex.clear_band_energies();
        vertex.clear_energy_gradient_at_bands();
    }
    m_bands.clear();

    std::unordered_map<std::size_t, std::size_t> local_index_from_node_tag;
    local_index_from_node_tag.reserve(m_node_tags.size());
    for (std::size_t local_index = 0; local_index < m_node_tags.size(); ++local_index) {
        local_index_from_node_tag.emplace(m_node_tags[local_index], local_index);
    }
    // std::cout << "Read gmsh views (band energy values) ..." << std::endl;
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
        std::cout << "Loading band " << count_band << ": " << name_view << std::endl;

        std::string              type;
        std::vector<std::size_t> tags;
        double                   time;
        int                      numComp;
        std::vector<double>      data_view;
        gmsh::view::getHomogeneousModelData(tag, 0, type, tags, data_view, time, numComp);
        if (tags.size() * static_cast<std::size_t>(numComp) != data_view.size()) {
            throw std::runtime_error("Band view '" + name_view + "' has inconsistent node tags and values.");
        }
        std::vector<double> ordered_data(m_list_vertices.size() * static_cast<std::size_t>(numComp));
        std::vector<bool>   node_was_set(m_list_vertices.size(), false);
        for (std::size_t source_index = 0; source_index < tags.size(); ++source_index) {
            const auto local_it = local_index_from_node_tag.find(tags[source_index]);
            if (local_it == local_index_from_node_tag.end()) {
                continue;
            }
            for (int component = 0; component < numComp; ++component) {
                ordered_data[numComp * local_it->second + component] = data_view[numComp * source_index + component];
            }
            node_was_set[local_it->second] = true;
        }
        if (std::find(node_was_set.begin(), node_was_set.end(), false) != node_was_set.end()) {
            throw std::runtime_error("Band view '" + name_view + "' does not provide every mesh node.");
        }

        // Num comp = 1  : band energy values at each vertex
        if (numComp == 1) {
            const auto minmax_band = std::minmax_element(ordered_data.begin(), ordered_data.end());
            const bool is_valence  = *minmax_band.second <= 0.1;
            m_bands.register_band(is_valence ? MeshParticleType::valence : MeshParticleType::conduction,
                                  *minmax_band.first,
                                  *minmax_band.second);
            count_band++;
            add_new_band_energies_to_vertices(ordered_data);
        } else if (numComp == 3) {
            // Band energy gradients at each vertex
            std::cout << "Reading band energy gradients for band " << count_band - 1 << std::endl;
            add_new_gradient_band_energies_to_vertices(ordered_data);
        } else {
            throw std::runtime_error("read_mesh_bands_from_msh_file: unsupported number of components per view: " +
                                     std::to_string(numComp));
        }
    }
    if (m_bands.total() != static_cast<std::size_t>(count_band)) {
        throw std::runtime_error("Band catalog count does not match loaded band views");
    }
    // PRINT INFO
    std::cout << "Number of bands loaded: " << get_number_bands_total() << std::endl;
    print_band_info();
    set_bands_in_right_order();
    if (nb_valence_bands >= 0 || nb_conduction_bands >= 0) {
        keep_only_bands(
            (nb_valence_bands >= 0) ? nb_valence_bands : m_bands.range(MeshParticleType::valence).count,
            (nb_conduction_bands >= 0) ? nb_conduction_bands : m_bands.range(MeshParticleType::conduction).count);
    }

    if (auto_shift_conduction_band) {
        auto_shift_conduction_band_energies();
    }
    if (set_positive_valence_band) {
        auto_set_positive_valence_band_energies();
    }

    compute_min_max_energies_at_tetras();
    compute_energy_gradient_at_tetras();

#pragma omp parallel for schedule(dynamic) num_threads(m_nb_threads_mesh_ops)
    for (auto&& tetra : m_list_tetrahedra) {
        tetra.pre_compute_sorted_slots_per_band();
    }

    // set_energy_gradient_at_vertices_by_averaging_tetras();
    recompute_tetra_ordered_energies(m_max_energy_global);

    print_band_info();
    fmt::print("Done reading band energies from mesh file.\n");
}

void MeshBZ::print_band_info() const {
    fmt::print("\nBand info:\n");
    fmt::print("{:<6} {:<11} {:<10} {:<12} {:<10} {:<10}\n",
               "Index",
               "Type",
               "LocalIdx",
               "GlobalIdx",
               "MinE(eV)",
               "MaxE(eV)");
    for (std::size_t i = 0; i < m_bands.info().size(); ++i) {
        const auto& band     = m_bands.info()[i];
        std::string type_str = (band.type == MeshParticleType::valence) ? "Valence" : "Conduction";
        fmt::print("{:<6} {:<11} {:<10} {:<12} {:<10.4f} {:<10.4f}\n",
                   i,
                   type_str,
                   band.local_index,
                   i,
                   m_bands.minima()[i],
                   m_bands.maxima()[i]);
    }
    fmt::print("\n");
}

/**
 * @brief Apply a scissor (energy window) to the band structure.
 * This operation shifts all conduction band energies by the scissor_value (in eV).
 *
 * @param scissor_value
 */
void MeshBZ::apply_scissor(double scissor_value) {
    std::cout << "Applying scissor of " << scissor_value << " eV to conduction bands ..." << std::endl;
#pragma omp parallel for schedule(dynamic) num_threads(m_nb_threads_mesh_ops)
    for (auto&& vtx : m_list_vertices) {
        const std::size_t nb_bands = vtx.get_number_bands();
        for (std::size_t i = 0; i < nb_bands; ++i) {
            if (m_bands.info()[i].type == MeshParticleType::conduction) {
                double new_energy = vtx.get_energy_at_band(i) + scissor_value;
                vtx.set_band_energy(i, new_energy);
            }
        }
    }
    const bool recompute_min_max = true;
    const bool recompute_grad    = false;
    const bool recompute_dos     = false;
    recompute_energies_data_and_sync(recompute_min_max, recompute_grad, recompute_dos, 0.0, 0.0);
    fmt::print("Done applying scissor.\n");
}

void MeshBZ::precompute_dos_tetra(double energy_step, double energy_max) {
    fmt::print("Precomputing DOS per tetrahedra with energy step = {:.3f} eV ...\n", energy_step);
#pragma omp parallel for schedule(dynamic) num_threads(m_nb_threads_mesh_ops)
    for (std::size_t i = 0; i < m_list_tetrahedra.size(); ++i) {
        m_list_tetrahedra[i].precompute_dos_on_energy_grid_per_band(energy_step, energy_max);
    }
}
void MeshBZ::set_energy_gradient_at_vertices_by_averaging_tetras() {
    fmt::print("Setting energy gradient at vertices by averaging tetrahedra gradients...\n");

    constexpr double  eps = 1e-18;
    const std::size_t nv  = m_list_vertices.size();

    // (Optional but recommended) pre-size per-vertex gradient storage
    for (std::size_t v = 0; v < nv; ++v) {
        const auto nb_bands = m_list_vertices[v].get_number_bands();
        m_list_vertices[v].resize_energy_gradient_at_bands(nb_bands);  // add this API
    }

#pragma omp parallel for schedule(dynamic, 64) num_threads(m_nb_threads_mesh_ops)
    for (std::size_t idx_vtx = 0; idx_vtx < nv; ++idx_vtx) {
        const std::size_t nb_bands = m_list_vertices[idx_vtx].get_number_bands();

        for (std::size_t idx_band = 0; idx_band < nb_bands; ++idx_band) {
            vector3 accum(0.0, 0.0, 0.0);
            double  wsum = 0.0;

            // incident tetrahedra of this vertex (read-only)
            for (auto t_idx : m_vertex_to_tetrahedra[idx_vtx]) {
                const auto&  T  = m_list_tetrahedra[t_idx];
                const double VT = std::abs(T.get_signed_volume());
                if (VT <= eps) {
                    continue;
                }

                const vector3 gT = T.get_gradient_energy_at_band(idx_band);  // P1 constant per tet
                accum += VT * gT;
                wsum += VT;
            }

            const vector3 g_i = (wsum > 0.0) ? (accum / wsum) : vector3(0.0, 0.0, 0.0);
            // Prefer indexed setter to avoid push_back reallocations
            m_list_vertices[idx_vtx].set_energy_gradient_at_band(idx_band, g_i);
        }
    }

    fmt::print("Done setting energy gradient at vertices.\n");
}

void MeshBZ::recompute_min_max_energies() {
    if (m_list_vertices.empty()) {
        throw std::logic_error("Cannot compute band extrema on an empty BZ mesh.");
    }
    m_bands.minima().clear();
    m_bands.maxima().clear();
    int nb_bands = m_list_vertices[0].get_number_bands();
    m_bands.minima().resize(nb_bands, std::numeric_limits<double>::max());
    m_bands.maxima().resize(nb_bands, std::numeric_limits<double>::lowest());
    for (auto&& vtx : m_list_vertices) {
        const auto& energies = vtx.get_band_energies();
        for (int i = 0; i < nb_bands; ++i) {
            if (energies[i] < m_bands.minima()[i]) {
                m_bands.minima()[i] = energies[i];
            }
            if (energies[i] > m_bands.maxima()[i]) {
                m_bands.maxima()[i] = energies[i];
            }
        }
    }
}

void MeshBZ::recompute_energies_data_and_sync(bool   recompute_min_max,
                                              bool   recompute_grad,
                                              bool   recompute_dos,
                                              double dos_energy_step,
                                              double dos_energy_max) {
    compute_min_max_energies_at_tetras();
    if (recompute_min_max) {
        recompute_min_max_energies();
    }
    if (recompute_grad) {
        compute_energy_gradient_at_tetras();
        set_energy_gradient_at_vertices_by_averaging_tetras();
    }
    if (recompute_dos) {
        precompute_dos_tetra(dos_energy_step, dos_energy_max);
    }
    // Check that the band info is still correct
    print_band_info();
    std::size_t nb_band_vtx = m_list_vertices[0].get_number_bands();
    if (nb_band_vtx != m_bands.info().size()) {
        throw std::runtime_error("The number of bands in the vertices does not match the band info. Abort.");
    }
    if (nb_band_vtx != m_bands.minima().size() || nb_band_vtx != m_bands.maxima().size()) {
        throw std::runtime_error("The number of bands in the vertices does not match the min/max band info. Abort.");
    }
    if (nb_band_vtx != get_number_bands_total()) {
        throw std::runtime_error(
            "The number of bands in the vertices does not match the total number of bands. Abort.");
    }
}

void MeshBZ::recompute_tetra_ordered_energies(double max_energy) {
    m_tetra_energy_index.rebuild(m_list_tetrahedra, get_number_bands_total(), max_energy, m_nb_threads_mesh_ops);
}

std::vector<std::size_t> MeshBZ::get_band_indices(MeshParticleType type) const { return m_bands.indices(type); }

/**
 * @brief Keep only a subset of bands.
 * Must be called after reading the band energies from the .msh file and after shifting/setting valence
 * absolute energies and so on.
 * @param required_nb_bands Number of bands to keep (from the top for valence bands, from the bottom for conduction
 * bands).
 */
void MeshBZ::keep_only_bands(std::size_t nb_valence_bands, std::size_t nb_conduction_bands) {
    if (nb_valence_bands + nb_conduction_bands > get_number_bands_total()) {
        fmt::print("Requested to keep {} valence bands and {} conduction bands, which is more than the total number of "
                   "bands ({}). Abort.\n",
                   nb_valence_bands,
                   nb_conduction_bands,
                   get_number_bands_total());
        throw std::runtime_error("Cannot keep more bands than available.");
    }
    // Valence
    std::vector<std::size_t> valence_indices = get_band_indices(MeshParticleType::valence);
    if (nb_valence_bands < valence_indices.size()) {
        int nb_valence_to_remove = valence_indices.size() - nb_valence_bands;
        std::cout << "Removing " << nb_valence_to_remove << " valence bands (keeping " << nb_valence_bands << ")\n";
        for (int i = 0; i < nb_valence_to_remove; ++i) {
            std::size_t band_index_to_remove = valence_indices[valence_indices.size() - 1 - i];
#pragma omp parallel for schedule(dynamic) num_threads(m_nb_threads_mesh_ops)
            for (auto&& vtx : m_list_vertices) {
                vtx.remove_band_energy(band_index_to_remove);
            }
            m_bands.info().erase(m_bands.info().begin() + band_index_to_remove);
            m_bands.range(MeshParticleType::valence).count--;
            m_bands.range(MeshParticleType::conduction).global_start_index--;
        }
    }
    // Conduction
    std::vector<std::size_t> conduction_indices = get_band_indices(MeshParticleType::conduction);
    if (nb_conduction_bands < conduction_indices.size()) {
        int nb_conduction_to_remove = conduction_indices.size() - nb_conduction_bands;
        std::cout << "Removing " << nb_conduction_to_remove << " conduction bands (keeping " << nb_conduction_bands
                  << ")\n";
        for (int i = 0; i < nb_conduction_to_remove; ++i) {
            std::size_t band_index_to_remove = conduction_indices[conduction_indices.size() - 1 - i];
#pragma omp parallel for schedule(dynamic) num_threads(m_nb_threads_mesh_ops)
            for (auto&& vtx : m_list_vertices) {
                vtx.remove_band_energy(band_index_to_remove);
            }
            m_bands.info().erase(m_bands.info().begin() + band_index_to_remove);
            m_bands.range(MeshParticleType::conduction).count--;
        }
    }
    m_bands.set_total(nb_valence_bands + nb_conduction_bands);
    // Recompute min/max band energies
    recompute_min_max_energies();
}

void MeshBZ::add_new_band_energies_to_vertices(const std::vector<double>& energies_at_vertices) {
    if (energies_at_vertices.size() != m_list_vertices.size()) {
        throw std::invalid_argument("The number of energy values does not match the number of vertices. Abort.");
    }
#pragma omp parallel for schedule(static) num_threads(m_nb_threads_mesh_ops)
    for (std::size_t index_vtx = 0; index_vtx < m_list_vertices.size(); ++index_vtx) {
        m_list_vertices[index_vtx].add_band_energy_value(energies_at_vertices[index_vtx]);
    }
}

void MeshBZ::append_band(MeshParticleType            type,
                         const std::vector<double>&  energies_at_vertices,
                         const std::vector<vector3>& gradients_at_vertices) {
    if (energies_at_vertices.size() != m_list_vertices.size() ||
        gradients_at_vertices.size() != m_list_vertices.size()) {
        throw std::invalid_argument("Band energies and gradients must match the number of mesh vertices.");
    }
    if (energies_at_vertices.empty()) {
        throw std::invalid_argument("Cannot append an empty band.");
    }

    const auto minmax = std::minmax_element(energies_at_vertices.begin(), energies_at_vertices.end());
    m_bands.register_band(type, *minmax.first, *minmax.second);

#pragma omp parallel for schedule(static) num_threads(m_nb_threads_mesh_ops)
    for (std::size_t index_vtx = 0; index_vtx < m_list_vertices.size(); ++index_vtx) {
        m_list_vertices[index_vtx].add_band_energy_value(energies_at_vertices[index_vtx]);
        m_list_vertices[index_vtx].add_band_energy_gradient(gradients_at_vertices[index_vtx]);
    }
}

void MeshBZ::add_new_gradient_band_energies_to_vertices(const std::vector<double>& gradients_at_vertices) {
    std::size_t nb_components = 3;  // Gradient has 3 components
    if (gradients_at_vertices.size() != m_list_vertices.size() * nb_components) {
        throw std::invalid_argument("The number of gradient values does not match the number of vertices. Abort.");
    }
#pragma omp parallel for schedule(static) num_threads(m_nb_threads_mesh_ops)
    for (std::size_t index_vtx = 0; index_vtx < m_list_vertices.size(); ++index_vtx) {
        vector3 gradient = {gradients_at_vertices[index_vtx * nb_components + 0],
                            gradients_at_vertices[index_vtx * nb_components + 1],
                            gradients_at_vertices[index_vtx * nb_components + 2]};
        m_list_vertices[index_vtx].add_band_energy_gradient(gradient);
    }
}

void        MeshBZ::compute_min_max_energies_at_tetras() {
#pragma omp parallel for schedule(dynamic) num_threads(m_nb_threads_mesh_ops)
    for (auto&& tetra : m_list_tetrahedra) {
        tetra.compute_min_max_energies_at_bands();
    }
}

void MeshBZ::auto_set_positive_valence_band_energies() {
    std::vector<std::size_t> valence_band_indices = get_band_indices(MeshParticleType::valence);
#pragma omp parallel for schedule(dynamic) num_threads(m_nb_threads_mesh_ops)
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
        if (m_bands.minima()[idx_band] < min_conduction) {
            min_conduction = m_bands.minima()[idx_band];
        }
    }
    std::cout << "Min conduction band energy: " << min_conduction << " eV\n";

    double band_gap = min_conduction - band_reference_eV;
    std::cout << "Computed band gap: " << band_gap << " eV\n";

    if (band_gap < 0.1) {
        std::cout << "Warning: computed band gap is very small or negative (" << band_gap
                  << " eV). No shift applied.\n";
        return;
    }

    // Shift all conduction bands so that the conduction band minimum aligns with the valence band maximum
    double shift_amount = -band_gap;
    std::cout << "Shifting conduction bands by " << shift_amount << " eV to align CBM with VBM.\n";
#pragma omp parallel for schedule(dynamic) num_threads(m_nb_threads_mesh_ops)
    for (auto&& vtx : m_list_vertices) {
        for (auto idx_band : conduction_band_indices) {
            double old_energy = vtx.get_energy_at_band(idx_band);
            double new_energy = old_energy + shift_amount;
            vtx.set_band_energy(idx_band, new_energy);
        }
    }

    // Recompute per-band min/max
    m_bands.minima().clear();
    m_bands.maxima().clear();
    for (std::size_t b = 0; b < m_list_vertices[0].get_number_bands(); ++b) {
        double bmin = std::numeric_limits<double>::max();
        double bmax = std::numeric_limits<double>::lowest();
        for (auto&& vtx : m_list_vertices) {
            double e = vtx.get_energy_at_band(b);
            bmin     = std::min(bmin, e);
            bmax     = std::max(bmax, e);
        }
        m_bands.minima().push_back(bmin);
        m_bands.maxima().push_back(bmax);
    }

    std::cout << "Post-shift band extrema:\n";
    for (std::size_t b = 0; b < m_bands.minima().size(); ++b) {
        std::cout << " Band " << b << ": min = " << m_bands.minima()[b] << " eV, max = " << m_bands.maxima()[b]
                  << " eV\n";
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
#pragma omp parallel for schedule(dynamic) num_threads(m_nb_threads_mesh_ops)
    for (auto&& vtx : m_list_vertices) {
        vtx.reverse_energies_order(global_idx_first, global_idx_last);
    }

    recompute_min_max_energies();
}

void        MeshBZ::compute_energy_gradient_at_tetras() {
#pragma omp parallel for schedule(dynamic) num_threads(m_nb_threads_mesh_ops)
    for (auto&& tetra : m_list_tetrahedra) {
        tetra.compute_gradient_energy_at_bands();
    }
}

vector3 MeshBZ::interpolate_energy_gradient_at_location(const vector3& location, const std::size_t& idx_band) const {
    const CanonicalK canonical = canonicalize_physical_k(location);
    Tetra*           tetra     = find_tetra_at_location(canonical.representative);
    if (tetra == nullptr) {
        throw std::runtime_error("Location is outside the mesh. Cannot interpolate energy gradient.");
    }
    return representative_vector_to_physical(
        tetra->interpolate_gradient_energy_at_band(canonical.representative, idx_band),
        canonical.signs);
}

double MeshBZ::compute_mesh_volume() const {
    double total_volume = 0.0;
#pragma omp parallel for schedule(dynamic) num_threads(m_nb_threads_mesh_ops) reduction(+ : total_volume)
    for (auto&& tetra : m_list_tetrahedra) {
        total_volume += std::fabs(tetra.get_signed_volume());
    }
    return total_volume;
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
    double      min_dist = std::numeric_limits<double>::max();
    std::size_t idx_min  = 0;
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

void MeshBZ::compute_band_structure_over_mesh(uepm::pseudopotential::BandStructure& band_structure, bool use_iwedge) {
    if (stores_positive_octant() && use_iwedge) {
        throw std::invalid_argument(
            "Cannot combine positive-octant storage with irreducible-wedge band computation");
    }
    const auto& full_list_vertices            = m_list_vertices;
    const auto& list_idx_irreducible_vertices = m_list_vtx_in_iwedge;

    const bool               compute_gradient = true;
    std::vector<std::size_t> list_vtx_used    = use_iwedge ? list_idx_irreducible_vertices : std::vector<std::size_t>{};
    const auto&              nb_vtx_used      = use_iwedge ? m_list_vtx_in_iwedge.size() : m_list_vertices.size();
    if (use_iwedge) {
        fmt::print("Computing band structure over the irreducible wedge with {} vertices...\n", nb_vtx_used);
    } else {
        fmt::print("Computing band structure over the full BZ with {} vertices...\n", nb_vtx_used);
        list_vtx_used = std::vector<std::size_t>(m_list_vertices.size());
        std::iota(list_vtx_used.begin(), list_vtx_used.end(), 0);
    }

    std::vector<Vector3D<double>> mesh_kpoints(nb_vtx_used);
    const double                  si_to_red = si_to_reduced_scale();
    std::cout << "Total number of vertices in the BZ mesh: " << m_list_vertices.size() << std::endl;
#pragma omp parallel for schedule(dynamic) num_threads(m_nb_threads_mesh_ops)
    for (std::size_t i = 0; i < nb_vtx_used; ++i) {
        const auto& vtx = full_list_vertices[list_vtx_used[i]];
        mesh_kpoints[i] = Vector3D<double>(vtx.get_position().x() * si_to_red,
                                           vtx.get_position().y() * si_to_red,
                                           vtx.get_position().z() * si_to_red);
    }
    std::cout << "Number of k-points in the irreducible BZ: " << mesh_kpoints.size() << std::endl;
    band_structure.set_kpoints(mesh_kpoints);
    band_structure.Compute_parallel(compute_gradient, m_nb_threads_mesh_ops);
    bool set_cond_band_zero = false;
    band_structure.AdjustValues(set_cond_band_zero);

    // Now set the computed energies to the mesh vertices

#pragma omp parallel for schedule(dynamic) num_threads(m_nb_threads_mesh_ops)
    for (std::size_t i = 0; i < nb_vtx_used; ++i) {
        const auto& energies_at_vtx  = band_structure.get_band_energies().at(i);
        const auto& gradients_at_vtx = band_structure.get_band_energy_gradients().at(i);
        // Convert Vector3D<double> to vector3
        std::vector<vector3> grad_vectors;
        for (const auto& grad : gradients_at_vtx) {
            grad_vectors.emplace_back(grad.X, grad.Y, grad.Z);
        }
        m_list_vertices[list_vtx_used[i]].set_all_band_energies(energies_at_vtx);
        m_list_vertices[list_vtx_used[i]].set_energy_gradient_at_bands(grad_vectors);
    }

    fmt::print("Band structure computed over with {} bands.\n", band_structure.get_number_of_bands());
    if (use_iwedge) {
        fmt::print("Distributing energies from irreducible wedge to full BZ...\n");
        distribute_energies_from_iw_wedge_to_full_bz();
    }
    fmt::print("Recomputing min/max energies...\n");
    recompute_min_max_energies();
    // Set band info
    m_bands.clear_classification();
    constexpr double eps              = 1e-6;
    const auto       minimum_energies = m_bands.minima();
    const auto       maximum_energies = m_bands.maxima();
    for (std::size_t band_index = 0; band_index < minimum_energies.size(); ++band_index) {
        const bool is_valence = minimum_energies[band_index] < eps;
        m_bands.classify_existing_band(is_valence ? MeshParticleType::valence : MeshParticleType::conduction);
    }
    print_band_info();
}

/**
 * @brief Apply the same symmetry operation to a vector in the BZ.
 * Un poco brute-force en attendant de faire mieux...
 *
 * @param k_iw
 * @param k_bz
 * @param v
 * @return vector3
 */
inline vector3 apply_same_symmetry_operation(const vector3& k_iw, const vector3& k_bz, const vector3& v) {
    constexpr double eps = 1e-8;

    std::array<double, 3> a = {k_iw.x(), k_iw.y(), k_iw.z()};
    std::array<double, 3> b = {k_bz.x(), k_bz.y(), k_bz.z()};

    std::array<int, 3>  perm   = {-1, -1, -1};
    std::array<int, 3>  sign   = {1, 1, 1};
    std::array<bool, 3> used_a = {false, false, false};

    double scale = 0.0;
    for (double x : a) {
        scale = std::max(scale, std::abs(x));
    }
    for (double x : b) {
        scale = std::max(scale, std::abs(x));
    }
    if (scale < 1.0) {
        scale = 1.0;
    }
    const double tol = eps * scale;

    // For each component of k_bz, find the matching component of k_iw up to |·|
    for (int jb = 0; jb < 3; ++jb) {
        bool found = false;
        for (int ia = 0; ia < 3; ++ia) {
            if (used_a[ia]) {
                continue;
            }
            if (std::abs(std::abs(a[ia]) - std::abs(b[jb])) <= tol) {
                perm[jb]   = ia;
                used_a[ia] = true;

                if (std::abs(a[ia]) > tol) {
                    const int sa = (a[ia] >= 0.0) ? 1 : -1;
                    const int sb = (b[jb] >= 0.0) ? 1 : -1;
                    sign[jb]     = sb * sa;
                } else {
                    // a[ia] ~ 0: sign of the mapping does not really matter; follow k_bz
                    sign[jb] = (b[jb] >= 0.0) ? 1 : -1;
                }
                found = true;
                break;
            }
        }
        if (!found) {
            // Could not represent the mapping as pure signed permutation; do nothing.
            return v;
        }
    }

    auto comp = [&](int idx) -> double {
        switch (idx) {
            case 0:
                return v.x();
            case 1:
                return v.y();
            default:
                return v.z();
        }
    };

    const double vx = sign[0] * comp(perm[0]);
    const double vy = sign[1] * comp(perm[1]);
    const double vz = sign[2] * comp(perm[2]);

    return vector3{vx, vy, vz};
}

void MeshBZ::distribute_energies_from_iw_wedge_to_full_bz() {
    std::cout << "Distributing band energies from irreducible wedge to full BZ ..." << std::endl;
    if (m_list_vtx_in_iwedge.empty()) {
        throw std::runtime_error("No vertices in irreducible wedge. Cannot distribute energies to full BZ.");
    }
    const auto& ref_vtx      = m_list_vertices[m_list_vtx_in_iwedge[0]];
    const auto& ref_energies = ref_vtx.get_band_energies();
    std::size_t nb_bands     = ref_energies.size();
    if (nb_bands == 0) {
        throw std::runtime_error(
            "No band energies found in the irreducible wedge vertex. Cannot distribute energies to full BZ.");
    }
#pragma omp parallel for schedule(dynamic) num_threads(m_nb_threads_mesh_ops)
    for (const auto& idx_vtx : m_list_vtx_in_iwedge) {
        const auto&                 vtx_iw           = m_list_vertices[idx_vtx];
        const auto&                 list_energies    = vtx_iw.get_band_energies();
        const std::vector<vector3>& energy_gradients = vtx_iw.get_energy_gradient_at_bands();
        if (list_energies.size() != nb_bands || energy_gradients.size() != nb_bands) {
            throw std::runtime_error(
                "Inconsistent number of bands in IW vertex. Cannot distribute energies to full BZ.");
        }
        // Get the 48 symmetry-equivalent k-points
        std::vector<std::size_t> sym_eq_indices = m_kstar_ibz_to_bz[idx_vtx];
        for (const auto& idx_sym : sym_eq_indices) {
            if (idx_sym == idx_vtx) {
                continue;  // already set
            }
            auto& vtx_full = m_list_vertices[idx_sym];
            vtx_full.set_all_band_energies(list_energies);
            // Transform gradients using the same symmetry operation
            std::vector<vector3> transformed_gradients;
            for (std::size_t i = 0; i < energy_gradients.size(); ++i) {
                transformed_gradients.push_back(
                    apply_same_symmetry_operation(vtx_iw.get_position(), vtx_full.get_position(), energy_gradients[i]));
            }
            vtx_full.set_energy_gradient_at_bands(transformed_gradients);
            // DEBUG
            // vector3              null_grad{0.0, 0.0, 0.0};
            // std::vector<vector3> null_grads(nb_bands, null_grad);
            // vtx_full.set_energy_gradient_at_bands(null_grads);
        }
    }
}

void MeshBZ::export_selected_bands_to_gmsh(const std::string& out_filename,
                                           std::size_t        nb_valence_to_export,
                                           std::size_t        nb_conduction_to_export,
                                           bool               highest_valence_as_band0,
                                           const std::string& model_name_or_msh_path,
                                           bool               write_gradients) const {
    const std::size_t nv = m_list_vertices.size();
    if (nv == 0 || m_node_tags.size() != nv) {
        throw std::runtime_error("Mesh vertices/tags not initialized or inconsistent.");
    }

    auto gather_band = [&](int global_band_index) -> std::vector<double> {
        std::vector<double> v(nv);
        for (std::size_t i = 0; i < nv; ++i) {
            const auto& e = m_list_vertices[i].get_band_energies();
            if (static_cast<int>(e.size()) <= global_band_index) {
                throw std::runtime_error("Vertex band vector too small for requested band index.");
            }
            v[i] = e[static_cast<std::size_t>(global_band_index)];
        }
        return v;
    };

    auto gather_grad_vector = [&](int global_band_index) -> std::vector<double> {
        std::vector<double> v;
        v.reserve(3 * nv);
        for (std::size_t i = 0; i < nv; ++i) {
            const auto& grads = m_list_vertices[i].get_energy_gradient_at_bands();
            if (static_cast<int>(grads.size()) <= global_band_index) {
                throw std::runtime_error("Vertex gradient vector too small for requested band index.");
            }
            const auto& g = grads[static_cast<std::size_t>(global_band_index)];
            v.push_back(g.x());
            v.push_back(g.y());
            v.push_back(g.z());
        }
        return v;
    };

    const int   verbose_level = 1;
    GmshSession guard(verbose_level);

    if (!model_name_or_msh_path.empty() && model_name_or_msh_path.size() > 4 &&
        (model_name_or_msh_path.ends_with(".msh") || model_name_or_msh_path.ends_with(".msh2") ||
         model_name_or_msh_path.ends_with(".msh4"))) {
        gmsh::open(model_name_or_msh_path);
    } else {
        gmsh::model::add(model_name_or_msh_path.empty() ? "bz_mesh" : model_name_or_msh_path);
    }

    std::string model_file_name;
    gmsh::model::getCurrent(model_file_name);
    gmsh::option::setNumber("Mesh.Binary", 1);

    bool write_mesh = true;
    if (std::filesystem::exists(out_filename)) {
        write_mesh = false;
    }
    std::size_t out_idx = 0;

    auto write_one_view = [&](const std::string& name, const std::vector<double>& vals) {
        const int data_tag = gmsh::view::add(name);
        gmsh::view::addHomogeneousModelData(data_tag, /*step=*/0, model_file_name, "NodeData", m_node_tags, vals);

        const int index_view = gmsh::view::getIndex(data_tag);
        gmsh::option::setNumber("View[" + std::to_string(index_view) + "].Visible", 0);
        gmsh::option::setNumber("PostProcessing.SaveMesh", write_mesh ? 1 : 0);
        gmsh::view::write(data_tag, out_filename, true);
        write_mesh = false;
    };

    auto write_one_vector_view = [&](const std::string& name, const std::vector<double>& vecVals /* size = 3*nv */) {
        const int data_tag = gmsh::view::add(name);
        // numComponents = 3 tells Gmsh this is a vector at each node
        gmsh::view::addHomogeneousModelData(data_tag,
                                            /*step=*/0,
                                            model_file_name,
                                            "NodeData",
                                            m_node_tags,
                                            vecVals,
                                            /*numComponents=*/3);

        const int index_view = gmsh::view::getIndex(data_tag);
        gmsh::option::setNumber("View[" + std::to_string(index_view) + "].Visible", 0);
        gmsh::option::setNumber("PostProcessing.SaveMesh", write_mesh ? 1 : 0);
        gmsh::view::write(data_tag, out_filename, true);
        write_mesh = false;
    };

    // Valence
    if (nb_valence_to_export > 0) {
        const int v_start = m_bands.range(MeshParticleType::valence).global_start_index;
        const int v_last  = v_start + static_cast<int>(m_bands.range(MeshParticleType::valence).count) - 1;

        if (highest_valence_as_band0) {
            for (int g = v_last; g >= v_last - static_cast<int>(nb_valence_to_export) + 1; --g) {
                write_one_view("band_" + std::to_string(out_idx++), gather_band(g));
                if (write_gradients) {
                    write_one_vector_view("grad_band_" + std::to_string(out_idx - 1), gather_grad_vector(g));
                }
            }
        } else {
            for (int g = v_start; g < v_start + static_cast<int>(nb_valence_to_export); ++g) {
                write_one_view("band_" + std::to_string(out_idx++), gather_band(g));
                if (write_gradients) {
                    write_one_vector_view("grad_band_" + std::to_string(out_idx - 1), gather_grad_vector(g));
                }
            }
        }
    } else {
        std::cout << "No valence bands requested for export." << std::endl;
    }

    // Conduction
    if (nb_conduction_to_export > 0) {
        for (std::size_t g = 0; g < nb_conduction_to_export; ++g) {
            int global_g = get_global_band_index(static_cast<int>(g), MeshParticleType::conduction);
            write_one_view("band_" + std::to_string(out_idx++), gather_band(global_g));
            if (write_gradients) {
                write_one_vector_view("grad_band_" + std::to_string(out_idx - 1), gather_grad_vector(global_g));
            }
        }
    } else {
        std::cout << "No conduction bands requested for export." << std::endl;
    }
}

}  // namespace uepm::mesh_bz
