/**
 * @file dielectric_mesh.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-05-17
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "dielectric_mesh.hpp"

#include <map>
#include <optional>
#include <stdexcept>
#include <unordered_map>

#include "gmsh.h"
#include "gmsh_guard.hpp"

namespace uepm::mesh_bz {

void DielectricMesh::read_dielectric_file(const std::string& filename) {
    std::cout << "Opening file " << filename << std::endl;
    GmshSession gmsh_session;
    gmsh::open(filename);
    std::vector<int> viewTags;
    gmsh::view::getTags(viewTags);
    if (viewTags.empty()) {
        throw std::runtime_error("The dielectric file contains no Gmsh views.");
    }

    std::unordered_map<std::size_t, std::size_t> local_index_from_node_tag;
    local_index_from_node_tag.reserve(m_node_tags.size());
    for (std::size_t local_index = 0; local_index < m_node_tags.size(); ++local_index) {
        local_index_from_node_tag.emplace(m_node_tags[local_index], local_index);
    }

    struct DielectricAtEnergy {
        std::optional<std::vector<double>> real;
        std::optional<std::vector<double>> imag;
    };
    std::map<double, DielectricAtEnergy> data_by_energy;

    for (auto&& tag : viewTags) {
        const int   index_view  = gmsh::view::getIndex(tag);
        std::string name_object = "View[" + std::to_string(index_view) + "].Name";
        std::string name_view;
        try {
            gmsh::option::getString(name_object, name_view);
        } catch (const std::exception& e) {
            std::cerr << e.what() << '\n';
        }

        const std::size_t separator = name_view.find_last_of('_');
        if (separator == std::string::npos || separator + 1 >= name_view.size()) {
            throw std::runtime_error("Cannot extract dielectric energy from view name '" + name_view + "'.");
        }
        const double energy = std::stod(name_view.substr(separator + 1));

        std::string              type;
        std::vector<std::size_t> tags;
        double                   time;
        int                      numComp;
        std::vector<double>      data_view;
        gmsh::view::getHomogeneousModelData(tag, 0, type, tags, data_view, time, numComp);
        if (numComp != 1 || tags.size() != data_view.size()) {
            throw std::runtime_error("Dielectric view '" + name_view + "' must contain one scalar value per tagged node.");
        }

        std::vector<double> values_by_local_node(m_list_vertices.size());
        std::vector<bool>   node_was_set(m_list_vertices.size(), false);
        for (std::size_t value_index = 0; value_index < tags.size(); ++value_index) {
            const auto local_it = local_index_from_node_tag.find(tags[value_index]);
            if (local_it == local_index_from_node_tag.end()) {
                throw std::runtime_error("Dielectric view references unknown node tag " + std::to_string(tags[value_index]) + ".");
            }
            values_by_local_node[local_it->second] = data_view[value_index];
            node_was_set[local_it->second]         = true;
        }
        if (std::find(node_was_set.begin(), node_was_set.end(), false) != node_was_set.end()) {
            throw std::runtime_error("Dielectric view '" + name_view + "' does not provide every mesh node.");
        }

        auto& values_at_energy = data_by_energy[energy];
        if (name_view.starts_with("eps_r")) {
            if (values_at_energy.real) {
                throw std::runtime_error("Duplicate real dielectric view at energy " + std::to_string(energy) + " eV.");
            }
            values_at_energy.real = std::move(values_by_local_node);
        } else if (name_view.starts_with("eps_i")) {
            if (values_at_energy.imag) {
                throw std::runtime_error("Duplicate imaginary dielectric view at energy " + std::to_string(energy) + " eV.");
            }
            values_at_energy.imag = std::move(values_by_local_node);
        } else {
            throw std::runtime_error("Unrecognized dielectric view name '" + name_view + "'.");
        }
    }

    m_energies.clear();
    m_dielectric_function.assign(m_list_vertices.size(), {});
    for (const auto& [energy, values] : data_by_energy) {
        if (!values.real || !values.imag) {
            throw std::runtime_error("Missing real or imaginary dielectric view at energy " + std::to_string(energy) + " eV.");
        }
        m_energies.push_back(energy);
    }
    for (std::size_t idx_node = 0; idx_node < m_list_vertices.size(); ++idx_node) {
        m_dielectric_function[idx_node].reserve(m_energies.size());
        for (const auto& [energy, values] : data_by_energy) {
            m_dielectric_function[idx_node].emplace_back((*values.real)[idx_node], (*values.imag)[idx_node]);
        }
    }
    std::cout << "Size of the dielectric function: " << m_dielectric_function.size() << " x " << m_dielectric_function[0].size()
              << std::endl;
    std::cout << "Dielectric function read." << std::endl;
}

std::pair<std::size_t, double> DielectricMesh::find_closest_energy(double energy) const {
    if (m_energies.empty()) {
        throw std::logic_error("Cannot interpolate an empty dielectric-energy table.");
    }
    if (m_energies.size() == 1) {
        return std::make_pair(0, 0.0);
    }
    if (energy <= m_energies.front()) {
        return std::make_pair(0, 0.0);
    }
    if (energy >= m_energies.back()) {
        return std::make_pair(m_energies.size() - 2, 1.0);
    }
    auto it = std::lower_bound(m_energies.begin(), m_energies.end(), energy);
    if (it == m_energies.begin()) {
        return std::make_pair(0, 0.0);
    }
    if (it == m_energies.end()) {
        return std::make_pair(m_energies.size() - 1, 0.0);
    }
    std::size_t idx = std::distance(m_energies.begin(), it);
    double      t   = (energy - m_energies[idx - 1]) / (m_energies[idx] - m_energies[idx - 1]);
    return std::make_pair(idx - 1, t);
}

complex_d DielectricMesh::interpolate_dielectric_function(const vector3& k, double energy) const {
    if (m_dielectric_function.size() != m_list_vertices.size() || m_energies.empty()) {
        throw std::logic_error("Dielectric data have not been loaded for this mesh.");
    }
    // The dielectric mesh is stored in the positive symmetry octant.
    auto   k_positive{vector3{std::abs(k.x()), std::abs(k.y()), std::abs(k.z())}};
    Tetra* p_tetra = find_tetra_at_location(k_positive);
    if (p_tetra == nullptr) {
        throw std::out_of_range("No dielectric tetrahedron contains k = " + std::to_string(k_positive.x()) + "," +
                                std::to_string(k_positive.y()) + "," + std::to_string(k_positive.z()) + ".");
    }
    std::array<double, 4> barycentric_coordinates = p_tetra->compute_barycentric_coordinates(k_positive);

    std::array<std::size_t, 4>        list_indices_vertices = p_tetra->get_list_indices_vertices();
    std::pair<std::size_t, double>    closest_energy        = find_closest_energy(energy);
    std::size_t                       idx_energy            = closest_energy.first;
    double                            t                     = closest_energy.second;
    std::vector<std::complex<double>> dielectric_function_low(4);
    std::vector<std::complex<double>> dielectric_function_high(4);
    for (std::size_t idx_vertex = 0; idx_vertex < 4; ++idx_vertex) {
        const std::size_t vertex_index = list_indices_vertices[idx_vertex];
        dielectric_function_low[idx_vertex] = m_dielectric_function[vertex_index][idx_energy];
        dielectric_function_high[idx_vertex] =
            m_dielectric_function[vertex_index][std::min(idx_energy + 1, m_energies.size() - 1)];
    }
    std::complex<double> dielectric_function_interpolated_low =
        p_tetra->interpolate_at_position(barycentric_coordinates, dielectric_function_low);
    std::complex<double> dielectric_function_interpolated_high =
        p_tetra->interpolate_at_position(barycentric_coordinates, dielectric_function_high);

    return (1 - t) * dielectric_function_interpolated_low + t * dielectric_function_interpolated_high;
}

}  // namespace uepm::mesh_bz
