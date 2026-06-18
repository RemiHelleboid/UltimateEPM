/**
 * @file bz_mesh_export.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-06-18
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#include <fstream>
#include <iomanip>
#include <stdexcept>

#include "bz_mesh.hpp"
#include "export_octree_vtu.hpp"

namespace uepm::mesh_bz {

namespace {

void write_vtk_scalars(std::ofstream&             output,
                       const std::string&         name,
                       const std::vector<double>& values,
                       std::size_t                expected_count) {
    if (values.size() != expected_count) {
        throw std::runtime_error("VTK scalar field '" + name + "' has an invalid element count");
    }
    output << "SCALARS " << name << " double 1\n";
    output << "LOOKUP_TABLE default\n";
    output << std::setprecision(8);
    for (double value : values) {
        output << value << "\n";
    }
}

void write_vtk_vectors(std::ofstream&              output,
                       const std::string&          name,
                       const std::vector<vector3>& values,
                       std::size_t                 expected_count) {
    if (values.size() != expected_count) {
        throw std::runtime_error("VTK vector field '" + name + "' has an invalid element count");
    }
    output << "VECTORS " << name << " double\n";
    output << std::setprecision(8);
    for (const vector3& value : values) {
        output << value.x() << " " << value.y() << " " << value.z() << "\n";
    }
}

}  // namespace

void MeshBZ::export_k_points_to_file(const std::string& filename) const {
    std::ofstream output(filename);
    if (!output) {
        throw std::invalid_argument("Could not open " + filename + " for writing");
    }
    for (const Vertex& vertex : m_list_vertices) {
        const vector3& position = vertex.get_position();
        output << position.x() << "," << position.y() << "," << position.z() << "\n";
    }
}

void MeshBZ::export_to_vtk(const std::string&        filename,
                           const MapStringToDoubles& point_scalars,
                           const MapStringToVectors& point_vectors,
                           const MapStringToDoubles& cell_scalars,
                           const MapStringToVectors& cell_vectors) const {
    const std::size_t number_of_points = m_list_vertices.size();
    const std::size_t number_of_cells  = m_list_tetrahedra.size();

    std::ofstream output(filename);
    if (!output) {
        throw std::runtime_error("Cannot open VTK file '" + filename + "' for writing");
    }

    output << "# vtk DataFile Version 4.2\n"
              "BZ mesh export\n"
              "ASCII\n"
              "DATASET UNSTRUCTURED_GRID\n";

    output << "POINTS " << number_of_points << " double\n";
    output << std::setprecision(8);
    for (const Vertex& vertex : m_list_vertices) {
        const vector3& position = vertex.get_position();
        output << position.x() << " " << position.y() << " " << position.z() << "\n";
    }

    output << "CELLS " << number_of_cells << " " << number_of_cells * 5 << "\n";
    for (const Tetra& tetra : m_list_tetrahedra) {
        const auto indices = tetra.get_list_indices_vertices();
        output << "4 " << indices[0] << " " << indices[1] << " " << indices[2] << " " << indices[3] << "\n";
    }

    output << "CELL_TYPES " << number_of_cells << "\n";
    for (std::size_t index = 0; index < number_of_cells; ++index) {
        output << "10\n";
    }

    if (!point_scalars.empty() || !point_vectors.empty()) {
        output << "POINT_DATA " << number_of_points << "\n";
        for (const auto& [name, values] : point_scalars) {
            write_vtk_scalars(output, name, values, number_of_points);
        }
        for (const auto& [name, values] : point_vectors) {
            write_vtk_vectors(output, name, values, number_of_points);
        }
    }

    if (!cell_scalars.empty() || !cell_vectors.empty()) {
        output << "CELL_DATA " << number_of_cells << "\n";
        for (const auto& [name, values] : cell_scalars) {
            write_vtk_scalars(output, name, values, number_of_cells);
        }
        for (const auto& [name, values] : cell_vectors) {
            write_vtk_vectors(output, name, values, number_of_cells);
        }
    }
}

void MeshBZ::export_energies_and_gradients_to_vtk(const std::string& filename) const {
    MapStringToDoubles point_scalars;
    MapStringToVectors point_vectors;

    if (m_list_vertices.empty()) {
        throw std::runtime_error("Cannot export band data from an empty mesh");
    }

    for (std::size_t band = 0; band < m_list_vertices.front().get_number_bands(); ++band) {
        auto& values = point_scalars["band_energy_" + std::to_string(band)];
        values.reserve(m_list_vertices.size());
        for (const Vertex& vertex : m_list_vertices) {
            values.push_back(vertex.get_energy_at_band(band));
        }
    }

    for (std::size_t band = 0; band < m_list_vertices.front().get_energy_gradient_at_bands().size(); ++band) {
        auto& values = point_vectors["band_grad_" + std::to_string(band)];
        values.reserve(m_list_vertices.size());
        for (const Vertex& vertex : m_list_vertices) {
            values.push_back(vertex.get_energy_gradient_at_band(band));
        }
    }

    export_to_vtk(filename, point_scalars, point_vectors, {}, {});
}

void MeshBZ::export_octree_to_vtu(const std::string& filename) const {
    if (!m_search_tree) {
        throw std::runtime_error("Octree not initialized; cannot export it");
    }
    write_octree_as_vtu(*m_search_tree, filename, true);
}

}  // namespace uepm::mesh_bz
