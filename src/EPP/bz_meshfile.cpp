/**
 * @file bz_mesh.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-07-07
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "bz_meshfile.hpp"

#include "BandStructure.h"
#include "gmsh.h"

void bz_mesh_points::add_k_point(Vector3D<double> kpoint) { m_kpoints.push_back(kpoint); }

void bz_mesh_points::add_k_point(double k_x, double k_y, double k_z) { m_kpoints.push_back(Vector3D<double>(k_x, k_y, k_z)); }

void bz_mesh_points::read_mesh() {
    std::cout << "Opening file " << m_filename << std::endl;
    gmsh::initialize();
    // gmsh::option::setNumber("General.Verbosity", 1);
    gmsh::model::add("bz_mesh");
    gmsh::open(m_filename);
    // std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords;
    std::vector<double> nodeParams;
    gmsh::model::mesh::reclassifyNodes();
    gmsh::model::mesh::getNodes(m_node_tags, nodeCoords, nodeParams, -1, -1, false, false);
    std::size_t size_nodes_tags        = m_node_tags.size();
    std::size_t size_nodes_coordinates = nodeCoords.size();
    if (size_nodes_coordinates != 3 * size_nodes_tags) {
        throw std::runtime_error("Number of coordinates is not 3 times the number of vertices. Abort.");
    }
    for (std::size_t index_vertex = 0; index_vertex < size_nodes_tags; ++index_vertex) {
        std::size_t node_tag = m_node_tags[index_vertex];
        add_k_point(nodeCoords[3 * index_vertex], nodeCoords[3 * index_vertex + 1], nodeCoords[3 * index_vertex + 2]);
    }
    gmsh::finalize();
    std::cout << "Number of k-points: " << m_kpoints.size() << std::endl;
}

void bz_mesh_points::add_band_on_mesh(const std::string& band_name, const std::vector<double>& band_values) {
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 99999999);
    gmsh::model::add("bz_mesh");
    gmsh::open(m_filename);

    std::string model_file_name;
    gmsh::model::getCurrent(model_file_name);

    int data_tag = gmsh::view::add(band_name);
    if (m_node_tags.size() != band_values.size()) {
        std::cout << "number of nodes: " << m_node_tags.size() << std::endl;
        std::cout << "number of values: " << band_values.size() << std::endl;
        throw std::runtime_error("Number of nodes and number of values are not the same. Abort.");
    }
    gmsh::view::addHomogeneousModelData(data_tag, 0, model_file_name, "NodeData", m_node_tags, band_values);
    const int   index_view             = gmsh::view::getIndex(data_tag);
    std::string name_object_visibility = "View[" + std::to_string(index_view) + "].Visible";
    gmsh::option::setNumber(name_object_visibility, 1);
    gmsh::view::write(data_tag, "band_mesh.msh", true);
    gmsh::finalize();
}

void bz_mesh_points::add_all_bands_on_mesh(const std::string& out_filename, const EmpiricalPseudopotential::BandStructure& my_band) {
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 0);
    gmsh::model::add("bz_mesh");
    gmsh::open(m_filename);

    std::string model_file_name;
    gmsh::model::getCurrent(model_file_name);

    for (int index_band = 0; index_band < my_band.get_number_of_bands(); ++index_band) {
        std::string         band_name   = "band_" + std::to_string(index_band);
        std::vector<double> band_values = my_band.get_band(index_band);
        int                 data_tag    = gmsh::view::add(band_name);
        if (m_node_tags.size() != band_values.size()) {
            std::cout << "number of nodes: " << m_node_tags.size() << std::endl;
            std::cout << "number of values: " << band_values.size() << std::endl;
            throw std::runtime_error("Number of nodes and number of values are not the same. Abort.");
        }
        gmsh::view::addHomogeneousModelData(data_tag, 0, model_file_name, "NodeData", m_node_tags, band_values);
        const int   index_view             = gmsh::view::getIndex(data_tag);
        std::string name_object_visibility = "View[" + std::to_string(index_view) + "].Visible";
        gmsh::option::setNumber(name_object_visibility, 0);
        gmsh::view::write(data_tag, out_filename, true);
    }
    gmsh::finalize();
}

/**
 * @brief Add band structure energies to gmsh mesh as views.
 * One view per band.
 *
 * The band structure is assumed to be given in band_values vector under the following format:
 * band_values[index_k_point * number_of_bands + index_band] = energy of the band with index index_band at k-point with index index_k_point.
 *
 * @param out_filename
 * @param band_values
 */
void bz_mesh_points::add_all_bands_on_mesh(const std::string& out_filename, const std::vector<double>& band_values, int number_bands) {
    if (band_values.size() != number_bands * m_node_tags.size()) {
        throw std::runtime_error("band_values vector is not the same size as the number of bands times the number of nodes. Abort.");
    }
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 0);
    gmsh::model::add("bz_mesh");
    gmsh::open(m_filename);

    std::string model_file_name;
    gmsh::model::getCurrent(model_file_name);

    // for (std::size_t index_value = 0; index_value < band_values.size(); ++index_value) {
    //     std::cout << "index_value: " << index_value << "  --->  " << band_values[index_value] << std::endl;
    // }

    for (int index_band = 0; index_band < number_bands; ++index_band) {
        std::string         band_name = "band_" + std::to_string(index_band);
        std::vector<double> current_band_values(m_node_tags.size());
        for (std::size_t index_node = 0; index_node < m_node_tags.size(); ++index_node) {
            current_band_values[index_node] = band_values[index_node * number_bands + index_band];
            // std::cout << "band_values[" << index_node << "]: " << band_values[index_node * number_bands + index_band] << std::endl;
        }
        int data_tag = gmsh::view::add(band_name);
        if (m_node_tags.size() != current_band_values.size()) {
            std::cout << "number of nodes: " << m_node_tags.size() << std::endl;
            std::cout << "number of values: " << current_band_values.size() << std::endl;
            throw std::runtime_error("Number of nodes and number of values are not the same. Abort.");
        }
        gmsh::view::addHomogeneousModelData(data_tag, 0, model_file_name, "NodeData", m_node_tags, current_band_values);
        const int   index_view             = gmsh::view::getIndex(data_tag);
        std::string name_object_visibility = "View[" + std::to_string(index_view) + "].Visible";
        gmsh::option::setNumber(name_object_visibility, 0);
        gmsh::view::write(data_tag, out_filename, true);
    }
}
