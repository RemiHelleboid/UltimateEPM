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

#include "bz_mesh.hpp"

#include "gmsh.h"

void bz_mesh::add_k_point(Vector3D<double> kpoint) { m_kpoints.push_back(kpoint); }

void bz_mesh::add_k_point(double k_x, double k_y, double k_z) { m_kpoints.push_back(Vector3D<double>(k_x, k_y, k_z)); }

void bz_mesh::read_mesh() {
    std::cout << "Opening file " << m_filename << std::endl;
    gmsh::initialize();
    // gmsh::option::setNumber("General.Verbosity", 1);
    gmsh::model::add("bz_mesh");
    gmsh::open(m_filename);
    // std::vector<std::size_t> nodeTags;
    std::vector<double>      nodeCoords;
    std::vector<double>      nodeParams;
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

void bz_mesh::add_band_on_mesh(const std::string& band_name, const std::vector<double> band_values) {    gmsh::initialize();
    // gmsh::option::setNumber("General.Verbosity", 1);
    gmsh::initialize();
    gmsh::model::add("bz_mesh");
      gmsh::model::add("another model");
    gmsh::open(m_filename);
    int data_tag = gmsh::view::add(band_name);
    std::cout << "debug line : " << __LINE__ <<std::endl;
    gmsh::view::addHomogeneousModelData(data_tag, 0, "bz_mesh", "NodeData", m_node_tags, band_values);
    std::cout << "debug line : " << __LINE__ <<std::endl;
    const int   index_view             = gmsh::view::getIndex(data_tag);
    std::string name_object_visibility = "View[" + std::to_string(index_view) + "].Visible";
    gmsh::option::setNumber(name_object_visibility, 1);
    gmsh::view::write(data_tag, "band_mesh.msh", true);
    gmsh::finalize();

}