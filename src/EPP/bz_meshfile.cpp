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

#include <filesystem>

#include "BandStructure.h"
#include "gmsh.h"
#include "rapidcsv.h"

void bz_mesh_points::add_k_point(Vector3D<double> kpoint) { m_kpoints.push_back(kpoint); }

void bz_mesh_points::add_k_point(double k_x, double k_y, double k_z) { m_kpoints.push_back(Vector3D<double>(k_x, k_y, k_z)); }

void bz_mesh_points::read_mesh_from_csv() {
    m_node_tags.clear();
    m_kpoints.clear();
    m_nb_points = 0;
    std::cout << "Opening file " << m_filename << std::endl;
    rapidcsv::Document  doc(m_filename, rapidcsv::LabelParams(-1, -1), rapidcsv::SeparatorParams(' '));
    std::vector<double> k_x = doc.GetColumn<double>(0);
    std::vector<double> k_y = doc.GetColumn<double>(1);
    std::vector<double> k_z = doc.GetColumn<double>(2);

    std::size_t size_k_x = k_x.size();
    std::size_t size_k_y = k_y.size();
    std::size_t size_k_z = k_z.size();

    if (size_k_x != size_k_y || size_k_x != size_k_z) {
        throw std::runtime_error("Number of k-points in x, y and z are not the same. Abort.");
    }

    for (std::size_t index_k = 0; index_k < size_k_x; ++index_k) {
        add_k_point(k_x[index_k], k_y[index_k], k_z[index_k]);
    }

    std::cout << "Number of k-points read from csv file: " << size_k_x << std::endl;
}

void bz_mesh_points::read_mesh() {
    std::cout << "Opening file " << m_filename << std::endl;
    if (m_filename.substr(m_filename.find_last_of(".") + 1) == "csv") {
        read_mesh_from_csv();
        return;
    }

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

void bz_mesh_points::add_all_bands_on_mesh(const std::string&                             out_filename,
                                           const uepm::pseudopotential::BandStructure& my_band,
                                           int                                            nb_valence_bands_to_export,
                                           int                                            nb_conduction_bands_to_export) {
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 0);
    gmsh::model::add("bz_mesh");
    gmsh::open(m_filename);

    std::string model_file_name;
    gmsh::model::getCurrent(model_file_name);
    bool first_view = true;

    // We want to inverse the order of the valence bands so that the highest valence band is band 0
    int nb_valence_bands = my_band.is_soc_enabled() ? 8 : 4;
    int nb_bands         = my_band.get_number_of_bands();

    // Sanity checks
    if (nb_valence_bands_to_export < 0 || nb_valence_bands_to_export > nb_valence_bands) {
        nb_valence_bands_to_export = nb_valence_bands;
        std::cout << "Warning: nb_valence_bands_to_export is not valid. Setting it to " << nb_valence_bands << std::endl;
    }
    if (nb_conduction_bands_to_export < 0 || nb_conduction_bands_to_export > (nb_bands - nb_valence_bands)) {
        nb_conduction_bands_to_export = nb_bands - nb_valence_bands;
        std::cout << "Warning: nb_conduction_bands_to_export is not valid. Setting it to " << (nb_bands - nb_valence_bands) << std::endl;
    }

    int count = 0;
    if (nb_valence_bands_to_export != 0) {
        for (int index_band = nb_valence_bands - 1; index_band >= nb_valence_bands - nb_valence_bands_to_export; --index_band) {
            std::string         band_name   = "band_" + std::to_string(count);
            count++;
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
            gmsh::option::setNumber("PostProcessing.SaveMesh", (first_view) ? 1 : 0);  // Save mesh only once
            gmsh::view::write(data_tag, out_filename, true);
            first_view = false;
        }
    }

    if (nb_conduction_bands_to_export != 0) {
        for (int index_band = nb_valence_bands; index_band < nb_valence_bands + nb_conduction_bands_to_export; ++index_band) {
            std::string         band_name   = "band_" + std::to_string(count);
            count++;
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
            gmsh::option::setNumber("PostProcessing.SaveMesh", (first_view) ? 1 : 0);  // Save mesh only once
            gmsh::view::write(data_tag, out_filename, true);
            first_view = false;
        }
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
        std::cout << "band_values.size(): " << band_values.size() << std::endl;
        std::cout << "number_bands: " << number_bands << std::endl;
        std::cout << "m_node_tags.size(): " << m_node_tags.size() << std::endl;
        std::cout << "m_kpts.size(): " << m_kpoints.size() << std::endl;
        std::cout << "number_bands * m_node_tags.size(): " << number_bands * m_node_tags.size() << std::endl;
        throw std::runtime_error("band_values vector is not the same size as the number of bands times the number of nodes. Abort.");
    }
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 0);
    gmsh::model::add("bz_mesh");
    gmsh::open(m_filename);

    std::string model_file_name;
    gmsh::model::getCurrent(model_file_name);

    for (int index_band = 0; index_band < number_bands; ++index_band) {
        std::string         band_name = "band_" + std::to_string(index_band);
        std::vector<double> current_band_values(m_node_tags.size());
        for (std::size_t index_node = 0; index_node < m_node_tags.size(); ++index_node) {
            current_band_values[index_node] = band_values[index_node * number_bands + index_band];
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

/**
 * @brief Export band structure energies to csv files (one file per band).
 *
 * The band structure is assumed to be given in band_values vector under the following format:
 * band_values[index_k_point * number_of_bands + index_band] = energy of the band with index index_band at k-point with index index_k_point.
 *
 * @param out_filename
 * @param band_values
 */
void bz_mesh_points::export_bands_as_csv(const std::vector<double>& band_values, int number_bands) {
    // if (band_values.size() != number_bands * m_node_tags.size()) {
    //     throw std::runtime_error("band_values vector is not the same size as the number of bands times the number of nodes. Abort.");
    // }
    for (int index_band = 0; index_band < number_bands; ++index_band) {
        std::string         band_name = "band_" + std::to_string(index_band);
        std::vector<double> current_band_values(m_node_tags.size());
        for (std::size_t index_node = 0; index_node < m_node_tags.size(); ++index_node) {
            current_band_values[index_node] = band_values[index_node * number_bands + index_band];
        }
        std::ofstream band_file(band_name + ".csv");
        for (std::size_t index_node = 0; index_node < m_node_tags.size(); ++index_node) {
            band_file << current_band_values[index_node] << std::endl;
        }
        band_file.close();
    }
}
