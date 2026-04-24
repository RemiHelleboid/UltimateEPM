/**
 * @file writer.cpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief
 * @version 0.1
 * @date 2021-07-23
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "writer.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "export_vector_to_csv.hpp"

namespace uepm {

namespace mesh {

void writer::export_as_msh(const mesh                     &myMesh,
                           const std::string              &filepath,
                           const std::vector<std::string> &datasets_to_export,
                           bool                            export_all_dataset) {
    /* Checking that the filepath is ok. */
    std::filesystem::path output_path(filepath);
    if (output_path.extension() != ".msh") {
        std::cout << "Output file : " << output_path << " has not the good extension for msh output." << std::endl;
        output_path.replace_extension(".msh");
        std::cout << "The new name of output file is : " << output_path << std::endl;
    }
    // std::cout << "ARMIN WILL EXPORT THE MESH." << std::endl;
    int bulk_element_type                  = 0;
    int contact_and_interface_element_type = 0;
    if (myMesh.m_dimension == 2) {
        bulk_element_type                  = 2;
        contact_and_interface_element_type = 1;
    } else if (myMesh.m_dimension == 3) {
        bulk_element_type                  = 4;
        contact_and_interface_element_type = 2;
    } else {
        assert("Dimension is neither 2 or 3.");
    }

    const int bulk_dimension              = myMesh.m_dimension;
    const int contact_interface_dimension = (myMesh.m_dimension == 3) ? 2 : 1;

    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 1);
    gmsh::model::add("mesh_discrete");

    std::vector<int>                      list_elements_tags(0);
    std::vector<std::size_t>              flatten_list_elements_vertices(0);
    std::vector<std::vector<std::size_t>> ListOfListRegionElementsTags(0);
    std::vector<std::vector<std::size_t>> ListOfListRegionFlattenElementVertices(0);

    std::vector<std::pair<int, int>> tags_gas_region;
    std::vector<int>                 ListTagsBulkGroups;

    for (const auto &bulk_region : myMesh.m_ListRegionsBulk) {
        std::string              region_name            = bulk_region.get_name();
        int                      entity_tag             = gmsh::model::addDiscreteEntity(bulk_dimension);
        std::vector<std::size_t> list_region_vertices   = bulk_region.get_unique_vertices_as_vector();
        std::size_t              number_vertices_region = list_region_vertices.size();
        std::vector<double>      flatten_vertices_coordinates(3 * number_vertices_region);
        std::vector<std::size_t> vertex_tag_list(number_vertices_region);
        for (std::size_t index_vertices = 0; index_vertices < number_vertices_region; ++index_vertices) {
            flatten_vertices_coordinates[3 * index_vertices]     = myMesh.m_ListVertices[list_region_vertices[index_vertices]].x();
            flatten_vertices_coordinates[3 * index_vertices + 1] = myMesh.m_ListVertices[list_region_vertices[index_vertices]].y();
            flatten_vertices_coordinates[3 * index_vertices + 2] = myMesh.m_ListVertices[list_region_vertices[index_vertices]].z();
            vertex_tag_list[index_vertices]                      = list_region_vertices[index_vertices] + 1;
        }
        gmsh::model::mesh::addNodes(bulk_dimension, entity_tag, vertex_tag_list, flatten_vertices_coordinates);
        std::vector<std::shared_ptr<element>> ListElements          = bulk_region.get_list_elements();
        std::vector<std::size_t>              ListRegionElementsTag = bulk_region.get_list_elements_index();
        std::for_each(ListRegionElementsTag.begin(), ListRegionElementsTag.end(), [](auto &&idx) { idx++; });
        std::vector<std::size_t> ListRegionFlattenedElementsVertices(0);
        for (const auto &elem : ListElements) {
            std::vector<std::size_t> list_element_vertices = elem->get_vertices_index();

            // int new_tag = list_elements_tags.size() + 1;
            // list_elements_tags.push_back(new_tag);
            // ListRegionElementsTag.push_back(new_tag);
            for (const auto &vtx_idx : list_element_vertices) {
                flatten_list_elements_vertices.push_back(vtx_idx + 1);
                ListRegionFlattenedElementsVertices.push_back(vtx_idx + 1);
            }
        }
        ListOfListRegionFlattenElementVertices.push_back(ListRegionFlattenedElementsVertices);
        ListOfListRegionElementsTags.push_back(ListRegionElementsTag);
        gmsh::model::mesh::addElements(bulk_dimension, entity_tag, {bulk_element_type}, {ListRegionElementsTag},
                                       {ListRegionFlattenedElementsVertices});

        int tag_phy_group = gmsh::model::addPhysicalGroup(bulk_dimension, {entity_tag});
        ListTagsBulkGroups.push_back(tag_phy_group);
        gmsh::model::setPhysicalName(bulk_dimension, tag_phy_group, region_name);
        if (region_name.find("Gas") != std::string::npos) {
            std::pair<int, int> gas_pair(bulk_dimension, tag_phy_group);
            tags_gas_region.push_back(gas_pair);
        }
    }

    for (const auto &interface_region : myMesh.m_ListRegionsInterface) {
        int                      entity_tag             = gmsh::model::addDiscreteEntity(contact_interface_dimension);
        std::vector<std::size_t> list_region_vertices   = interface_region.get_unique_vertices_as_vector();
        std::size_t              number_vertices_region = list_region_vertices.size();
        std::vector<double>      flatten_vertices_coordinates(3 * number_vertices_region);
        std::vector<std::size_t> vertex_tag_list(number_vertices_region);
        for (std::size_t index_vertices = 0; index_vertices < number_vertices_region; ++index_vertices) {
            flatten_vertices_coordinates[3 * index_vertices]     = myMesh.m_ListVertices[list_region_vertices[index_vertices]].x();
            flatten_vertices_coordinates[3 * index_vertices + 1] = myMesh.m_ListVertices[list_region_vertices[index_vertices]].y();
            flatten_vertices_coordinates[3 * index_vertices + 2] = myMesh.m_ListVertices[list_region_vertices[index_vertices]].z();
            vertex_tag_list[index_vertices]                      = list_region_vertices[index_vertices] + 1;
        }
        gmsh::model::mesh::addNodes(contact_interface_dimension, entity_tag, vertex_tag_list, flatten_vertices_coordinates);
        std::vector<std::shared_ptr<element>> ListElements = interface_region.get_list_elements();
        std::vector<std::size_t>              ListRegionElementsTag(0);
        std::vector<std::size_t>              ListRegionFlattenedElementsVertices(0);
        for (const auto &elem : ListElements) {
            int new_tag = list_elements_tags.size() + 1;
            list_elements_tags.push_back(new_tag);
            ListRegionElementsTag.push_back(new_tag);
            std::vector<std::size_t> list_element_vertices = elem->get_vertices_index();
            for (const auto &vtx_idx : list_element_vertices) {
                flatten_list_elements_vertices.push_back(vtx_idx + 1);
                ListRegionFlattenedElementsVertices.push_back(vtx_idx + 1);
            }
        }
        ListOfListRegionFlattenElementVertices.push_back(ListRegionFlattenedElementsVertices);
        ListOfListRegionElementsTags.push_back(ListRegionElementsTag);

        gmsh::model::mesh::addElements(contact_interface_dimension, entity_tag, {contact_and_interface_element_type},
                                       {ListRegionElementsTag}, {ListRegionFlattenedElementsVertices});

        std::string region_name   = interface_region.get_name();
        int         tag_phy_group = gmsh::model::addPhysicalGroup(contact_interface_dimension, {entity_tag});
        gmsh::model::setPhysicalName(contact_interface_dimension, tag_phy_group, region_name);
    }

    for (const auto &interface_region : myMesh.m_ListRegionsContact) {
        int                      entity_tag             = gmsh::model::addDiscreteEntity(contact_interface_dimension);
        std::vector<std::size_t> list_region_vertices   = interface_region.get_unique_vertices_as_vector();
        std::size_t              number_vertices_region = list_region_vertices.size();
        std::vector<double>      flatten_vertices_coordinates(3 * number_vertices_region);
        std::vector<std::size_t> vertex_tag_list(number_vertices_region);
        for (std::size_t index_vertices = 0; index_vertices < number_vertices_region; ++index_vertices) {
            flatten_vertices_coordinates[3 * index_vertices]     = myMesh.m_ListVertices[list_region_vertices[index_vertices]].x();
            flatten_vertices_coordinates[3 * index_vertices + 1] = myMesh.m_ListVertices[list_region_vertices[index_vertices]].y();
            flatten_vertices_coordinates[3 * index_vertices + 2] = myMesh.m_ListVertices[list_region_vertices[index_vertices]].z();
            vertex_tag_list[index_vertices]                      = list_region_vertices[index_vertices] + 1;
        }
        gmsh::model::mesh::addNodes(contact_interface_dimension, entity_tag, vertex_tag_list, flatten_vertices_coordinates);
        std::vector<std::shared_ptr<element>> ListElements = interface_region.get_list_elements();
        std::vector<std::size_t>              ListRegionElementsTag(0);
        std::vector<std::size_t>              ListRegionFlattenedElementsVertices(0);
        for (const auto &elem : ListElements) {
            int new_tag = list_elements_tags.size() + 1;
            list_elements_tags.push_back(new_tag);
            ListRegionElementsTag.push_back(new_tag);
            std::vector<std::size_t> list_element_vertices = elem->get_vertices_index();
            for (const auto &vtx_idx : list_element_vertices) {
                flatten_list_elements_vertices.push_back(vtx_idx + 1);
                ListRegionFlattenedElementsVertices.push_back(vtx_idx + 1);
            }
        }
        ListOfListRegionFlattenElementVertices.push_back(ListRegionFlattenedElementsVertices);
        ListOfListRegionElementsTags.push_back(ListRegionElementsTag);

        gmsh::model::mesh::addElements(contact_interface_dimension, entity_tag, {contact_and_interface_element_type},
                                       {ListRegionElementsTag}, {ListRegionFlattenedElementsVertices});

        std::string region_name   = interface_region.get_name();
        int         tag_phy_group = gmsh::model::addPhysicalGroup(contact_interface_dimension, {entity_tag});
        gmsh::model::setPhysicalName(contact_interface_dimension, tag_phy_group, region_name);
    }

    // The following line is commented to avoid a bug in the export of the views (data).
    // We should check if their is a way to remove the duplicate nodes without this bug.
    // gmsh::model::mesh::removeDuplicateNodes();

    gmsh::model::geo::synchronize();

    /* If no dataset to export, we just write the mesh. */
    // std::cout << "Datasets to export : " << datasets_to_export.size() << std::endl;
    // std::cout << "Export all dataset : " << export_all_dataset << std::endl;
    if (datasets_to_export.empty() && !export_all_dataset) {
        std::cout << "No dataset to export, only the mesh will be exported." << std::endl;
        gmsh::write(output_path);
    }

    /***  If necessary, we export the dataset, the mesh will be exported with the first dataset
     *   and then the others dataset will be written on the same file.
     ***/
    else {
        // std::cout << "Exporting datasets ... " << std::endl;
        std::vector<std::size_t> global_vertex_tag_list(myMesh.m_ListVertices.size());
        std::transform(myMesh.m_ListVertices.begin(), myMesh.m_ListVertices.end(), global_vertex_tag_list.begin(),
                       [&](const auto vtx) { return (vtx.get_index() + 1); });

        // Handling scalar data
        std::vector<std::string> list_name_dataset_to_export;
        if (export_all_dataset) {
            // std::All datasets will be exported << "All datasets will be exported" << std::endl;
            list_name_dataset_to_export = myMesh.get_all_functions_names();
        } else {
            list_name_dataset_to_export = datasets_to_export;
        }

        for (const auto &dataname : list_name_dataset_to_export) {
            // std::cout << "Exporting dataset :" << dataname << std::endl;
            if (!myMesh.scalar_dataset_exists(dataname)) {
                // std::cout << "Error: scalar datatset " << dataname << " does not exists, it can't be exported." << std::endl;
                continue;
            }
            auto my_function = myMesh.get_sp_scalar_function(dataname);
            if (my_function->get_location_type() == DataLocationType::vertex) {
                auto vertex_index_values_of_function = myMesh.get_vertices_index_value_of_scalar_function(dataname);
                int  data_tag                        = gmsh::view::add(dataname);

                std::vector<std::size_t> index_vtx(vertex_index_values_of_function.first);
                std::for_each(index_vtx.begin(), index_vtx.end(), [](auto &idx) { ++idx; });
                gmsh::view::addHomogeneousModelData(data_tag, 0, "mesh_discrete", "NodeData", index_vtx,
                                                    vertex_index_values_of_function.second);
                const int   index_view             = gmsh::view::getIndex(data_tag);
                std::string name_object_visibility = "View[" + std::to_string(index_view) + "].Visible";
                gmsh::option::setNumber(name_object_visibility, 1);
                gmsh::view::write(data_tag, output_path, true);
            } else if (my_function->get_location_type() == DataLocationType::cell) {
                auto cell_index_values_of_function = myMesh.get_cells_index_value_of_scalar_function(dataname);
                int  data_tag                      = gmsh::view::add(dataname);

                std::vector<std::size_t> index_vtx(cell_index_values_of_function.first);
                std::for_each(index_vtx.begin(), index_vtx.end(), [](auto &idx) { ++idx; });
                gmsh::view::addHomogeneousModelData(data_tag, 0, "mesh_discrete", "ElementData", index_vtx,
                                                    cell_index_values_of_function.second);
                const int   index_view             = gmsh::view::getIndex(data_tag);
                std::string name_object_visibility = "View[" + std::to_string(index_view) + "].Visible";
                gmsh::option::setNumber(name_object_visibility, 1);
                gmsh::view::write(data_tag, output_path, true);
            }
        }

        // Handling vector data
        for (auto &&dataname : list_name_dataset_to_export) {
            if (!myMesh.vector_dataset_exists(dataname)) {
                // std::cout << "Error: vector datatset " << dataname << " does not exists, it can't be exported." << std::endl;
                continue;
            }
            // std::cout << "Exporting: vector datatset " << dataname << std::endl;
            std::vector<vector3>     ScalarDataValues                = myMesh.get_all_vector_dataset_values(dataname);
            auto                     vertex_index_values_of_function = myMesh.get_vertices_index_value_of_vector_function(dataname);
            std::vector<std::size_t> index_vtx(vertex_index_values_of_function.first);
            std::for_each(index_vtx.begin(), index_vtx.end(), [](auto &idx) { ++idx; });

            std::vector<double> flattened_list_data;
            // flattened_list_data.reserve(3 * vertex_index_values_of_function.first.size());
            for (auto &&vector_data : vertex_index_values_of_function.second) {
                flattened_list_data.push_back(vector_data.x());
                flattened_list_data.push_back(vector_data.y());
                flattened_list_data.push_back(vector_data.z());
            }
            int data_tag = gmsh::view::add(dataname);
            gmsh::view::addHomogeneousModelData(data_tag, 0, "mesh_discrete", "NodeData", index_vtx, flattened_list_data);
            const int   index_view             = gmsh::view::getIndex(data_tag);
            std::string name_object_visibility = "View[" + std::to_string(index_view) + "].Visible";
            gmsh::option::setNumber(name_object_visibility, 0);
            gmsh::view::write(data_tag, output_path, true);
        }
    }
    // gmsh::fltk::run();
    gmsh::finalize();
    // std::cout << "FINAL : MSH FILE EXPORTED : " << output_path << std::endl;
}

}  // namespace mesh

}  // namespace uepm