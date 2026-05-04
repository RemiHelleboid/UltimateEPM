/**
 * @file msh_file.cpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief source code file for msh_file class.
 * @version 0.1
 * @date 2021-08-10
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "msh_file.hpp"

#include <algorithm>
#include <functional>
#include <map>
#include <numeric>
#include <regex>
#include <set>
#include <tuple>
#include <utility>

#include "dataset.hpp"
#include "element.hpp"
#include "element1d.hpp"
#include "element2d.hpp"
#include "element3d.hpp"
#include "export_vector_to_csv.hpp"
#include "file.hpp"
#include "gmsh.h"

namespace uepm {

namespace file {

std::string find_material_name_from_region_name(const std::string& region_name) {
    const std::string delimiter     = "_";
    std::string       material_name = region_name.substr(0, region_name.find(delimiter));
    return material_name;
}

msh_file::msh_file(const std::string& filepath) : file(filepath) {}

/**
 * @brief Open msh file
 *
 */
void msh_file::open_file() {
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 1);
    // IF_PLOG(plog::debug){
    //     gmsh::option::setNumber("General.Verbosity", 5);
    // }
    LOG_INFO << "OPEN MSH FILE : " << m_file_path << "\n";
    gmsh::open(m_file_path);
    std::string name;
    gmsh::model::getCurrent(name);
    int dimension = gmsh::model::getDimension();
    m_Mesh.set_dimension(dimension);
}

/**
 * @brief read the mesh data from msh file
 *
 * You shall not understand this function.
 *
 */
void msh_file::read_mesh() {
    this->open_file();
    int dimension = m_Mesh.get_dimension();
    m_Mesh.set_dimension(dimension);

    /* READING OF VERTICES */
    std::vector<std::size_t> nodeTags;
    std::vector<double>      nodeCoords;
    std::vector<double>      nodeParams;
    gmsh::model::mesh::reclassifyNodes();
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, -1, -1, false, false);
    std::size_t size_nodes_tags        = nodeTags.size();
    std::size_t size_nodes_coordinates = nodeCoords.size();
    if (size_nodes_coordinates != 3 * size_nodes_tags) {
        throw std::runtime_error("Number of coordinates is not 3 times the number of vertices. Abort.");
    }
    for (std::size_t index_vertex = 0; index_vertex < size_nodes_tags; ++index_vertex) {
        std::size_t node_tag = nodeTags[index_vertex];
        m_Mesh.add_vertex(nodeCoords[3 * index_vertex], nodeCoords[3 * index_vertex + 1], nodeCoords[3 * index_vertex + 2], node_tag - 1);
    }
    m_Mesh.sort_vertices();

    /* READING OF REGIONS AS PHYSICAL GROUPS */
    std::map<int, std::tuple<int, int>> entities_for_region;
    int                                 number_region = -1;
    std::vector<std::pair<int, int>>    PhysicalGroups;
    std::map<std::pair<int, int>, int>  PhysicalGroupsToDimTag;
    gmsh::model::getPhysicalGroups(PhysicalGroups);
    /* The order matter because the interface region will require the bulk regions,
       so the bulk region must be added first. */
    std::sort(PhysicalGroups.begin(), PhysicalGroups.end(), std::greater<>());
    int nb_no_name_group = 0;
    for (auto&& phy_grp : PhysicalGroups) {
        int index_new_region = ++number_region;
        int dim_phy_group    = phy_grp.first;
        int tag_phy_group    = phy_grp.second;

        std::vector<int> tags;
        gmsh::model::getEntitiesForPhysicalGroup(dim_phy_group, tag_phy_group, tags);
        PhysicalGroupsToDimTag.insert({{dim_phy_group, tag_phy_group}, index_new_region});
        LOG_INFO << "DIM PHY GROUP : " << dim_phy_group;
        LOG_INFO << "TAG PHY GROUP : " << tag_phy_group;
        LOG_INFO << "TAGS SIZE PHY GROUP : " << tags[0];

        std::string name_phy_group;
        gmsh::model::getPhysicalName(dim_phy_group, tag_phy_group, name_phy_group);

        std::string name_region_bulk    = name_phy_group;
        std::string name_region_contact = name_phy_group;
        if (name_phy_group.empty()) {
            name_region_bulk    = "Silicon_" + std::to_string(++nb_no_name_group);
            name_region_contact = "contact";
        }
        std::string material_name = find_material_name_from_region_name(name_phy_group);
        if (material_name.empty()) {
            material_name = "Silicon";
        }
        LOG_INFO << "READ REGION : " << name_region_bulk;
        /* Handling bulk region */
        if (dim_phy_group == dimension) {
            mesh::region_bulk new_bulk_region(dimension, name_region_bulk, index_new_region, material_name);
            m_Mesh.add_bulk_region(new_bulk_region);
        } else if (dim_phy_group == dimension - 1) {
            auto it_plus_sign = name_phy_group.find("+");
            /* Check if it is an interface region */
            if (it_plus_sign != std::string::npos) {
                std::string         region_bulk_0   = name_phy_group.substr(0, it_plus_sign);
                std::string         region_bulk_1   = name_phy_group.substr(it_plus_sign + 1, std::string::npos);
                const mesh::region* p_region_bulk_0 = m_Mesh.get_p_region(region_bulk_0);
                const mesh::region* p_region_bulk_1 = m_Mesh.get_p_region(region_bulk_1);
                if (p_region_bulk_0 != nullptr && p_region_bulk_1 != nullptr) {
                    int                    index_region_bulk_0 = p_region_bulk_0->get_index();
                    int                    index_region_bulk_1 = p_region_bulk_1->get_index();
                    mesh::region_interface new_interface_region(dimension,
                                                                name_phy_group,
                                                                index_new_region,
                                                                index_region_bulk_0,
                                                                index_region_bulk_1);
                    m_Mesh.add_interface_region(new_interface_region);
                } else {
                    std::cout << "REGION : " << name_phy_group << " cannot be added." << std::endl;
                    continue;
                    // throw std::invalid_argument("One bulk region of interface region was not found. Abort.
                    // ");
                }
            }
            /* Third case : contact region */
            else {
                mesh::region_contact new_contact_region(dimension, name_region_contact, index_new_region, -1);
                m_Mesh.add_contact_region(new_contact_region);
            }
        }
    }
    IF_PLOG(plog::verbose) { m_Mesh.print_regions_info(); }

    /* READING OF GEOMETRY */
    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities);
    std::string entity_name;

    for (auto e : entities) {
        // Dimension and tag of the entity:
        int region_index = -1;
        int dim          = e.first;
        int tag          = e.second;

        std::vector<int> physical_tags;
        gmsh::model::getPhysicalGroupsForEntity(dim, tag, physical_tags);
        if (physical_tags.size() > 1) {
            // for (auto&& phy_tag : physical_tags) {
            //     std::cout << phy_tag << ", ";
            // }
            // std::cout << std::endl;
            LOG_ERROR << "An entity belongs to multiple physical group. Fatal error.";
            // throw std::runtime_error("An element belongs to multiple physical group : error.");
        }
        if (physical_tags.empty()) {
            LOG_WARNING << "An entity does not belong to any physical group. It is skipped.";
            // throw std::runtime_error("An element belongs to no physical group : error.");
            continue;
        }
        auto MyRegionTag = PhysicalGroupsToDimTag.find({dim, physical_tags[0]});
        if (MyRegionTag != PhysicalGroupsToDimTag.end()) {
            region_index = MyRegionTag->second;
            LOG_INFO << " MY REGION TAG : " << region_index;
        } else {
            LOG_FATAL << " MY REGION TAG : " << region_index;
            throw std::runtime_error("Error when adding element to region. Region Unknown.");
        }
        gmsh::model::getEntityName(dim, tag, entity_name);
        mesh::region* p_current_region = m_Mesh.get_p_region(region_index);
        if (p_current_region == nullptr) {
            throw std::runtime_error("NO REGION");
        }

        /* READING OF ELEMENTS */
        std::vector<int>                      elemTypes;
        std::vector<std::vector<std::size_t>> elemTags;
        std::vector<std::vector<std::size_t>> elemNodeTags;

        gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);
        if (elemNodeTags.empty()) {
            continue;
        }

        std::vector<std::size_t> ElementNodeTags = elemNodeTags[0];
        std::vector<std::size_t> MyElementTag    = elemTags[0];
        LOG_DEBUG << "ElementTag size : " << MyElementTag.size();
        LOG_DEBUG << "ElementTag 0 : " << MyElementTag[0];


        if (dim == 1) {
            std::size_t element_index_1d = 0;
            for (std::size_t idx_vtx = 0; idx_vtx < ElementNodeTags.size(); idx_vtx += 2) {
                mesh::vertex* p_vtxA = m_Mesh.get_p_vertex(ElementNodeTags[idx_vtx] - 1);
                mesh::vertex* p_vtxB = m_Mesh.get_p_vertex(ElementNodeTags[idx_vtx + 1] - 1);
                if (p_vtxB == nullptr || p_vtxA == nullptr) {
                    throw std::runtime_error("Bad vertex index for construct element1d.");
                }
                std::shared_ptr<mesh::element1d> sp_element =
                    std::make_shared<mesh::element1d>(MyElementTag[element_index_1d++] - 1, p_vtxA, p_vtxB);
                sp_element->set_region_index(p_current_region->get_index());
                p_current_region->add_element(sp_element);
            }
        } else if (dim == 2) {
            std::size_t element_index_2d = 0;
            for (std::size_t idx_vtx = 0; idx_vtx < ElementNodeTags.size(); idx_vtx += 3) {
                mesh::vertex* p_vtxA = m_Mesh.get_p_vertex(ElementNodeTags[idx_vtx] - 1);
                mesh::vertex* p_vtxB = m_Mesh.get_p_vertex(ElementNodeTags[idx_vtx + 1] - 1);
                mesh::vertex* p_vtxC = m_Mesh.get_p_vertex(ElementNodeTags[idx_vtx + 2] - 1);
                if (p_vtxB == nullptr || p_vtxA == nullptr || p_vtxC == nullptr) {
                    std::cout << "One node not found : " << ElementNodeTags[idx_vtx] << std::endl;
                    std::cout << "One node not found : " << ElementNodeTags[idx_vtx + 1] << std::endl;
                    std::cout << "One node not found : " << ElementNodeTags[idx_vtx + 2] << std::endl;
                    throw std::runtime_error("NO " + std::to_string(ElementNodeTags[idx_vtx]));
                }
                std::shared_ptr<mesh::element2d> sp_element =
                    std::make_shared<mesh::element2d>(MyElementTag[element_index_2d++] - 1, p_vtxA, p_vtxB, p_vtxC);
                sp_element->set_region_index(p_current_region->get_index());
                p_current_region->add_element(sp_element);
                // LOG_INFO << "Add element2D (triangle) to " << index_mesh_region;
            }

        } else if (dim == 3) {
            std::size_t element_index_3d = 0;
            for (std::size_t idx_vtx = 0; idx_vtx < ElementNodeTags.size(); idx_vtx += 4) {
                mesh::vertex* p_vtxA = m_Mesh.get_p_vertex(ElementNodeTags[idx_vtx] - 1);
                mesh::vertex* p_vtxB = m_Mesh.get_p_vertex(ElementNodeTags[idx_vtx + 1] - 1);
                mesh::vertex* p_vtxC = m_Mesh.get_p_vertex(ElementNodeTags[idx_vtx + 2] - 1);
                mesh::vertex* p_vtxD = m_Mesh.get_p_vertex(ElementNodeTags[idx_vtx + 3] - 1);
                if (p_vtxB == nullptr || p_vtxA == nullptr || p_vtxC == nullptr || p_vtxD == nullptr) {
                    throw std::runtime_error("ERROR");
                }
                std::shared_ptr<mesh::element3d> sp_element =
                    std::make_shared<mesh::element3d>(MyElementTag[element_index_3d++] - 1, p_vtxA, p_vtxB, p_vtxC, p_vtxD);
                sp_element->set_region_index(p_current_region->get_index());
                p_current_region->add_element(sp_element);
            }
        } else {
            throw std::runtime_error("Dimension not supported.");
        }
        p_current_region->compute_unique_vertices();
    }
    IF_PLOG(plog::verbose) { m_Mesh.print_regions_info(); }
    gmsh::finalize();
    m_Mesh.build_search_tree();
}

// /**
//  * @brief Read the "views" from msh file
//  *
//  */
// void msh_file::read_states() {
//     this->open_file();
//     std::cout << "Read gmsh views ..." << std::endl;
//     std::vector<int> viewTags;
//     gmsh::view::getTags(viewTags);
//     for (auto&& tag : viewTags) {
//         const int   index_view  = gmsh::view::getIndex(tag);
//         std::string name_object = "View[" + std::to_string(index_view) + "].Name";
//         std::string name_view;
//         try {
//             gmsh::option::getString(name_object, name_view);
//         } catch (const std::exception& e) {
//             LOG_ERROR << "ERROR : A gmsh view could not be read. Skipping";
//             std::cerr << e.what() << '\n';
//         }
//         // std::cout << "NAME OF THE VIEW READ: " << name_view << std::endl;
//         std::string              type;
//         std::vector<std::size_t> tags;
//         double                   time;
//         int                      numComp;
//         std::vector<double>      data_view;
//         gmsh::view::getHomogeneousModelData(tag, 0, type, tags, data_view, time, numComp);
//         mesh::DataLocationType data_location_type = mesh::msh_to_armin_data_location_type.at(type);
//         LOG_INFO << "Read view type : " << type;
//         LOG_INFO << "Read view tag : " << tag;
//         LOG_INFO << "Number component :  " << numComp;
//         LOG_INFO << "Size vector data :  " << data_view.size();
//         LOG_INFO << "Number vertices : " << m_Mesh.get_nb_vertices();
//         std::transform(tags.begin(), tags.end(), tags.begin(), [&](auto vtx) { return --vtx; });

//         if (data_location_type == mesh::DataLocationType::vertex) {
//             if (numComp == 1) {
//                 m_Mesh.create_scalar_datasets_from_idx_vertex_and_values(name_view, data_location_type, tags, data_view);
//             }
//         } else if (data_location_type == mesh::DataLocationType::cell) {
//             if (numComp == 1) {
//                 // Scalar dataset defined on elements
//                 m_Mesh.create_scalar_datasets_from_idx_cells_and_values(name_view, data_location_type, tags, data_view);
//             }
//         }
//     }
//     gmsh::finalize();
//     m_Mesh.re_index_datasets();

//     // m_Mesh.add_scalar_data_to_all_vertices();
//     // m_Mesh.add_vector_data_to_all_vertices();
//     // m_Mesh.create_gradient_function("ElectrostaticPotential", "ArminElectricField");
// }

void msh_file::read_states(const std::vector<std::string>& list_dataset_to_import) {
    this->open_file();
    std::vector<int> viewTags;
    gmsh::view::getTags(viewTags);
    for (auto&& tag : viewTags) {
        const int   index_view  = gmsh::view::getIndex(tag);
        std::string name_object = "View[" + std::to_string(index_view) + "].Name";
        std::string name_view;
        try {
            gmsh::option::getString(name_object, name_view);
        } catch (const std::exception& e) {
            LOG_ERROR << "ERROR : A gmsh view could not be read. Skipping";
            std::cerr << e.what() << '\n';
        }

        if (!list_dataset_to_import.empty() &&
            std::find(list_dataset_to_import.begin(), list_dataset_to_import.end(), name_view) == list_dataset_to_import.end()) {
            // If the name of the view is not in the list of dataset to import, we skip it.
            // If the list is empty, we import all the views.
            continue;
        }
        std::string              type;
        std::vector<std::size_t> tags;
        double                   time;
        int                      numComp;
        std::vector<double>      data_view;
        gmsh::view::getHomogeneousModelData(tag, 0, type, tags, data_view, time, numComp);
        mesh::DataLocationType data_location_type = mesh::msh_to_armin_data_location_type.at(type);
        LOG_INFO << "Read view type : " << type;
        LOG_INFO << "Read view tag : " << tag;
        LOG_INFO << "Number component :  " << numComp;
        LOG_INFO << "Size vector data :  " << data_view.size();
        LOG_INFO << "Number vertices : " << m_Mesh.get_nb_vertices();
        std::transform(tags.begin(), tags.end(), tags.begin(), [&](auto vtx) { return --vtx; });

        if (data_location_type == mesh::DataLocationType::vertex) {
            if (numComp == 1) {
                m_Mesh.create_scalar_datasets_from_idx_vertex_and_values(name_view, data_location_type, tags, data_view);
            } else if (numComp == 2) {
                std::vector<mesh::vector3> data_view_vector2;
                for (std::size_t i = 0; i < data_view.size(); i += 2) {
                    data_view_vector2.push_back(mesh::vector3(data_view[i], data_view[i + 1], 0.0));
                }
                m_Mesh.create_vector_datasets_from_idx_vertex_and_values(name_view, data_location_type, tags, data_view_vector2);
            } else if (numComp == 3) {
                std::vector<mesh::vector3> data_view_vector3;
                for (std::size_t i = 0; i < data_view.size(); i += 3) {
                    data_view_vector3.push_back(mesh::vector3(data_view[i], data_view[i + 1], data_view[i + 2]));
                }
                m_Mesh.create_vector_datasets_from_idx_vertex_and_values(name_view, data_location_type, tags, data_view_vector3);
            }
        } else if (data_location_type == mesh::DataLocationType::cell) {
            if (numComp == 1) {
                // Scalar dataset defined on elements
                m_Mesh.create_scalar_datasets_from_idx_cells_and_values(name_view, data_location_type, tags, data_view);
            }
        }
    }
    gmsh::finalize();
    m_Mesh.re_index_datasets();

    m_Mesh.add_scalar_data_to_all_vertices();
    m_Mesh.add_vector_data_to_all_vertices();
}

}  // namespace file

}  // namespace uepm
