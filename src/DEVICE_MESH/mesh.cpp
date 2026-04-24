/**
 * @file mesh.cpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief
 * @version 0.1
 * @date 2021-07-10
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "mesh.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <utility>

#include <plog/Log.h>


// #include "config.h"
#include "export_vector_to_csv.hpp"
#include "grid_data.hpp"
#include "octree_node.hpp"
#include "quadtree_node.hpp"
#include "tree_node.hpp"
#include "vector.hpp"

#ifdef USE_OPENMP_ACCELERATION
#include <omp.h>
#endif

#pragma omp declare reduction(merge : std::vector <int> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp declare reduction(merge : std::vector <std::size_t> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp declare reduction(merge : std::vector <double> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

namespace uepm {

namespace mesh {

void mesh::reset_all() {
    m_dimension   = 0;
    m_nb_regions  = 0;
    m_nb_vertices = 0;
    m_ListVertices.clear();
    m_p_search_tree         = nullptr;
    m_transformation_matrix = no_transformation_matrix;
    m_translation_vector    = no_translation_vector;
    m_ListRegionsBulk.clear();
    m_ListRegionsInterface.clear();
    m_ListRegionsContact.clear();
    m_list_scalar_functions.clear();
    m_list_vector_functions.clear();
}

void mesh::add_vertex(double x, double y) {
    vertex vtx(m_ListVertices.size(), x, y, 0.0);
    m_ListVertices.push_back(vtx);
}

void mesh::add_vertex(double x, double y, double z) {
    const vertex vtx(m_ListVertices.size(), x, y, z);
    m_ListVertices.push_back(vtx);
}

void mesh::add_vertex(double x, double y, double z, std::size_t index_vertex) {
    const vertex vtx(index_vertex, x, y, z);
    m_ListVertices.push_back(vtx);
}

void mesh::sort_vertices() {
    auto sort_by_index = [&](const vertex &vtxA, const vertex &vtxB) { return (vtxA.get_index() < vtxB.get_index()); };
    std::sort(m_ListVertices.begin(), m_ListVertices.end(), sort_by_index);
}

bbox mesh::get_bounding_box() const {
    std::vector<double> X_coords(m_ListVertices.size());
    std::vector<double> Y_coords(m_ListVertices.size());
    std::vector<double> Z_coords(m_ListVertices.size());
    std::transform(m_ListVertices.begin(), m_ListVertices.end(), X_coords.begin(), [&](const auto &p_vtx) { return p_vtx.x(); });
    std::transform(m_ListVertices.begin(), m_ListVertices.end(), Y_coords.begin(), [&](const auto &p_vtx) { return p_vtx.y(); });
    std::transform(m_ListVertices.begin(), m_ListVertices.end(), Z_coords.begin(), [&](const auto &p_vtx) { return p_vtx.z(); });
    const double x_min = *std::min_element(X_coords.begin(), X_coords.end());
    const double x_max = *std::max_element(X_coords.begin(), X_coords.end());
    const double y_min = *std::min_element(Y_coords.begin(), Y_coords.end());
    const double y_max = *std::max_element(Y_coords.begin(), Y_coords.end());
    const double z_min = *std::min_element(Z_coords.begin(), Z_coords.end());
    const double z_max = *std::max_element(Z_coords.begin(), Z_coords.end());
    return bbox(x_min, x_max, y_min, y_max, z_min, z_max);
}

void mesh::build_search_tree() {
    // After function call
    auto         start             = std::chrono::high_resolution_clock::now();
    const double dilatation_factor = 1.005;
    if (m_dimension == 2) {
        LOG_INFO << "BUILDING QUADTREE ... ";
        bbox primary_bbox = get_bounding_box();
        primary_bbox.dilate(dilatation_factor);
        m_p_search_tree = std::make_unique<quadtree_node>(get_list_p_bulk_element(), primary_bbox, 0);
        LOG_INFO << "BUILDING QUADTREE : DONE.";
    } else if (m_dimension == 3) {
        LOG_INFO << "BUILDING OCTREE ... ";
        bbox primary_bbox = get_bounding_box();
        primary_bbox.dilate(dilatation_factor);
        bool is_root    = true;
        m_p_search_tree = std::make_unique<octree_node>(get_list_p_bulk_element(), primary_bbox, is_root);
        LOG_INFO << "BUILDING OCTREE : DONE.";
    }
    auto stop     = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    LOG_INFO << "Time to compute the search tree : " << duration.count() << " ms \n";
}

bbox mesh::compute_bounding_box() const {
    std::vector<double> X_coords(m_ListVertices.size());
    std::vector<double> Y_coords(m_ListVertices.size());
    std::vector<double> Z_coords(m_ListVertices.size());
    std::transform(m_ListVertices.begin(), m_ListVertices.end(), X_coords.begin(), [&](const auto &p_vtx) { return p_vtx.x(); });
    std::transform(m_ListVertices.begin(), m_ListVertices.end(), Y_coords.begin(), [&](const auto &p_vtx) { return p_vtx.y(); });
    std::transform(m_ListVertices.begin(), m_ListVertices.end(), Z_coords.begin(), [&](const auto &p_vtx) { return p_vtx.z(); });
    const double x_min = *std::min_element(X_coords.begin(), X_coords.end());
    const double x_max = *std::max_element(X_coords.begin(), X_coords.end());
    const double y_min = *std::min_element(Y_coords.begin(), Y_coords.end());
    const double y_max = *std::max_element(Y_coords.begin(), Y_coords.end());
    const double z_min = *std::min_element(Z_coords.begin(), Z_coords.end());
    const double z_max = *std::max_element(Z_coords.begin(), Z_coords.end());
    return bbox(x_min, x_max, y_min, y_max, z_min, z_max);
}

// Regions

const region &mesh::get_region_of_index(std::size_t search_index) const {
    auto same_index = [&](auto const &reg) { return reg.get_index() == search_index; };
    auto it_bulk    = std::find_if(m_ListRegionsBulk.begin(), m_ListRegionsBulk.end(), same_index);
    if (it_bulk != m_ListRegionsBulk.end()) {
        return *it_bulk;
    }
    auto it_contact = std::find_if(m_ListRegionsContact.begin(), m_ListRegionsContact.end(), same_index);
    if (it_contact != m_ListRegionsContact.end()) {
        return *it_contact;
    }
    auto it_interface = std::find_if(m_ListRegionsInterface.begin(), m_ListRegionsInterface.end(), same_index);
    if (it_interface != m_ListRegionsInterface.end()) {
        return *it_interface;
    }
    std::cerr << "No region found for index : " << search_index << ". Return nullptr." << std::endl;
    throw std::runtime_error("Region does not exist.");
}

region *mesh::get_p_region(std::size_t search_index) {
    region *p_region = nullptr;
    for (region &reg_bulk : m_ListRegionsBulk) {
        if (reg_bulk.get_index() == search_index) {
            p_region = &reg_bulk;
            return (p_region);
        }
    }
    for (region &reg_inter : m_ListRegionsInterface) {
        if (reg_inter.get_index() == search_index) {
            p_region = &reg_inter;
            return (p_region);
        }
    }
    for (region &reg_contact : m_ListRegionsContact) {
        if (reg_contact.get_index() == search_index) {
            p_region = &reg_contact;
            return (p_region);
        }
    }
    std::cerr << "No region found for index : " << search_index << ". Return nullptr." << std::endl;
    return nullptr;
}

const region_bulk *mesh::get_p_region_bulk(std::size_t search_index) const {
    for (const region_bulk &reg_bulk : m_ListRegionsBulk) {
        if (reg_bulk.get_index() == search_index) {
            return &reg_bulk;
        }
    }
    std::cerr << "No bulk region found for index : " << search_index << ". Return nullptr." << std::endl;
    return nullptr;
}

region_interface *mesh::get_p_interface_region_from_bulk_indices(std::size_t index_bulk_1, std::size_t index_bulk_2) {
    auto check_same_bulk_region_indices = [&](const region_interface &interface_reg) {
        const uint index_1     = interface_reg.get_bulk_0();
        const uint index_2     = interface_reg.get_bulk_1();
        const bool condition_1 = (index_1 == index_bulk_1 && index_2 == index_bulk_2);
        const bool condition_2 = (index_2 == index_bulk_1 && index_1 == index_bulk_2);
        return (condition_1 || condition_2);
    };
    auto it_interface_region = std::find_if(m_ListRegionsInterface.begin(), m_ListRegionsInterface.end(), check_same_bulk_region_indices);
    if (it_interface_region != m_ListRegionsInterface.end()) {
        return &(*it_interface_region);
    }
    return nullptr;
}

const region *mesh::get_p_region(const std::string &region_name) const {
    for (const region &reg_bulk : m_ListRegionsBulk) {
        if (reg_bulk.get_name() == region_name) {
            return &reg_bulk;
        }
    }
    for (const region &reg_inter : m_ListRegionsInterface) {
        if (reg_inter.get_name() == region_name) {
            return &reg_inter;
        }
    }
    for (const region &reg_contact : m_ListRegionsContact) {
        if (reg_contact.get_name() == region_name) {
            return &reg_contact;
        }
    }
    std::cerr << "No region found for index : " << region_name << ". Return nullptr." << std::endl;
    return (nullptr);
}

region *mesh::get_p_region(const std::string &region_name) {
    for (region &reg_bulk : m_ListRegionsBulk) {
        if (reg_bulk.get_name() == region_name) {
            return &reg_bulk;
        }
    }
    for (region &reg_inter : m_ListRegionsInterface) {
        if (reg_inter.get_name() == region_name) {
            return &reg_inter;
        }
    }
    for (region &reg_contact : m_ListRegionsContact) {
        if (reg_contact.get_name() == region_name) {
            return &reg_contact;
        }
    }
    std::cerr << "No region found for index : " << region_name << ". Return nullptr." << std::endl;
    return (nullptr);
}

std::vector<region const *> mesh::get_all_p_region() const {
    std::vector<region const *> vector_all_p_region;
    for (const auto &region : m_ListRegionsBulk) {
        vector_all_p_region.push_back(&region);
    }
    for (const auto &region : m_ListRegionsContact) {
        vector_all_p_region.push_back(&region);
    }
    for (const auto &region : m_ListRegionsInterface) {
        vector_all_p_region.push_back(&region);
    }
    return vector_all_p_region;
}

std::vector<region_bulk const *> mesh::get_all_p_bulk_region() const {
    std::vector<region_bulk const *> vector_all_p_region;
    for (const auto &region_bulk : m_ListRegionsBulk) {
        vector_all_p_region.push_back(&region_bulk);
    }
    return vector_all_p_region;
}

void mesh::remove_region(std::size_t index_region) {
    const auto *p_region_to_remove = get_p_region(index_region);
    std::cout << "REMOVING REGION : " << p_region_to_remove->get_name() << std::endl;
    const auto region_type = p_region_to_remove->get_region_type();
    if (region_type == RegionType::bulk) {
        LOG_DEBUG << "REMOVE BULK REGION";
        std::erase_if(m_ListRegionsBulk, [&](auto p_region) { return p_region.get_index() == index_region; });
    } else if (region_type == RegionType::interface) {
        LOG_DEBUG << "REMOVE INTERFACE REGION";
        std::erase_if(m_ListRegionsInterface, [&](auto p_region) { return p_region.get_index() == index_region; });
    } else if (region_type == RegionType::contact) {
        LOG_DEBUG << "REMOVE CONTACT REGION";
        std::erase_if(m_ListRegionsContact, [&](auto p_region) { return p_region.get_index() == index_region; });
    }
}

void mesh::compute_interface_region_betwwen_two_bulks(std::size_t                     index_bulk_1,
                                                      std::size_t                     index_bulk_2,
                                                      const std::vector<std::size_t> &list_index_element_region_2_to_check) {
    std::vector<std::size_t> list_interface_elements;
    std::vector<std::size_t> list_potential_interface_element_region_1;

    auto *p_region_1 = get_p_region(index_bulk_1);
    auto *p_region_2 = get_p_region(index_bulk_2);

    LOG_DEBUG << "COMPUTE THE INTEFRACE BETWEEN : " << p_region_1->get_name() << " AND  " << p_region_2->get_name();
    const std::string new_region_interface = p_region_1->get_name() + "+" + p_region_2->get_name();
    region_interface  new_interface_region(2, new_region_interface, get_nb_regions(), index_bulk_1, index_bulk_2);

    const auto               set_region_1_vertices = p_region_1->get_unique_vertices();
    const auto               set_region_2_vertices = p_region_2->get_unique_vertices();
    std::vector<std::size_t> common_vertices;

    std::set_intersection(set_region_1_vertices.begin(), set_region_1_vertices.end(), set_region_2_vertices.begin(),
                          set_region_2_vertices.end(), std::back_inserter(common_vertices));

    std::set set_common_vertices(common_vertices.begin(), common_vertices.end());

    std::vector<vector3> list_vtx_interface{};
    for (auto &&vtx_idx : common_vertices) {
        list_vtx_interface.push_back(*get_p_vertex(vtx_idx));
    }
    utils::export_vector_postion_to_csv("INTERFACE.csv", "", list_vtx_interface);
    std::vector<sp_element> element_region_1_to_check(list_index_element_region_2_to_check.size());
    std::transform(list_index_element_region_2_to_check.begin(), list_index_element_region_2_to_check.end(),
                   element_region_1_to_check.begin(), [&](auto &&idx_elem) { return p_region_2->get_p_element(idx_elem); });

    std::size_t counter_new_element = 0;
    for (auto &&idx_element : list_index_element_region_2_to_check) {
        auto                     sptr_element     = p_region_2->get_p_element(idx_element);
        auto                     list_vtx_element = sptr_element->get_vertices_index();
        std::set                 set_vtx_element(list_vtx_element.begin(), list_vtx_element.end());
        std::vector<std::size_t> interfaces_vtx;
        std::set_intersection(set_vtx_element.begin(), set_vtx_element.end(), set_common_vertices.begin(), set_common_vertices.end(),
                              std::back_inserter(interfaces_vtx));
        if (interfaces_vtx.size() == 3) {
            std::cout << "\rCounter detected element interface : " << ++counter_new_element << std::flush;
            vertex                    *p_vtx_1               = get_p_vertex(interfaces_vtx[0]);
            vertex                    *p_vtx_2               = get_p_vertex(interfaces_vtx[1]);
            vertex                    *p_vtx_3               = get_p_vertex(interfaces_vtx[2]);
            std::shared_ptr<element2d> new_interface_element = std::make_shared<element2d>(counter_new_element, p_vtx_1, p_vtx_2, p_vtx_3);
            new_interface_region.add_element(new_interface_element);
            list_interface_elements.push_back(idx_element);
        }
    }
    new_interface_region.compute_unique_vertices();
    new_interface_region.print_info();
    this->add_interface_region(new_interface_region);
}

std::vector<std::size_t> mesh::get_idx_bulk_elements_adjacent_to_contact_region(const std::string &region_name) const {
    std::vector<std::size_t> list_idx_bulk_elements_adjacent_to_contact_region;

    auto *p_region = get_p_region(region_name);
    if (!p_region) {
        std::cerr << "The region " << region_name << " does not exist." << std::endl;
        return list_idx_bulk_elements_adjacent_to_contact_region;
    }
    if (p_region->get_region_type() != RegionType::contact) {
        std::cerr << "The region " << region_name << " is not a contact region." << std::endl;
        return list_idx_bulk_elements_adjacent_to_contact_region;
    }

    const std::set<unsigned int> list_vertex_contact_region = p_region->get_unique_vertices();
    auto                         list_bulk_elements         = get_list_bulk_element();
    std::size_t                  nb_bulk_elements           = list_bulk_elements.size();
    for (std::size_t index_bulk_element = 0; index_bulk_element < nb_bulk_elements; ++index_bulk_element) {
        auto                      p_element           = list_bulk_elements[index_bulk_element];
        const auto                list_vertex_element = p_element->get_vertices_index();
        std::set<unsigned int>    set_vertex_element(list_vertex_element.begin(), list_vertex_element.end());
        std::vector<unsigned int> common_vertices;
        std::set_intersection(list_vertex_contact_region.begin(), list_vertex_contact_region.end(), set_vertex_element.begin(),
                              set_vertex_element.end(), std::back_inserter(common_vertices));
        if (common_vertices.size() > 0) {
            list_idx_bulk_elements_adjacent_to_contact_region.push_back(index_bulk_element);
        }
    }
    return list_idx_bulk_elements_adjacent_to_contact_region;
}

void mesh::compare_bulks_region_with_bounding_box() const {
    for (auto &&p_region : get_all_p_bulk_region()) {
        bbox   region_bbox        = p_region->compute_bounding_box();
        double bbox_volume        = region_bbox.get_surface();
        double real_volume_region = compute_region_volume(p_region->get_name());
        std::cout << "REGION NAME      : " << p_region->get_name() << std::endl;
        std::cout << "REGION VOLUME      : " << real_volume_region << std::endl;
        std::cout << "REGION BBOX VOLUME : " << bbox_volume << std::endl;
        std::cout << "RELATIVE ERROR : " << 100.0 * fabs(real_volume_region - bbox_volume) / real_volume_region << "%" << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
    }
}

//  Element
std::vector<sp_element> mesh::get_list_bulk_element() const {
    std::vector<sp_element> list_bulk_elements;
    for (const auto &bulk_region : m_ListRegionsBulk) {
        std::vector<sp_element> list_sp_elements = bulk_region.get_list_elements();
        list_bulk_elements.insert(list_bulk_elements.end(), list_sp_elements.begin(), list_sp_elements.end());
    }
    return list_bulk_elements;
}

std::vector<element *> mesh::get_list_p_bulk_element() const {
    const auto             list_sp_bulk_element = get_list_bulk_element();
    std::vector<element *> list_p_elements(list_sp_bulk_element.size());
    std::transform(list_sp_bulk_element.begin(), list_sp_bulk_element.end(), list_p_elements.begin(),
                   [](const auto &sp_elem) { return sp_elem.get(); });
    return list_p_elements;
}

void mesh::transfer_element_to_other_region(sp_element p_element, region *origin_region, region *new_region) {
    auto s_ptr_element = origin_region->get_p_element(p_element->get_index());
    new_region->add_element(s_ptr_element);
    origin_region->remove_element(s_ptr_element->get_index());
    new_region->compute_unique_vertices();
    origin_region->compute_unique_vertices();
}

void mesh::transfer_elements_to_other_region(std::vector<std::size_t> list_element_indexes, region *origin_region, region *new_region) {
    for (auto &&idx_element : list_element_indexes) {
        auto s_ptr_element = origin_region->get_p_element(idx_element);
        new_region->add_element(s_ptr_element);
        origin_region->remove_element(idx_element);
    }
    new_region->compute_unique_vertices();
    origin_region->compute_unique_vertices();
}

//  Datasets

std::vector<sp_scalar_dataset> mesh::get_list_scalar_datasets() const {
    std::vector<sp_scalar_dataset> list_sp_scalar_dataset;
    for (auto &&scalar_function : m_list_scalar_functions) {
        std::vector<sp_scalar_dataset> list_function_dataset = scalar_function->get_list_sp_datasets();
        list_sp_scalar_dataset.insert(list_sp_scalar_dataset.end(), list_function_dataset.begin(), list_function_dataset.end());
    }
    return list_sp_scalar_dataset;
}

std::vector<sp_vector_dataset> mesh::get_list_vector_datasets() const {
    std::vector<sp_vector_dataset> list_sp_vector_dataset;
    for (auto &&vector_function : m_list_vector_functions) {
        std::vector<sp_vector_dataset> list_function_dataset = vector_function->get_list_sp_datasets();
        list_sp_vector_dataset.insert(list_sp_vector_dataset.end(), list_function_dataset.begin(), list_function_dataset.end());
    }
    return list_sp_vector_dataset;
}

void mesh::create_scalar_datasets_from_idx_vertex_and_values(const std::string              &dataset_name,
                                                             DataLocationType                data_location_type,
                                                             const std::vector<std::size_t> &list_index_vertices,
                                                             const std::vector<double>      &data_values) {
    LOG_DEBUG << "CREATING SCALAR DATASET : " << dataset_name;
    for (const auto &region : m_ListRegionsBulk) {
        LOG_DEBUG << "WORKING ON REGION: " << region.get_name();
        const int                index_region_validity  = region.get_index();
        const int                index_dataset          = get_total_number_dataset() + 1;
        const DataType           Dataset_datatype       = DataType::scalar;
        const int                data_dimension         = 1;
        const auto               region_unique_vertices = region.get_unique_vertices_as_vector();
        std::vector<double>      dataset_values{};
        std::vector<std::size_t> dataset_index_element{};
        dataset_values.reserve(data_values.size());
        dataset_index_element.reserve(data_values.size());
        // std::size_t first_element_index =
        //     std::distance(list_index_vertices.begin(),
        //                   std::find(list_index_vertices.begin(), list_index_vertices.end(), region_unique_vertices[0]));
        bool right_region = true;
        for (auto &&vtx_index : region_unique_vertices) {
            // If a vertex of the region is not in the dataset indexes, the dataset is not defined on this region.
            // We then just exit the function.
            auto it_value_idx = std::find(list_index_vertices.begin(), list_index_vertices.end(), vtx_index);
            if (it_value_idx == list_index_vertices.end()) {
                LOG_DEBUG << "DATASET : " << dataset_name << " NOT DEFINED ON REGION: " << region.get_name();
                right_region = false;
                break;
            }
            std::size_t index_value = std::distance(list_index_vertices.begin(), it_value_idx);
            dataset_values.push_back(data_values[index_value]);
            dataset_index_element.push_back(index_value);
        }
        if (!right_region) {
            continue;
        }
        sp_scalar_dataset NewDataset =
            std::make_shared<dataset<double>>(dataset_name, index_dataset, index_region_validity, dataset_values, dataset_index_element,
                                              Dataset_datatype, data_location_type, data_dimension);
        add_scalar_dataset(NewDataset);
        add_scalar_data_to_vertices(*NewDataset);
    }
}
void mesh::create_vector_datasets_from_idx_vertex_and_values(const std::string              &dataset_name,
                                                             DataLocationType                data_location_type,
                                                             const std::vector<std::size_t> &list_index_vertices,
                                                             const std::vector<vector3>     &data_values) {
    LOG_DEBUG << "CREATING VECTOR DATASET : " << dataset_name;
    for (const auto &region : m_ListRegionsBulk) {
        LOG_DEBUG << "WORKING ON REGION: " << region.get_name();
        const int                index_region_validity  = region.get_index();
        const int                index_dataset          = get_total_number_dataset() + 1;
        const DataType           Dataset_datatype       = DataType::vector;
        const int                data_dimension         = m_dimension;
        const auto               region_unique_vertices = region.get_unique_vertices_as_vector();
        std::vector<vector3>     dataset_values{};
        std::vector<std::size_t> dataset_index_element{};
        dataset_values.reserve(data_values.size());
        dataset_index_element.reserve(data_values.size());
        // std::size_t first_element_index =
        //     std::distance(list_index_vertices.begin(),
        //                   std::find(list_index_vertices.begin(), list_index_vertices.end(), region_unique_vertices[0]));
        bool right_region = true;
        for (auto &&vtx_index : region_unique_vertices) {
            // If a vertex of the region is not in the dataset indexes, the dataset is not defined on this region.
            // We then just exit the function.
            auto it_value_idx = std::find(list_index_vertices.begin(), list_index_vertices.end(), vtx_index);
            if (it_value_idx == list_index_vertices.end()) {
                LOG_DEBUG << "DATASET : " << dataset_name << " NOT DEFINED ON REGION: " << region.get_name();
                right_region = false;
                break;
            }
            std::size_t index_value = std::distance(list_index_vertices.begin(), it_value_idx);
            dataset_values.push_back(data_values[index_value]);
            dataset_index_element.push_back(index_value);
        }
        if (!right_region) {
            continue;
        }
        sp_vector_dataset NewDataset =
            std::make_shared<dataset<vector3>>(dataset_name, index_dataset, index_region_validity, dataset_values, dataset_index_element,
                                               Dataset_datatype, data_location_type, data_dimension);
        add_vector_dataset(NewDataset);
        add_vector_data_to_vertices(*NewDataset);
    }
}

void mesh::create_scalar_datasets_from_idx_cells_and_values(const std::string              &dataset_name,
                                                            DataLocationType                data_location_type,
                                                            const std::vector<std::size_t> &list_index_cells,
                                                            const std::vector<double>      &data_values) {
    for (const auto &region : m_ListRegionsBulk) {
        LOG_DEBUG << "WORKING ON REGION: " << region.get_name();
        const int           index_region_validity = region.get_index();
        const int           index_dataset         = get_total_number_dataset() + 1;
        const DataType      Dataset_datatype      = DataType::scalar;
        const int           data_dimension        = 1;
        const auto          region_element_idx    = region.get_list_elements_index();
        std::vector<double> dataset_values{};
        dataset_values.reserve(data_values.size());
        std::vector<std::size_t> dataset_index_element{};
        dataset_index_element.reserve(data_values.size());
        bool right_region = true;

        std::size_t counter = 0;
        std::size_t first_element_index =
            std::distance(list_index_cells.begin(), std::find(list_index_cells.begin(), list_index_cells.end(), region_element_idx[0]));
        for (auto &&element_index : region_element_idx) {
            // If a vertex of the region is not in the dataset indexes, the dataset is not defined on this region.
            // We then just exit the function.
            auto it_value_idx =
                std::find(list_index_cells.begin() + first_element_index + counter - 1, list_index_cells.end(), element_index);
            counter++;
            if (it_value_idx == list_index_cells.end()) {
                LOG_DEBUG << "DATASET : " << dataset_name << " NOT DEFINED ON REGION: " << region.get_name();
                LOG_DEBUG << region.get_name() << "  IDX ELEMENTS : " << std::distance(list_index_cells.begin(), it_value_idx);
                right_region = false;
                break;
            }

            dataset_index_element.push_back(element_index);
            dataset_values.push_back(data_values[std::distance(list_index_cells.begin(), it_value_idx)]);
        }
        if (!right_region) {
            continue;
        }
        sp_scalar_dataset NewDataset =
            std::make_shared<dataset<double>>(dataset_name, index_dataset, index_region_validity, dataset_values, dataset_index_element,
                                              Dataset_datatype, data_location_type, data_dimension);
        add_scalar_dataset(NewDataset);
        add_scalar_data_to_elements(*NewDataset);
    }
}

void mesh::add_scalar_dataset(sp_scalar_dataset new_dataset) {
    const std::string  dataset_name = new_dataset->get_name();
    sp_scalar_function sp_function  = get_sp_scalar_function(dataset_name);
    if (sp_function == nullptr) {
        DataLocationType   data_location_type = new_dataset->get_data_location_type();
        sp_scalar_function sp_new_function    = std::make_shared<function<double>>(dataset_name, DataType::scalar, data_location_type);
        m_list_scalar_functions.push_back(sp_new_function);
        sp_function = sp_new_function;
    }
    sp_function->add_dataset(new_dataset);
}

void mesh::add_vector_dataset(sp_vector_dataset new_dataset) {
    const std::string  dataset_name = new_dataset->get_name();
    sp_vector_function sp_function  = get_sp_vector_function(dataset_name);
    if (sp_function == nullptr) {
        DataLocationType   data_location_type = new_dataset->get_data_location_type();
        sp_vector_function sp_new_function    = std::make_shared<function<vector3>>(dataset_name, DataType::vector, data_location_type);
        m_list_vector_functions.push_back(sp_new_function);
        sp_function = sp_new_function;
    }
    sp_function->add_dataset(new_dataset);
}

void mesh::add_scalar_data_to_vertices(const dataset<double> &dtset) {
    if (dtset.get_data_location_type() != DataLocationType::vertex) {
        LOG_WARNING << "DATASET " << dtset.get_name() << " IS DEFINED ON ELEMENTS, CANNOT BE ADD TO VERTICES.";
        return;
    }
    const std::string                          dataset_name     = dtset.get_name();
    const std::vector<double>                  dataset_values   = dtset.get_values();
    const std::shared_ptr<std::vector<double>> p_dataset_values = dtset.get_p_values();
    region                                    *p_myregion       = get_p_region(dtset.get_index_region_validity());
    if (p_myregion == nullptr) {
        return;
    }
    std::vector<std::size_t> list_vertices = p_myregion->get_unique_vertices_as_vector();
    if (dataset_values.size() != list_vertices.size()) {
        std::cerr << "Mismatch between number of vertices and number of values for the following dataset :" << std::endl;
        std::cerr << "Name            : " << dataset_name << std::endl;
        std::cerr << "Region validity : " << dtset.get_index_region_validity() << std::endl;
        std::cerr << "The region will not be added to the mesh model.\n " << std::endl;
        return;
    }
    std::size_t value_counter = 0;
    for (const auto &vtx_idx : list_vertices) {
        m_ListVertices[vtx_idx].add_scalar_data(dataset_name, &(p_dataset_values->at(value_counter)));
        ++value_counter;
    }
}

void mesh::add_scalar_data_to_elements(const dataset<double> &dtset) {
    if (dtset.get_data_location_type() != DataLocationType::cell) {
        LOG_WARNING << "DATASET " << dtset.get_name() << " IS DEFINED ON ELEMENTS, CANNOT BE ADD TO VERTICES.";
        return;
    }
    const std::string                          dataset_name     = dtset.get_name();
    const std::vector<double>                  dataset_values   = dtset.get_values();
    const std::shared_ptr<std::vector<double>> p_dataset_values = dtset.get_p_values();
    // const std::vector<std::size_t>            &list_elements_idx = dtset.get_index_geometry_elements();
    region *p_myregion = get_p_region(dtset.get_index_region_validity());
    if (p_myregion == nullptr) {
        return;
    }
    auto list_p_elements = p_myregion->get_list_elements();
    if (dataset_values.size() != list_p_elements.size()) {
        std::cerr << "Missmatch between number of elements and number of values for the following dataset :" << std::endl;
        std::cerr << "Name            : " << dataset_name << std::endl;
        std::cerr << "Region validity : " << dtset.get_index_region_validity() << std::endl;
        std::cerr << "The region will not be added to the mesh model.\n " << std::endl;
        return;
    }
    for (std::size_t idx_element = 0; idx_element < list_p_elements.size(); ++idx_element) {
        list_p_elements[idx_element]->add_scalar_data(dataset_name, &(p_dataset_values->at(idx_element)));
    }
}

void mesh::add_doping_concentration_to_vertices(const std::string &doping_fieldname) {
    for (auto &&vtx : m_ListVertices) {
        vtx.set_doping_concentration(vtx.get_scalar_data(doping_fieldname));
    }
}

void mesh::add_electric_field_to_vertices(const std::string &electric_field_fieldname, double factor) {
    for (auto &&vtx : m_ListVertices) {
        vtx.set_electric_field(factor * vtx.get_vector_data(electric_field_fieldname));
    }
}

void mesh::add_diffusion_gradient_to_vertices(const std::string &e_gradient_fieldname, const std::string &h_gradient_fieldname) {
    for (auto &&vtx : m_ListVertices) {
        vtx.set_e_grad_diffusion(vtx.get_vector_data(e_gradient_fieldname));
        vtx.set_h_grad_diffusion(vtx.get_vector_data(h_gradient_fieldname));
    }
}

void mesh::add_space_charge_to_vertices(const std::string &space_charge_fieldname) {
    for (auto &&vtx : m_ListVertices) {
        // DEBUG
        double space_charge = vtx.get_scalar_data(space_charge_fieldname);
        vtx.set_space_charge(space_charge);
    }
}

void mesh::add_charge_density_to_vertices(const std::string &electron_charge_density_fieldname,
                                          const std::string &hole_charge_density_fieldname) {
    for (auto &&vtx : m_ListVertices) {
        vtx.set_charge_density(vtx.get_scalar_data(hole_charge_density_fieldname) - vtx.get_scalar_data(electron_charge_density_fieldname));
    }
}

void mesh::add_vector_data_to_vertices(const dataset<vector3> &dtset) {
    if (dtset.get_data_location_type() != DataLocationType::vertex) {
        LOG_WARNING << "DATASET " << dtset.get_name() << " IS DEFINED ON ELEMENTS, CANNOT BE ADD TO VERTICES.";
    }
    const std::string                     dataset_name     = dtset.get_name();
    std::shared_ptr<std::vector<vector3>> p_dataset_values = dtset.get_p_values();
    region                               *p_myregion       = get_p_region(dtset.get_index_region_validity());
    if (p_myregion == nullptr) {
        return;
    }
    std::vector<std::size_t> list_vertices = p_myregion->get_unique_vertices_as_vector();

    std::size_t value_counter = 0;
    for (auto &&vtx_idx : list_vertices) {
        m_ListVertices[vtx_idx].add_vector_data(dataset_name, &(p_dataset_values->at(value_counter)));
        ++value_counter;
    }
}

void mesh::add_vector_data_to_elements(const dataset<vector3> &dtset) {
    if (dtset.get_data_location_type() != DataLocationType::cell) {
        LOG_WARNING << "DATASET " << dtset.get_name() << " IS DEFINED ON ELEMENTS, CANNOT BE ADD TO VERTICES.";
    }
    const std::string                     dataset_name     = dtset.get_name();
    std::shared_ptr<std::vector<vector3>> p_dataset_values = dtset.get_p_values();
    region                               *p_myregion       = get_p_region(dtset.get_index_region_validity());
    if (p_myregion == nullptr) {
        return;
    }
    auto list_p_elements = p_myregion->get_list_elements();
    for (std::size_t idx_element = 0; idx_element < list_p_elements.size(); ++idx_element) {
        list_p_elements[idx_element]->add_vector_data(dataset_name, &(p_dataset_values->at(idx_element)));
    }
}

void mesh::add_scalar_data_to_all_vertices() {
    for (const auto &dtset : get_list_scalar_datasets()) {
        if (dtset->get_data_location_type() == DataLocationType::vertex) {
            add_scalar_data_to_vertices(*dtset);
        }
    }
    add_doping_concentration_to_vertices("DopingConcentration");
    add_space_charge_to_vertices("SpaceCharge");
    add_charge_density_to_vertices("eDensity", "hDensity");
}

void mesh::add_vector_data_to_all_vertices() {
    for (const auto &dtset : get_list_vector_datasets()) {
        if (dtset->get_data_location_type() == DataLocationType::vertex) {
            add_vector_data_to_vertices(*dtset);
        }
    }
    add_electric_field_to_vertices("ElectricFieldV");
}

vertex *mesh::get_p_vertex_safe(std::size_t index) {
    auto it_vtx = std::find_if(m_ListVertices.begin(), m_ListVertices.end(), [&](const vertex &vtx) { return (vtx.get_index() == index); });
    return it_vtx != m_ListVertices.end() ? &(*it_vtx) : nullptr;
}

double mesh::get_scalar_data_at_vertex(std::size_t idx_vertex, const std::string &dataset_name) const {
    return m_ListVertices[idx_vertex].get_scalar_data(dataset_name);
}

double mesh::get_scalar_data_at_element(std::size_t idx_element, const std::string &dataset_name) const {
    return get_list_bulk_element()[idx_element]->get_scalar_data(dataset_name);
}

// TO CHANGE
std::vector<double> mesh::get_all_scalar_dataset_values(const std::string &dataset_name) const {
    std::vector<double> values(m_ListVertices.size());
    std::transform(m_ListVertices.begin(), m_ListVertices.end(), values.begin(),
                   [&](const vertex &vtx) { return vtx.get_scalar_data(dataset_name); });
    return (values);
}

std::pair<std::vector<std::size_t>, std::vector<double>> mesh::get_vertices_index_value_of_scalar_function(
    const std::string &function_name) const {
    std::vector<std::size_t> list_global_vertex_index{};
    std::vector<double>      list_values{};
    sp_scalar_function       shared_pointer_function = get_sp_scalar_function(function_name);
    // for (auto const p_dataset : shared_pointer_function->get_list_sp_datasets()) {
    for (auto &&p_dataset : shared_pointer_function->get_list_sp_datasets()) {
        auto index_region       = p_dataset->get_index_region_validity();
        auto region_dataset     = get_region_of_index(index_region);
        auto list_vertex_region = region_dataset.get_unique_vertices_as_vector();
        list_global_vertex_index.insert(list_global_vertex_index.end(), list_vertex_region.begin(), list_vertex_region.end());
        std::vector<double> dataset_values(list_vertex_region.size());
        std::transform(list_vertex_region.begin(), list_vertex_region.end(), dataset_values.begin(),
                       [&](auto const vtx_idx) { return get_scalar_data_at_vertex(vtx_idx, function_name); });
        list_values.insert(list_values.end(), dataset_values.begin(), dataset_values.end());
    }

    return {list_global_vertex_index, list_values};
}

std::pair<std::vector<std::size_t>, std::vector<double>> mesh::get_cells_index_value_of_scalar_function(
    const std::string &function_name) const {
    std::vector<std::size_t> list_global_cells_index{};
    std::vector<double>      list_values{};
    sp_scalar_function       shared_pointer_function = get_sp_scalar_function(function_name);
    for (auto &&p_dataset : shared_pointer_function->get_list_sp_datasets()) {
        auto index_region           = p_dataset->get_index_region_validity();
        auto region_dataset         = get_region_of_index(index_region);
        auto list_sp_element_region = region_dataset.get_list_elements();

        auto list_index_elements_dataset = p_dataset->get_index_geometry_elements();
        list_global_cells_index.insert(list_global_cells_index.end(), list_index_elements_dataset.begin(),
                                       list_index_elements_dataset.end());

        std::vector<double> dataset_values = p_dataset->get_values();
        // std::transform(list_sp_element_region.begin(), list_sp_element_region.end(), dataset_values.begin(), [&](auto const pp_element) {
        //     return pp_element->get_scalar_data(function_name);
        // });
        list_values.insert(list_values.end(), dataset_values.begin(), dataset_values.end());
    }

    return {list_global_cells_index, list_values};
}

std::pair<std::vector<std::size_t>, std::vector<vector3>> mesh::get_vertices_index_value_of_vector_function(
    const std::string &function_name) const {
    std::vector<std::size_t> list_global_vertex_index{};
    std::vector<vector3>     list_values{};
    sp_vector_function       shared_pointer_function = get_sp_vector_function(function_name);
    for (auto &&p_dataset : shared_pointer_function->get_list_sp_datasets()) {
        auto index_region       = p_dataset->get_index_region_validity();
        auto region_dataset     = get_region_of_index(index_region);
        auto list_vertex_region = region_dataset.get_unique_vertices_as_vector();
        list_global_vertex_index.insert(list_global_vertex_index.end(), list_vertex_region.begin(), list_vertex_region.end());
        std::vector<vector3> dataset_values(list_vertex_region.size());
        std::transform(list_vertex_region.begin(), list_vertex_region.end(), dataset_values.begin(),
                       [&](auto const vtx_idx) { return get_vector_data_at_vertex(vtx_idx, function_name); });
        list_values.insert(list_values.end(), dataset_values.begin(), dataset_values.end());
    }
    return {list_global_vertex_index, list_values};
}

vector3 mesh::get_vector_data_at_vertex(std::size_t idx_vertex, const std::string &dataset_name) const {
    return m_ListVertices[idx_vertex].get_vector_data(dataset_name);
}

std::vector<vector3> mesh::get_all_vector_dataset_values(const std::string &dataset_name) const {
    std::vector<vector3> values(m_ListVertices.size());
    std::transform(m_ListVertices.begin(), m_ListVertices.end(), values.begin(),
                   [&](const vertex &vtx) { return vtx.get_vector_data(dataset_name); });
    return (values);
}

bool mesh::scalar_dataset_exists(const std::string &dataset_name) const {
    auto m_list_scalar_datasets = get_list_scalar_datasets();
    auto it_scalar_dataset      = std::find_if(m_list_scalar_datasets.begin(), m_list_scalar_datasets.end(), [&](const auto &sp_dataset) {
        if (sp_dataset == nullptr) {
            return false;
        }
        return sp_dataset->get_name() == dataset_name;
    });
    return it_scalar_dataset != m_list_scalar_datasets.end();
}

bool mesh::vector_dataset_exists(const std::string &dataset_name) const {
    auto m_list_vector_datasets = get_list_vector_datasets();
    auto it_scalar_dataset      = std::find_if(m_list_vector_datasets.begin(), m_list_vector_datasets.end(), [&](const auto &sp_dataset) {
        if (sp_dataset == nullptr) {
            return false;
        }
        return sp_dataset->get_name() == dataset_name;
    });
    return it_scalar_dataset != m_list_vector_datasets.end();
}

std::vector<int> mesh::get_vertices_belonging_to_other_regions(unsigned int idx_region) {
    /* First let's check if the region exists. */
    region *p_main_region = get_p_region(idx_region);
    if (p_main_region == nullptr) {
        throw std::runtime_error("Region does not exist.");
    }
    std::vector<int> MultipleRegionsVertex;
    auto             SetMainRegionUniqueVertices = p_main_region->get_unique_vertices();
    std::vector<int> MainRegionUniqueVertices(0);
    std::copy(SetMainRegionUniqueVertices.begin(), SetMainRegionUniqueVertices.end(), std::back_inserter(MainRegionUniqueVertices));
    for (auto &&bulk_region : m_ListRegionsBulk) {
        if (bulk_region.get_index() == idx_region) {
            continue;
        }
        auto             SetRegionUniqueVertices = bulk_region.get_unique_vertices();
        std::vector<int> RegionUniqueVertices(0);
        std::copy(SetRegionUniqueVertices.begin(), SetRegionUniqueVertices.end(), std::back_inserter(RegionUniqueVertices));
        std::vector<int> v_intersection;
        std::set_intersection(MainRegionUniqueVertices.begin(), MainRegionUniqueVertices.end(), RegionUniqueVertices.begin(),
                              RegionUniqueVertices.end(), std::back_inserter(v_intersection));

        std::copy(v_intersection.begin(), v_intersection.end(), std::back_inserter(MultipleRegionsVertex));
    }
    return MultipleRegionsVertex;
}

void mesh::re_index_datasets() {
    int new_dataset_index = -1;
    for (auto &scalar_dataset : get_list_scalar_datasets()) {
        scalar_dataset->set_index(++new_dataset_index);
    }
    for (auto &vector_dataset : get_list_vector_datasets()) {
        vector_dataset->set_index(++new_dataset_index);
    }
}

//  Functions methods

std::vector<std::string> mesh::get_scalar_functions_name() const {
    std::vector<std::string> list_scalar_functions;
    for (auto const &p_scalar_function : m_list_scalar_functions) {
        list_scalar_functions.push_back(p_scalar_function->get_name());
    }
    return list_scalar_functions;
}

std::vector<std::string> mesh::get_vector_functions_name() const {
    std::vector<std::string> list_vector_functions;
    for (auto const &p_vector_function : m_list_vector_functions) {
        list_vector_functions.push_back(p_vector_function->get_name());
    }
    return list_vector_functions;
}

std::vector<std::string> mesh::get_all_functions_names() const {
    std::set<std::string> set_functions_name;
    for (auto &&scalar_function : m_list_scalar_functions) {
        set_functions_name.insert(scalar_function->get_name());
    }
    for (auto &&scalar_function : m_list_vector_functions) {
        set_functions_name.insert(scalar_function->get_name());
    }
    std::vector<std::string> list_functions(set_functions_name.begin(), set_functions_name.end());
    return list_functions;
}

sp_scalar_function mesh::get_sp_scalar_function(const std::string &name) const {
    auto it_function = std::find_if(m_list_scalar_functions.begin(), m_list_scalar_functions.end(),
                                    [&](const auto &sp_func_scalar) { return (sp_func_scalar->get_name() == name); });
    if (it_function != m_list_scalar_functions.end()) {
        return *it_function;
    }
    return nullptr;
}

sp_vector_function mesh::get_sp_vector_function(const std::string &name) const {
    auto it_function = std::find_if(m_list_vector_functions.begin(), m_list_vector_functions.end(),
                                    [&](const auto &sp_func_vector) { return (sp_func_vector->get_name() == name); });
    if (it_function != m_list_vector_functions.end()) {
        return *it_function;
    }
    return nullptr;
}

void mesh::remove_scalar_function(const std::string &name) {
    LOG_INFO << "Removing the function : " << name;
    auto sp_function_to_remove = get_sp_scalar_function(name);
    sp_function_to_remove->remove_all_datasets();
    m_list_scalar_functions.erase(std::remove_if(m_list_scalar_functions.begin(), m_list_scalar_functions.end(),
                                                 [&](const auto &sp_scalar_func) { return (sp_scalar_func->get_name() == name); }),
                                  m_list_scalar_functions.end());
    re_index_datasets();
}

void mesh::remove_vector_function(const std::string &name) {
    LOG_INFO << "Removing the function : " << name;
    auto sp_function_to_remove = get_sp_vector_function(name);
    sp_function_to_remove->remove_all_datasets();
    m_list_vector_functions.erase(std::remove_if(m_list_vector_functions.begin(), m_list_vector_functions.end(),
                                                 [&](const auto &sp_vector_func) { return (sp_vector_func->get_name() == name); }),
                                  m_list_vector_functions.end());
    re_index_datasets();
}

void mesh::remove_all_scalar_function() {
    std::vector<std::string> list_func_name;
    for (auto &&scalar_func : m_list_scalar_functions) {
        list_func_name.push_back(scalar_func->get_name());
    }
    for (auto &&name : list_func_name) {
        remove_scalar_function(name);
    }
}

void mesh::remove_all_vector_function() {
    std::vector<std::string> list_func_name;
    for (auto &&vector_func : m_list_vector_functions) {
        list_func_name.push_back(vector_func->get_name());
    }
    for (auto &&name : list_func_name) {
        remove_vector_function(name);
    }
}

bool mesh::scalar_function_exists(const std::string &name) const { return (get_sp_scalar_function(name) != nullptr); }

bool mesh::vector_function_exists(const std::string &name) const { return (get_sp_vector_function(name) != nullptr); }

void mesh::add_scalar_function_to_vertices(const std::string &name) {
    auto my_function      = get_sp_scalar_function(name);
    auto list_sp_datasets = my_function->get_list_sp_datasets();
    for (const auto &my_sp_dataset : list_sp_datasets) {
        add_scalar_data_to_vertices(*my_sp_dataset);
    }
}

void mesh::creates_all_scalar_functions_from_datasets() {
    for (auto &scalar_dataset : get_list_scalar_datasets()) {
        std::string        dataset_name       = scalar_dataset->get_name();
        sp_scalar_function sp_ptr_function    = get_sp_scalar_function(dataset_name);
        DataLocationType   data_location_type = scalar_dataset->get_data_location_type();
        if (sp_ptr_function == nullptr) {
            sp_scalar_function my_new_sp_function = std::make_shared<function<double>>(dataset_name, DataType::scalar, data_location_type);
            my_new_sp_function->add_dataset(scalar_dataset);
            m_list_scalar_functions.push_back(my_new_sp_function);
        } else {
            sp_ptr_function->add_dataset(scalar_dataset);
        }
    }
}

void mesh::creates_all_vector_functions_from_datasets() {
    LOG_DEBUG << "CREATE FUNCTIONS ";
    for (auto &vector_dataset : get_list_vector_datasets()) {
        LOG_DEBUG << "CREATE FUNCTIONS ";
        std::string        dataset_name       = vector_dataset->get_name();
        sp_vector_function sp_ptr_function    = get_sp_vector_function(dataset_name);
        DataLocationType   data_location_type = vector_dataset->get_data_location_type();
        if (sp_ptr_function == nullptr) {
            sp_vector_function my_new_sp_function = std::make_shared<function<vector3>>(dataset_name, DataType::scalar, data_location_type);
            my_new_sp_function->add_dataset(vector_dataset);
            m_list_vector_functions.push_back(my_new_sp_function);
        } else {
            sp_ptr_function->add_dataset(vector_dataset);
        }
    }
}

void mesh::creates_all_functions_from_datasets() {
    creates_all_scalar_functions_from_datasets();
    creates_all_vector_functions_from_datasets();
    for (auto &&sp_sc_dt : m_list_scalar_functions) {
        LOG_INFO << "SCALAR FUNCTION NAME : " << sp_sc_dt->get_name();
    }
    for (auto &&sp_vt_dt : m_list_vector_functions) {
        LOG_INFO << "VECTOR FUNCTION NAME : " << sp_vt_dt->get_name();
    }
}

void mesh::create_scalar_function_from_function(const std::string &new_name, const std::string &origin_function_name) {
    sp_scalar_function sp_origin_function = get_sp_scalar_function(origin_function_name);
    if (sp_origin_function == nullptr) {
        LOG_ERROR << "The function " << origin_function_name << " does not exist";
        throw std::runtime_error("The function " + origin_function_name + " does not exist");
    }
    function           new_function    = sp_origin_function->get_function_copy(new_name);
    sp_scalar_function sp_new_function = std::make_shared<function<double>>(new_function);
    m_list_scalar_functions.push_back(sp_new_function);
    add_scalar_function_to_vertices(new_function.get_name());
    re_index_datasets();
}

void mesh::create_scalar_function_from_values_on_vertex(const std::string &function_name, const std::vector<double> &data_values) {
    // If the function already exists, it is first completely removed.
    if (get_sp_scalar_function(function_name) != nullptr) {
        LOG_WARNING << "Creating a new scalar function with a name" << function_name
                    << " already existing. Original function is overwriten.";
        remove_scalar_function(function_name);
    }
    DataLocationType   data_location_type = DataLocationType::vertex;
    sp_scalar_function my_new_function    = std::make_shared<function<double>>(function_name, DataType::scalar, data_location_type);
    for (auto &&p_region : m_ListRegionsBulk) {
        const int                index_region_validity = p_region.get_index();
        const int                index_dataset         = get_total_number_dataset() + 1;
        const DataType           Dataset_datatype      = DataType::scalar;
        const int                data_dimension        = 1;
        std::vector<double>      dataset_values{};
        std::vector<std::size_t> dataset_index_vtx{};
        for (auto &&vtx_index : p_region.get_unique_vertices()) {
            dataset_values.push_back(data_values[vtx_index]);
            dataset_index_vtx.push_back(vtx_index);
        }
        sp_scalar_dataset NewDataset =
            std::make_shared<dataset<double>>(function_name, index_dataset, index_region_validity, dataset_values, dataset_index_vtx,
                                              Dataset_datatype, data_location_type, data_dimension);
        add_scalar_data_to_vertices(*NewDataset);
        my_new_function->add_dataset(NewDataset);
    }
    m_list_scalar_functions.push_back(my_new_function);
    // re_index_datasets();
}

void mesh::create_vector_function_from_values_on_vertex(const std::string &function_name, const std::vector<vector3> &data_values) {
    auto ptr_function = get_sp_vector_function(function_name);
    if (ptr_function != nullptr) {
        remove_vector_function(function_name);
    }
    DataLocationType   data_location_type = DataLocationType::vertex;
    sp_vector_function my_new_function    = std::make_shared<function<vector3>>(function_name, DataType::vector, data_location_type);
    for (auto &&ptr_region : get_all_p_region()) {
        const int                index_region_validity = ptr_region->get_index();
        const int                index_dataset         = get_total_number_dataset() + 1;
        const DataType           Dataset_datatype      = DataType::vector;
        const int                data_dimension        = m_dimension;
        std::vector<vector3>     dataset_values{};
        std::vector<std::size_t> dataset_index_element{};
        for (auto &&vtx_index : ptr_region->get_unique_vertices()) {
            dataset_values.push_back(data_values[vtx_index]);
            dataset_index_element.push_back(vtx_index);
        }
        sp_vector_dataset NewDataset =
            std::make_shared<dataset<vector3>>(function_name, index_dataset, index_region_validity, dataset_values, dataset_index_element,
                                               Dataset_datatype, data_location_type, data_dimension);
        add_vector_data_to_vertices(*NewDataset);
        my_new_function->add_dataset(NewDataset);
    }
    m_list_vector_functions.push_back(my_new_function);
    re_index_datasets();
}

void mesh::create_scalar_function_from_values_on_element(const std::string &function_name, const std::vector<double> &values) {
    // If the function already exists, it is first completely removed.
    if (get_sp_scalar_function(function_name) != nullptr) {
        LOG_WARNING << "Creating a new scalar function with a name" << function_name
                    << " already existing. Original function is overwritten.";
        remove_scalar_function(function_name);
    }
    DataLocationType   data_location_type = DataLocationType::cell;
    sp_scalar_function my_new_function    = std::make_shared<function<double>>(function_name, DataType::scalar, data_location_type);
    for (auto &&p_region : m_ListRegionsBulk) {
        const int                index_region_validity = p_region.get_index();
        const int                index_dataset         = get_total_number_dataset() + 1;
        const DataType           Dataset_datatype      = DataType::scalar;
        const int                data_dimension        = 1;
        std::vector<double>      dataset_values{};
        std::vector<std::size_t> dataset_index_elements{};
        // for (auto &&element_idx : p_region.get_list_elements_index()) {
        //     dataset_values.push_back(values[element_idx]);
        //     dataset_index_elements.push_back(element_idx);
        // }
        std::size_t nb_elements = p_region.get_list_elements_index().size();
        for (std::size_t index_element = 0; index_element < nb_elements; ++index_element) {
            dataset_values.push_back(values[index_element]);
            dataset_index_elements.push_back(index_element);
        }
        sp_scalar_dataset NewDataset =
            std::make_shared<dataset<double>>(function_name, index_dataset, index_region_validity, dataset_values, dataset_index_elements,
                                              Dataset_datatype, data_location_type, data_dimension);
        add_scalar_data_to_elements(*NewDataset);
        my_new_function->add_dataset(NewDataset);
    }
    m_list_scalar_functions.push_back(my_new_function);
}

void mesh::create_vector_function_from_values_on_element(const std::string &name, const std::vector<vector3> &values) {
    // If the function already exists, it is first completely removed.
    if (get_sp_vector_function(name) != nullptr) {
        LOG_WARNING << "Creating a new vector function with a name" << name << " already existing. Original function is overwritten.";
        remove_vector_function(name);
    }
    DataLocationType   data_location_type = DataLocationType::cell;
    sp_vector_function my_new_function    = std::make_shared<function<vector3>>(name, DataType::vector, data_location_type);
    for (const auto &p_region : m_ListRegionsBulk) {
        const int      index_region_validity = p_region.get_index();
        const int      index_dataset         = get_total_number_dataset() + 1;
        const DataType Dataset_datatype      = DataType::vector;
        const int      data_dimension        = m_dimension;

        std::vector<vector3>     dataset_values{};
        std::vector<std::size_t> dataset_index_elements{};
        for (auto &&element_idx : p_region.get_list_elements_index()) {
            dataset_values.push_back(values[element_idx]);
            dataset_index_elements.push_back(element_idx);
        }
        sp_vector_dataset NewDataset =
            std::make_shared<dataset<vector3>>(name, index_dataset, index_region_validity, dataset_values, dataset_index_elements,
                                               Dataset_datatype, data_location_type, data_dimension);
        add_vector_data_to_elements(*NewDataset);
        my_new_function->add_dataset(NewDataset);
    }
    m_list_vector_functions.push_back(my_new_function);
    re_index_datasets();
}

void mesh::create_test_functions(const std::string &name_test_function) {
    std::vector<double>  scalar_values;
    std::vector<vector3> vector_values;
    const double         frequency = 13.0;
    const double         power     = 3.0;
    for (auto &&vtx : m_ListVertices) {
        double  value = pow(cos(frequency * vtx.norm()), power) - pow(sin(frequency * vtx.norm()), power);
        vector3 vector_value(vtx);
        scalar_values.push_back(value);
        vector_values.push_back(vector_value);
    }
    const std::string name_dataset_test        = name_test_function;
    const std::string name_dataset_vector_test = name_test_function + "_vector";
    create_scalar_function_from_values_on_vertex(name_dataset_test, scalar_values);
    create_vector_function_from_values_on_vertex(name_dataset_vector_test, vector_values);
    re_index_datasets();
}

void mesh::create_test_function_on_elements(const std::string &name_test_function) {
    auto                 list_bulk_elements = get_list_bulk_element();
    std::vector<double>  scalar_values;
    std::vector<vector3> vector_values;
    // const double         frequency = 13.0;
    // const double         power     = 3.0;
    // DEBUG
    std::size_t nb_elements = list_bulk_elements.size();
    std::cout << "Number of elements : " << nb_elements << std::endl;
    std::size_t counter = 0;
    for (std::size_t index_element = 0; index_element < nb_elements; ++index_element) {
        auto        p_element   = list_bulk_elements[index_element];
        std::size_t idx_element = p_element->get_index();
        if (idx_element != index_element + 1) {
            std::cout << "A : " << idx_element << std::endl;
            std::cout << "B : " << index_element << std::endl;
            ++counter;

            // throw std::runtime_error("Index element is not correct.");
        }
    }
    std::cout << "Nb of issues : " << counter << std::endl;
    std::cout << "End of test." << std::endl;

    for (auto &&p_element : list_bulk_elements) {
        double  value = p_element->get_barycenter().x();
        vector3 vector_value(p_element->get_barycenter());
        scalar_values.push_back(value);
        vector_values.push_back(vector_value);
    }
    const std::string name_dataset_test        = name_test_function;
    const std::string name_dataset_vector_test = name_test_function + "_vector";
    create_scalar_function_from_values_on_element(name_dataset_test, scalar_values);
    create_vector_function_from_values_on_element(name_dataset_vector_test, vector_values);
    re_index_datasets();
}

void mesh::create_null_scalar_function(const std::string &name) {
    std::vector<double> scalar_values(get_nb_vertices());
    create_scalar_function_from_values_on_vertex(name, scalar_values);
    re_index_datasets();
}

void mesh::create_null_vector_function(const std::string &name) {
    const vector3        null_vector{0.0, 0.0, 0.0};
    std::vector<vector3> vector_values(get_nb_vertices(), null_vector);
    create_vector_function_from_values_on_vertex(name, vector_values);
    re_index_datasets();
}

void mesh::create_gradient_function(const std::string &scalar_field, const std::string &new_name) {
    // fmt::print("Creating gradient function from scalar field {}\n", scalar_field);
    const vector3           null_vector{0.0, 0.0, 0.0};
    std::vector<vector3>    vector_values(get_nb_vertices(), null_vector);
    std::vector<double>     vector_average_renormalization(get_nb_vertices(), 0.0);
    std::vector<double>     vector_norm_values(get_nb_vertices(), 0.0);
    std::vector<sp_element> list_p_element = this->get_list_bulk_element();

    sp_scalar_function        scalar_func        = get_sp_scalar_function(scalar_field);
    std::vector<unsigned int> list_valid_regions = scalar_func->get_list_valid_regions_index();

    for (auto &&p_element : list_p_element) {
        if (std::find(list_valid_regions.begin(), list_valid_regions.end(), p_element->get_region_index()) == list_valid_regions.end()) {
            // std::cout << "Element " << p_element->get_index() << " is not in the valid region." << std::endl;
            continue;
        }
        std::vector<vertex *> p_vertices_list     = p_element->get_vertices();
        vector3               gradient_at_element = p_element->compute_gradient(scalar_field);
        // If the gradient is NaN, we set it to 0.
        double norm_gradient = gradient_at_element.norm();
        if (std::isnan(norm_gradient)) {
            gradient_at_element = null_vector;
        }
        for (std::size_t index_row = 0; index_row < p_vertices_list.size(); ++index_row) {
            vector_values[p_vertices_list[index_row]->get_index()] += gradient_at_element;
            vector_average_renormalization[p_vertices_list[index_row]->get_index()] += 1.0;
        }
    }
    for (std::size_t index_vtx = 0; index_vtx < get_nb_vertices(); ++index_vtx) {
        if (vector_average_renormalization[index_vtx] == 0.0) {
            vector_values[index_vtx] = null_vector;
            continue;
        }
        vector_values[index_vtx] *= (1.0 / vector_average_renormalization[index_vtx]);
        vector_norm_values[index_vtx] = vector_values[index_vtx].norm();
    }
    constexpr double micron_to_cm = -1e4;
    for (std::size_t index_vtx = 0; index_vtx < get_nb_vertices(); ++index_vtx) {
        vector_values[index_vtx] *= micron_to_cm;
        vector_norm_values[index_vtx] *= std::abs(micron_to_cm);
    }
    create_vector_function_from_values_on_vertex(new_name, vector_values);
    create_scalar_function_from_values_on_vertex(new_name + "_norm", vector_norm_values);
    re_index_datasets();
}

void mesh::create_scalar_function_from_grid_data(const std::string &new_name,
                                                 const grid_data   &my_grid_data,
                                                 const std::string &interpolation_method) {
    std::vector<double> scalar_values;
    scalar_values.reserve(get_nb_vertices());
    for (auto &&vtx : m_ListVertices) {
        if (interpolation_method == "linear") {
            scalar_values.push_back(my_grid_data.multi_linear_interpolate_scalar_at_location(vtx));
        } else if (interpolation_method == "nearest") {
            scalar_values.push_back(my_grid_data.nearest_neighbor_interpolate_scalar_at_location(vtx));
        } else {
            throw std::invalid_argument("Interpolation method not implemented.");
        }
    }
    std::cout << "Creating scalar function from grid data" << std::endl;
    create_scalar_function_from_values_on_vertex(new_name, scalar_values);
    re_index_datasets();
}

void mesh::create_scalar_function_from_csv_grid_data(const std::string &new_name,
                                                     const std::string &filename,
                                                     const std::string &interpolation_method) {
    grid_data my_grid_data{};
    my_grid_data.load_scalar_from_csv(filename, new_name);
    my_grid_data.print_grid_data_info();

    std::vector<double> scalar_values;
    scalar_values.reserve(get_nb_vertices());
    for (auto &&vtx : m_ListVertices) {
        if (interpolation_method == "linear") {
            scalar_values.push_back(my_grid_data.multi_linear_interpolate_scalar_at_location(vtx));
        } else if (interpolation_method == "nearest") {
            scalar_values.push_back(my_grid_data.nearest_neighbor_interpolate_scalar_at_location(vtx));
        } else {
            throw std::invalid_argument("Interpolation method not implemented.");
        }
    }
    create_scalar_function_from_values_on_vertex(new_name, scalar_values);
    re_index_datasets();
}

void mesh::create_density_function_from_list_positions_gaussian(const std::string          &new_fieldname,
                                                                const std::vector<vector3> &list_positions,
                                                                double                      gaussian_factor) {
    std::cout << "Creating density function from list positions gaussian" << std::endl;
    std::cout << "Number of positions: " << list_positions.size() << std::endl;
    std::cout << "Gaussian factor: " << gaussian_factor << std::endl;
    std::vector<double> scalar_values;
    scalar_values.reserve(get_nb_vertices());
    for (auto &&vtx : m_ListVertices) {
        double value = 0.0;
        for (const auto &pos : list_positions) {
            value += exp(-pow((vtx - pos).norm(), 2.0) * (1.0 / gaussian_factor));
        }
        scalar_values.push_back(value);
    }
    std::cout << "Creating scalar function from list positions gaussian: DONE." << std::endl;
    create_scalar_function_from_values_on_vertex(new_fieldname, scalar_values);
    re_index_datasets();
}

void mesh::addition_scalar_function_to_function(const std::string &name_function_to_add, const std::string &name_function_to_add_to) {
    auto sp_function_to_add    = get_sp_scalar_function(name_function_to_add);
    auto sp_function_to_add_to = get_sp_scalar_function(name_function_to_add_to);
    sp_function_to_add_to->add(*sp_function_to_add);
}

void mesh::multiplication_scalar_function_to_function(const std::string &name_new_function,
                                                      const std::string &name_function_1,
                                                      const std::string &name_function_2) {
    create_scalar_function_from_function(name_new_function, name_function_1);
    auto sp_function_new = get_sp_scalar_function(name_new_function);
    auto sp_function_2   = get_sp_scalar_function(name_function_2);
    sp_function_new->multiply(*sp_function_2);

    m_list_scalar_functions.push_back(sp_function_new);
    re_index_datasets();
}

void mesh::apply_scalar_function_to_function(const std::function<double(double)> &scalar_function,
                                             const std::string                   &name_function_to_apply_to) {
    auto sp_function_to_apply_to = get_sp_scalar_function(name_function_to_apply_to);
    sp_function_to_apply_to->apply_function(scalar_function);
}

// void mesh::apply_scalar_function_to_function(const std::string &name_new_function, const std::string &name_function_to_copy) {
//     auto sp_function_to_copy = get_sp_scalar_function(name_function_to_copy);
//     auto sp_function_new     = std::make_shared<scalar_function>(name_new_function, *sp_function_to_copy);
//     m_list_scalar_functions.push_back(sp_function_new);
//     re_index_datasets();
// }

/**
 * @brief Create a density function from a list of positions.
 *
 * The density is computed by counting, for each element, the number of positions that are inside the element.
 * The density is then normalized by the number of elements and the volume/surface of the element.
 *
 * @param new_fieldname
 * @param list_positions
 */
void mesh::create_density_function_from_list_positions_element_method(const std::string          &new_fieldname,
                                                                      const std::vector<vector3> &list_positions,
                                                                      double                      conversion_factor) {
    // std::cout << "Creating density function from list of positions" << std::endl;
    // std::cout << "Nb positions: " << list_positions.size() << std::endl;

    std::map<element *, int> map_element_nb_positions_inside;

    for (auto &&pos : list_positions) {
        vector3  position  = (m_dimension == 3) ? pos : pos.to_2d();
        element *p_element = find_element_at_location(position);
        if (p_element != nullptr) {
            if (map_element_nb_positions_inside.find(p_element) == map_element_nb_positions_inside.end()) {
                map_element_nb_positions_inside[p_element] = 0;
            }
            map_element_nb_positions_inside[p_element] += 1;
        }
    }

    std::vector<double> scalar_values(get_nb_vertices(), 0.0);
    std::vector<int>    count_per_vertices(get_nb_vertices(), 0);
    for (const auto &pair : map_element_nb_positions_inside) {
        std::vector<vertex *> p_vertices_list = pair.first->get_vertices();
        for (const auto &p_vtx : p_vertices_list) {
            scalar_values[p_vtx->get_index()] += pair.second / fabs(pair.first->get_measure());
            count_per_vertices[p_vtx->get_index()] += 1;
        }
    }
    // If a vertex "received" density from several elements, we average the density.
    for (std::size_t index_vtx = 0; index_vtx < get_nb_vertices(); ++index_vtx) {
        scalar_values[index_vtx] /= (count_per_vertices[index_vtx] > 0 ? count_per_vertices[index_vtx] : 1);
        scalar_values[index_vtx] *= conversion_factor;

        // DEBUG
        // if (scalar_values[index_vtx] > 2.0e19) {
        //     scalar_values[index_vtx] = 2.0e19;
        // }
        // if (std::abs(scalar_values[index_vtx] - 1.0e18) < 1e17 && m_ListVertices[index_vtx].x() > 0.5) {
        //     scalar_values[index_vtx] = 1.0e18;
        // }
    }
    create_scalar_function_from_values_on_vertex(new_fieldname, scalar_values);
    // const double integral_new_field = integrate_over_mesh(new_fieldname);
    // fmt::print("Integral of the field {} is {}\n", new_fieldname, integral_new_field);
}

// POISSON RELATED FUNCTIONS

void mesh::create_space_charge_function_from_vtx_values(const std::string &new_fieldname) {
    std::vector<double> space_charge_values(get_nb_vertices(), 0.0);
    for (auto &&vtx : m_ListVertices) {
        space_charge_values[vtx.get_index()] = vtx.get_space_charge();
    }
    create_scalar_function_from_values_on_vertex(new_fieldname, space_charge_values);
}

void mesh::compute_n_charge_on_elements(const std::vector<vector3> &list_positions) {
    auto list_bulk_elements = get_list_bulk_element();
    for (auto &&p_element : list_bulk_elements) {
        // std::cout << p_element->get_index() << std::endl;
        p_element->set_n_charge(0.0);
    }
    for (auto &&pos : list_positions) {
        vector3  position  = (m_dimension == 3) ? pos : pos.to_2d();
        element *p_element = find_element_at_location(position);
        // std::cout << "Electron X: " << p_element->get_barycenter().x() << std::endl;
        if (p_element != nullptr) {
            p_element->add_n_charge(1.0);
        }
    }
}

void mesh::compute_p_charge_on_elements(const std::vector<vector3> &list_positions) {
    auto list_bulk_elements = get_list_bulk_element();
    for (auto &&p_element : list_bulk_elements) {
        p_element->set_p_charge(0.0);
    }
    for (auto &&pos : list_positions) {
        vector3  position  = (m_dimension == 3) ? pos : pos.to_2d();
        element *p_element = find_element_at_location(position);
        if (p_element != nullptr) {
            p_element->add_p_charge(1.0);
        }
    }
}

void mesh::reset_all_space_charge_on_vertex() {
    for (auto &&vtx : m_ListVertices) {
        vtx.set_space_charge(0.0);
        vtx.set_charge_density(0.0);
    }
}

void mesh::recompute_space_charge() {
    for (auto &&vtx : m_ListVertices) {
        // space_charge = p - n + N_D - N_A
        double space_charge = vtx.get_charge_density() + vtx.get_doping_concentration();
        vtx.set_space_charge(space_charge);
    }
}

void mesh::update_charge_density_from_elements_values(int nb_iter_windows) {
    auto                list_bulk_elements = get_list_bulk_element();
    std::vector<double> n_density_values(list_bulk_elements.size());
    std::vector<double> p_density_values(list_bulk_elements.size());
    std::vector<double> doping_value(list_bulk_elements.size());
    std::vector<double> space_charge_element_values(list_bulk_elements.size());

    constexpr double conversion_factor = 1e12;

    for (std::size_t i = 0; i < list_bulk_elements.size(); ++i) {
        const auto &ptr_element        = list_bulk_elements[i];
        n_density_values[i]            = conversion_factor * (ptr_element->get_n_charge() / ptr_element->get_measure()) / nb_iter_windows;
        p_density_values[i]            = conversion_factor * (ptr_element->get_p_charge() / ptr_element->get_measure()) / nb_iter_windows;
        doping_value[i]                = ptr_element->integrate_scalar("DopingConcentration");
        space_charge_element_values[i] = p_density_values[i] - n_density_values[i] + doping_value[i];
    }
    const std::string n_density_field_name = "Poisson_nDensity";
    const std::string p_density_field_name = "Poisson_pDensity";
    create_scalar_function_from_values_on_element(n_density_field_name, n_density_values);
    create_scalar_function_from_values_on_element(p_density_field_name, p_density_values);
    create_scalar_function_from_values_on_element("DopingConcentration", doping_value);
    create_scalar_function_from_values_on_element("ADMCElementSpaceCharge", space_charge_element_values);
}

void mesh::update_charge_density_from_vertex_values() {
    std::vector<double> charge_density_values(m_ListVertices.size());
    std::vector<double> space_charge_values(m_ListVertices.size());
    std::transform(m_ListVertices.begin(), m_ListVertices.end(), charge_density_values.begin(),
                   [](const auto &vtx) { return vtx.get_charge_density(); });
    std::transform(m_ListVertices.begin(), m_ListVertices.end(), space_charge_values.begin(),
                   [](const auto &vtx) { return vtx.get_space_charge(); });

    create_scalar_function_from_values_on_vertex("MC_Density", charge_density_values);
    create_scalar_function_from_values_on_vertex("MC_SpaceCharge", space_charge_values);
}

void mesh::convert_element_function_to_vertex_function_scalar(const std::string &name_element_function,
                                                              const std::string &name_new_vertex_function) {
    auto                list_sp_bulk_elements = get_list_bulk_element();
    std::vector<double> scalar_values(get_nb_vertices(), 0.0);
    std::vector<double> volume_per_vertices(get_nb_vertices(), 0);
    for (auto &&p_element : list_sp_bulk_elements) {
        double                element_value   = p_element->get_scalar_data(name_element_function);
        std::vector<vertex *> p_vertices_list = p_element->get_vertices();
        double                volume_element  = std::fabs(p_element->get_measure());
        for (const auto &p_vtx : p_vertices_list) {
            scalar_values[p_vtx->get_index()] += element_value * volume_element;
            volume_per_vertices[p_vtx->get_index()] += volume_element;
        }
    }
    // If a vertex "received" density from several elements, we average the density.
    for (std::size_t index_vtx = 0; index_vtx < get_nb_vertices(); ++index_vtx) {
        scalar_values[index_vtx] /= (volume_per_vertices[index_vtx] > 0 ? volume_per_vertices[index_vtx] : 1);
    }
    create_scalar_function_from_values_on_vertex(name_new_vertex_function, scalar_values);
}

/**
 * @brief Convert the charge on elements into a charge on vertices.
 * The factor is used to potentially convert the charge from a unit to another OR to divide by the number of iterations
 * during which the charge was accumulated.
 *
 * @param factor
 */
void mesh::convert_charge_on_element_into_charge_at_vtx(double factor) {
    auto                list_sp_bulk_elements = get_list_bulk_element();
    std::vector<double> scalar_values(get_nb_vertices(), 0.0);
    std::vector<double> sum_volume_per_vertices(get_nb_vertices(), 0);
    constexpr double    conversion_factor = 1e12;

    for (auto &&p_element : list_sp_bulk_elements) {
        double                volume_element    = std::fabs(p_element->get_measure()) / conversion_factor;
        double                element_n_density = (p_element->get_n_charge() * factor) / volume_element;
        double                element_p_density = (p_element->get_p_charge() * factor) / volume_element;
        std::vector<vertex *> p_vertices_list   = p_element->get_vertices();
        for (const auto &p_vtx : p_vertices_list) {
            scalar_values[p_vtx->get_index()] += (element_p_density - element_n_density) * volume_element;
            sum_volume_per_vertices[p_vtx->get_index()] += volume_element;
        }
    }

    // If a vertex "received" density from several elements, we average the density.
    for (std::size_t index_vtx = 0; index_vtx < get_nb_vertices(); ++index_vtx) {
        scalar_values[index_vtx] /= (sum_volume_per_vertices[index_vtx] > 0 ? sum_volume_per_vertices[index_vtx] : 1);
        m_ListVertices[index_vtx].set_charge_density(scalar_values[index_vtx]);
    }
    create_scalar_function_from_values_on_vertex("MCC_Density", scalar_values);
    recompute_space_charge();
}

//  Interpolation

element *mesh::find_element_at_location(const vector3 &location) const { return m_p_search_tree->find_element_at_location(location); }

double mesh::interpolate_scalar_at_location(const std::string &fieldname, const vector3 &location) const {
    element *p_location_element = m_p_search_tree->find_element_at_location(location);
    if (p_location_element != nullptr) {
        return p_location_element->interpolate_scalar_at_location(fieldname, location);
    }
    LOG_WARNING << "No element found at position : " << location;
    return 0.0;
}

vector3 mesh::interpolate_vector_at_location(const std::string &fieldname, const vector3 &location) const {
    element *p_location_element = m_p_search_tree->find_element_at_location(location);
    if (p_location_element != nullptr) {
        return p_location_element->interpolate_vector_at_location(fieldname, location);
    }
    LOG_WARNING << "No element found at position : " << location;
    return vector3(0.0, 0.0, 0.0);
}

vector3 mesh::interpolate_gradient_at_location(const std::string &fieldname, const vector3 &location) const {
    element *p_location_element = m_p_search_tree->find_element_at_location(location);
    if (p_location_element != nullptr) {
        return p_location_element->compute_gradient(fieldname);
    }
    LOG_WARNING << "No element found at position : " << location;
    return vector3(0.0, 0.0, 0.0);
}

const region_bulk *mesh::get_p_region_at_location(const vector3 &location) {
    element *p_location_element = m_p_search_tree->find_element_at_location(location);
    if (p_location_element != nullptr) {
        int region_index = p_location_element->get_region_index();
        return get_p_region_bulk(region_index);
    }
    LOG_WARNING << "No region found at position : " << location;
    return nullptr;
}

std::string mesh::get_region_name_at_location(const vector3 &location) {
    const region_bulk *region_bulk_ptr = get_p_region_at_location(location);
    return region_bulk_ptr != nullptr ? region_bulk_ptr->get_name() : "";
}

std::string mesh::get_material_name_at_location(const vector3 &location) {
    const region_bulk *region_bulk_ptr = get_p_region_at_location(location);
    return region_bulk_ptr != nullptr ? region_bulk_ptr->get_material() : "";
}

std::string mesh::get_material_name_at_element(element *p_element_location) const {
    if (p_element_location != nullptr) {
        int region_index = p_element_location->get_region_index();
        return get_p_region_bulk(region_index)->get_material();
    }
    return "";
}

std::pair<vector3, double> mesh::get_argmax_max_of_function(const std::string &fieldname) const {
    vector3 position_max{};
    double  max_function = -std::numeric_limits<double>::infinity();
    for (const auto &vtx : m_ListVertices) {
        double value_at_vtx = vtx.get_scalar_data(fieldname);
        if (value_at_vtx > max_function) {
            position_max = vtx;
            max_function = value_at_vtx;
        }
    }
    return {position_max, max_function};
}

std::pair<vector3, double> mesh::get_argmin_min_of_function(const std::string &fieldname) const {
    vector3 position_min{};
    double  min_function = std::numeric_limits<double>::max();
    for (const auto &vtx : m_ListVertices) {
        double value_at_vtx = vtx.get_scalar_data(fieldname);
        if (value_at_vtx > min_function) {
            position_min = vtx;
            min_function = value_at_vtx;
        }
    }
    return {position_min, min_function};
}

// Find elements

std::vector<element *> mesh::find_elements_overlapping_box(const bbox &my_box) const {
    return m_p_search_tree->find_elements_overlapping_box(my_box);
}

// Intersection

/**
 * @brief Find the intersection location between a line and the mesh boundaries.
 *
 * @param point_origin
 * @param point_end
 * @return std::optional<vector3>
 */
std::optional<vector3> mesh::find_location_line_boundaries_intersection(const vector3 &point_origin, const vector3 &point_end) const {
    auto *elemnt_ptA = find_element_at_location(point_origin);
    auto *elemnt_ptB = find_element_at_location(point_end);
    if (elemnt_ptA == nullptr) {
        LOG_ERROR << point_origin;
        return std::nullopt;
        // throw std::runtime_error("Error: origin point already out of the device mesh.");
    }
    if (elemnt_ptB != nullptr) {
        return std::nullopt;
    }

    const std::size_t max_iter          = 1000;
    constexpr double  epsilon           = 1.0e-15;
    vector3           point_A           = point_origin;
    vector3           point_B           = point_end;
    vector3           point_middle      = middle(point_A, point_B);
    auto             *elemnt_pt_middle  = find_element_at_location(point_middle);
    double            distance_A_middle = distance(point_A, point_middle);
    std::size_t       iter              = 0;
    while ((++iter < max_iter) && (distance_A_middle > epsilon || elemnt_pt_middle == nullptr)) {
        elemnt_pt_middle = find_element_at_location(point_middle);
        if (elemnt_pt_middle == nullptr) {
            // The middle point is outside of the device
            point_B = point_middle;
        } else {
            point_A = point_middle;
        }
        point_middle      = middle(point_A, point_B);
        distance_A_middle = distance(point_A, point_middle);
    }
    if ((distance_A_middle <= epsilon && elemnt_pt_middle != nullptr)) {
        return point_A;
    } else {
        return std::nullopt;
    }
}

/**
 * @brief Find the intersection location between a line and the boundary of the region in which the origin point belongs.
 *
 * @param point_origin
 * @param point_end
 * @return std::optional<vector3>
 */
std::optional<vector3> mesh::find_location_line_region_interface_intersection(const vector3 &point_origin, const vector3 &point_end) const {
    auto *elemnt_ptA = find_element_at_location(point_origin);
    auto *elemnt_ptB = find_element_at_location(point_end);
    if (elemnt_ptA == nullptr) {
        return std::nullopt;
        // throw std::runtime_error("Error: origin point already out of the device mesh.");
    } else if (elemnt_ptB == nullptr) {
        return std::nullopt;
    }
    const int index_region_start_point = elemnt_ptA->get_region_index();
    const int index_region_end_point   = elemnt_ptB->get_region_index();

    if (index_region_start_point == index_region_end_point) {
        // std::cout << "Both points of the lines are aleady in the same region, cannot find an intersection with an interface.\n";
        return std::nullopt;
    }

    const std::size_t max_iter                  = 1000;
    constexpr double  epsilon                   = 1.0e-15;
    vector3           point_A                   = point_origin;
    vector3           point_B                   = point_end;
    vector3           point_middle              = middle(point_A, point_B);
    auto             *elemnt_pt_middle          = find_element_at_location(point_middle);
    int               index_region_middle_point = elemnt_pt_middle->get_region_index();
    double            distance_A_middle         = distance(point_A, point_middle);
    std::size_t       iter                      = 0;
    while ((++iter < max_iter) && (distance_A_middle > epsilon || index_region_middle_point != index_region_start_point)) {
        elemnt_pt_middle = find_element_at_location(point_middle);
        // index_region_middle_point = elemnt_pt_middle->get_region_index();
        if (elemnt_pt_middle == nullptr || elemnt_pt_middle->get_region_index() != index_region_start_point) {
            // The middle point is outside of the device
            point_B = point_middle;
        } else {
            point_A = point_middle;
        }
        point_middle      = middle(point_A, point_B);
        distance_A_middle = distance(point_A, point_middle);
    }
    if ((distance_A_middle <= epsilon && elemnt_pt_middle != nullptr && index_region_middle_point == index_region_start_point)) {
        return point_A;
    }
    return std::nullopt;
}

/**
 * @brief Find the element, and the location of the intersection between a line and the mesh boundaries.
 *
 * @param point_A
 * @param point_B
 * @return std::optional<std::pair<sp_element, vector3>>
 */
std::optional<std::pair<sp_element, vector3>> mesh::find_element_face_line_boundaries_intersection(const vector3 &point_A,
                                                                                                   const vector3 &point_B) const {
    auto location_intersection = find_location_line_boundaries_intersection(point_A, point_B);
    if (!location_intersection.has_value()) {
        return std::nullopt;
    }
    auto                         *elemnt_intersection          = find_element_at_location(location_intersection.value());
    std::map<sp_element, vector3> faces_positions_intersection = elemnt_intersection->compute_element_line_intersection(point_A, point_B);
    if (faces_positions_intersection.empty()) {
        return std::nullopt;
    }
    if (faces_positions_intersection.size() == 1) {
        return *(faces_positions_intersection.begin());
    }
    double                         distance_to_intersection             = 1e100;
    auto                           it_smallest_distance_to_intersection = faces_positions_intersection.begin();
    std::pair<sp_element, vector3> result_intersection;
    for (auto it = faces_positions_intersection.begin(); it != faces_positions_intersection.end(); ++it) {
        double new_distance = distance(it->second, location_intersection.value());
        if (new_distance < distance_to_intersection) {
            it_smallest_distance_to_intersection = it;
            result_intersection                  = *it;
            distance_to_intersection             = new_distance;
        }
    }
    return result_intersection;
}

std::optional<std::pair<sp_element, vector3>> mesh::find_element_face_line_region_interface_intersection(const vector3 &point_A,
                                                                                                         const vector3 &point_B) const {
    auto location_intersection = find_location_line_region_interface_intersection(point_A, point_B);
    if (!location_intersection.has_value()) {
        return std::nullopt;
    }
    auto                         *elemnt_intersection          = find_element_at_location(location_intersection.value());
    std::map<sp_element, vector3> faces_positions_intersection = elemnt_intersection->compute_element_line_intersection(point_A, point_B);
    if (faces_positions_intersection.empty()) {
        return std::nullopt;
    }
    if (faces_positions_intersection.size() == 1) {
        return *(faces_positions_intersection.begin());
    }
    double                         distance_to_intersection = 1e100;
    std::pair<sp_element, vector3> result_intersection;
    for (auto it = faces_positions_intersection.begin(); it != faces_positions_intersection.end(); ++it) {
        double new_distance = distance(it->second, location_intersection.value());
        if (new_distance < distance_to_intersection) {
            result_intersection      = *it;
            distance_to_intersection = new_distance;
        }
    }
    return result_intersection;
}

std::optional<std::pair<sp_element, vector3>> mesh::find_line_first_intersection(const vector3 &point_A, const vector3 &point_B) const {
    assert(0);
    if (point_A == point_B) {
        return std::nullopt;
    }
    return std::nullopt;
}

//  Integration

double mesh::integrate_over_mesh(const std::string &fieldname) const {
    double integral{0.0};
    for (auto &&p_bulk_element : get_list_p_bulk_element()) {
        // integral += fabs(p_bulk_element->get_measure()) *
        //             p_bulk_element->interpolate_scalar_at_location(fieldname, p_bulk_element->get_barycenter());
        double to_add = fabs(p_bulk_element->get_measure()) *
                        p_bulk_element->interpolate_scalar_at_location(fieldname, p_bulk_element->get_barycenter());
        integral += to_add;
    }
    return integral;
}

double mesh::integrate_over_mesh_element_data(const std::string &fieldname) const {
    double integral{0.0};
    for (auto &&p_bulk_element : get_list_p_bulk_element()) {
        integral += fabs(p_bulk_element->get_measure()) * p_bulk_element->get_scalar_data(fieldname);
    }
    return integral;
}

double mesh::integrate_over_region(const std::string &fieldname, const std::string &region_name) const {
    const auto *p_region = this->get_p_region(region_name);
    if (p_region == nullptr) {
        LOG_ERROR << "ERROR : THE REGION IS UKNOWN, IMPOSSIBLE TO COMPUTE THE INTEGRAL, WILL RETURN 0. : " << region_name;
        return 0.0;
    }
    double integral{0.0};
    for (auto &&p_bulk_element : p_region->get_list_elements()) {
        integral += fabs(p_bulk_element->get_measure()) *
                    p_bulk_element->interpolate_scalar_at_location(fieldname, p_bulk_element->get_barycenter());
    }
    return integral;
}

double mesh::integrate_over_material(const std::string &fieldname, const std::string &material_name) const {
    double integral{0.0};
    for (auto &&p_bulk_element : get_list_p_bulk_element()) {
        if (get_material_name_at_element(p_bulk_element) == material_name) {
            integral += fabs(p_bulk_element->get_measure()) *
                        p_bulk_element->interpolate_scalar_at_location(fieldname, p_bulk_element->get_barycenter());
        }
    }
    return integral;
}

double mesh::compute_total_mesh_volume() const {
    double total_volume{0.0};
    for (auto &&p_bulk_element : get_list_p_bulk_element()) {
        total_volume += fabs(p_bulk_element->get_measure());
    }
    return total_volume;
}

double mesh::compute_region_volume(const std::string &region_name) const {
    const auto *p_region = this->get_p_region(region_name);
    if (p_region == nullptr) {
        LOG_ERROR << "ERROR : THE REGION IS UKNOWN, IMPOSSIBLE TO COMPUTE THE INTEGRAL, WILL RETURN 0. : " << region_name;
        return 0.0;
    }
    double total_volume{0.0};
    for (auto &&p_bulk_element : p_region->get_list_elements()) {
        total_volume += fabs(p_bulk_element->get_measure());
    }
    return total_volume;
}

double mesh::compute_material_volume(const std::string &material_name) const {
    double total_volume{0.0};
    for (auto &&p_bulk_element : get_list_p_bulk_element()) {
        if (get_material_name_at_element(p_bulk_element) == material_name) {
            total_volume += fabs(p_bulk_element->get_measure());
        }
    }
    return total_volume;
}

// IO

void mesh::print_regions_info() const {
    std::cout << "BULK REGIONS : " << std::endl;
    for (const auto &my_region : m_ListRegionsBulk) {
        my_region.print_info();
    }
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "INTERFACE REGIONS : " << std::endl;
    for (const auto &my_region : m_ListRegionsInterface) {
        my_region.print_info();
    }
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "CONTACT REGIONS : " << std::endl;
    for (const auto &my_region : m_ListRegionsContact) {
        my_region.print_info();
    }
    std::cout << "---------------------------------------------" << std::endl;
}

void mesh::print_datasets_info() const {
    std::cout << "DATASET INFO " << std::endl;
    std::cout << "SCALAR DATASET :  " << std::endl;
    for (auto &&sp_dataset : get_list_scalar_datasets()) {
        std::cout << "DATASET NAME              : " << sp_dataset->get_name() << std::endl;
        std::cout << "DATASET REGION VALIDITY   : " << sp_dataset->get_index_region_validity() << std::endl;
        std::cout << "DATASET INDEX             : " << sp_dataset->get_index() << std::endl;
        std::cout << "--------------" << std::endl;
    }
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "SCALAR DATASET :  " << std::endl;
    for (auto &&vector_dtset : get_list_vector_datasets()) {
        std::cout << "DATASET NAME              : " << vector_dtset->get_name() << std::endl;
        std::cout << "DATASET REGION VALIDITY   : " << vector_dtset->get_index_region_validity() << std::endl;
        std::cout << "DATASET INDEX             : " << vector_dtset->get_index() << std::endl;
        std::cout << "--------------" << std::endl;
    }
    std::cout << "---------------------------------------------" << std::endl;
}

void mesh::print_functions_info() const {
    std::cout << "FUNCTIONS INFO " << std::endl;
    for (const auto &function_scalar : m_list_scalar_functions) {
        std::cout << "Function : " << function_scalar->get_name() << " --->  " << function_scalar->get_datatype() << std::endl;
    }
    for (const auto &function_vector : m_list_vector_functions) {
        std::cout << "Function : " << function_vector->get_name() << " --->  " << function_vector->get_datatype() << std::endl;
    }
}

void mesh::export_vertices_to_csv(const std::string &filename) const {
    std::ofstream file_csv(filename);
    file_csv << "X,Y,Z\n";
    for (const auto &vtx : m_ListVertices) {
        file_csv << vtx.x() << ',';
        file_csv << vtx.y() << ',';
        file_csv << vtx.z() << "\n";
    }
    file_csv.close();
}

void mesh::export_vertices_data_to_csv(const std::string &filename, const std::string &dataset_name) const {
    std::ofstream file_csv(filename);
    file_csv << "X,Y,Z," << dataset_name << std::endl;
    for (const auto &vtx : m_ListVertices) {
        file_csv << vtx.x() << ',';
        file_csv << vtx.y() << ',';
        file_csv << vtx.z() << ',';
        file_csv << vtx.get_scalar_data(dataset_name) << "\n";
    }
    file_csv.close();
}

void mesh::export_all_vertices_data_to_csv(const std::string &filename) const {
    std::ofstream file_csv(filename);
    file_csv << "X,Y,Z,";
    auto list_function = get_all_functions_names();
    for (auto &&dataset_name : list_function) {
        file_csv << dataset_name << ",";
    }
    file_csv << "BlankColumn" << std::endl;
    for (const auto &bulk_region : get_all_p_bulk_region()) {
        for (const auto &vtx_index : bulk_region->get_unique_vertices_as_vector()) {
            auto vtx = m_ListVertices[vtx_index];
            file_csv << vtx.x() << ',';
            file_csv << vtx.y() << ',';
            file_csv << vtx.z() << ',';
            for (auto &&dataset_name : list_function) {
                file_csv << vtx.get_scalar_data(dataset_name) << ",";
            }
            file_csv << "0\n";
        }
    }
    file_csv.close();
}

std::vector<vector3> mesh::generate_mesh_grid(double dx, double dy, double dz) const {
    bbox                 my_box = get_bounding_box();
    std::size_t          N_x    = dx == 0.0 ? 0 : abs(int((my_box.get_x_max() - my_box.get_x_min()) / dx));
    std::size_t          N_y    = dy == 0.0 ? 0 : abs(int((my_box.get_y_max() - my_box.get_y_min()) / dy));
    std::vector<vector3> mesh_grid;
    if (m_dimension == 2) {
        std::cout << "Nx : " << int(my_box.get_x_max() - my_box.get_x_min()) << "\n";
        std::cout << "Nx : " << int(my_box.get_x_max() - my_box.get_x_min()) << "\n";
        std::cout << "Nx : " << N_x << "\n";
        std::cout << "Ny : " << N_y << "\n";
        mesh_grid = my_box.generate_mesh_grid_2d(N_x, N_y);
    } else if (m_dimension == 3) {
        std::size_t N_z = dz == 0.0 ? 0 : int((my_box.get_z_max() - my_box.get_z_min()) / dz);
        mesh_grid       = my_box.generate_mesh_grid_3d(N_x, N_y, N_z);
    }
    return mesh_grid;
}

std::vector<vector3> mesh::generate_mesh_grid(int N_x, int N_y, int N_z) const {
    bbox                 my_box = get_bounding_box();
    std::vector<vector3> mesh_grid;
    if (m_dimension == 2) {
        mesh_grid = my_box.generate_mesh_grid_2d(N_x, N_y);
    } else if (m_dimension == 3) {
        mesh_grid = my_box.generate_mesh_grid_3d(N_x, N_y, N_z);
    }
    return mesh_grid;
}

// void mesh::export_on_grid(const std::string &filename, double dx, double dy, double dz) {
//     auto          list_elements_bulk = get_list_p_bulk_element();
//     std::ofstream file_grid(filename);
//     file_grid << "X,Y,Z,";
//     auto list_function = get_all_functions_names();
//     for (auto &&dataset_name : list_function) {
//         file_grid << dataset_name << ",";
//     }
//     file_grid << "BlankColumn" << std::endl;
//     for (const auto &bulk_element : list_elements_bulk) {
//         auto point = bulk_element->get_barycenter();
//         file_grid << point.x() << ',';
//         file_grid << point.y() << ',';
//         file_grid << point.z() << ',';
//         for (auto &&field_name : list_function) {
//             file_grid << interpolate_scalar_at_location(field_name, point) << ",";
//         }
//         file_grid << "0\n";
//     }
//     file_grid.close();
// }

// void mesh::export_on_grid(const std::string &filename, double dx, double dy, double dz) {
//     std::vector<vector3> grid_mesh = generate_mesh_grid(dx, dy, dz);

//     std::ofstream file_grid(filename);
//     file_grid << "X,Y,Z,";
//     auto list_function = get_all_functions_names();
//     for (auto &&dataset_name : list_function) {
//         file_grid << dataset_name << ",";
//     }
//     file_grid << "BlankColumn" << std::endl;
//     for (const auto &point : grid_mesh) {
//         file_grid << point.x() << ',';
//         file_grid << point.y() << ',';
//         file_grid << point.z() << ',';
//         for (auto &&field_name : list_function) {
//             file_grid << interpolate_scalar_at_location(field_name, point) << ",";
//         }
//         file_grid << "0\n";
//     }
//     file_grid.close();
// }

void mesh::export_on_grid(const std::string &filename, double dx, double dy, double dz) const {
    std::vector<vector3> grid_mesh = generate_mesh_grid(dx, dy, dz);

    std::ofstream file_grid(filename);
    file_grid << "X,Y,Z,";
    auto list_function = get_all_functions_names();
    for (auto &&dataset_name : list_function) {
        file_grid << dataset_name << ",";
    }
    file_grid << "BlankColumn" << std::endl;
    for (const auto &point : grid_mesh) {
        file_grid << point.x() << ',';
        file_grid << point.y() << ',';
        file_grid << point.z() << ',';
        for (auto &&field_name : list_function) {
            file_grid << interpolate_scalar_at_location(field_name, point) << ",";
        }
        file_grid << "0\n";
    }
    file_grid.close();
}

void mesh::export_z_cut(const std::string &filename, double x_const, double y_const, double dz) const {
    // epsilon is used to avoid to be on the boundary of the mesh (1e-6mum=0.001nm)
    constexpr double     epsilon = 1.0e-6;
    const double         z_min   = get_bounding_box().get_z_min() + epsilon;
    const double         z_max   = get_bounding_box().get_z_max() - epsilon;
    std::size_t          Nz      = static_cast<int>((z_max - z_min) / dz);
    bbox                 z_line_box(x_const, x_const, y_const, y_const, z_min, z_max);
    std::vector<vector3> z_line_grid = z_line_box.generate_mesh_grid_3d(1, 1, Nz);
    std::ofstream        file_export(filename);
    file_export << "Z,";
    auto list_function = get_all_functions_names();
    for (auto &&dataset_name : list_function) {
        file_export << dataset_name << ",";
    }
    file_export << "BlankColumn" << std::endl;
    for (const auto &point : z_line_grid) {
        file_export << point.z() << ',';
        for (auto &&field_name : list_function) {
            file_export << interpolate_scalar_at_location(field_name, point) << ",";
        }
        file_export << "0\n";
    }
    file_export.close();
}

void mesh::export_y_cut(const std::string &filename, double x_const, double z_const, double dy) const {
    // epsilon is used to avoid to be on the boundary of the mesh (1e-6mum=0.001nm)
    constexpr double     epsilon = 1.0e-6;
    const double         y_min   = get_bounding_box().get_y_min() + epsilon;
    const double         y_max   = get_bounding_box().get_y_max() - epsilon;
    std::size_t          Ny      = static_cast<int>((y_max - y_min) / dy);
    bbox                 y_line_box(x_const, x_const, y_min, y_max, z_const, z_const);
    std::vector<vector3> y_line_grid = y_line_box.generate_mesh_grid_3d(1, Ny, 1);
    std::ofstream        file_export(filename);
    file_export << "Y,";
    auto list_function = get_all_functions_names();
    for (auto &&dataset_name : list_function) {
        file_export << dataset_name << ",";
    }
    file_export << "BlankColumn" << std::endl;
    for (const auto &point : y_line_grid) {
        file_export << point.y() << ',';
        for (auto &&field_name : list_function) {
            file_export << interpolate_scalar_at_location(field_name, point) << ",";
        }
        file_export << "0\n";
    }
    file_export.close();
}

void mesh::export_x_cut(const std::string &filename, double y_const, double z_const, double dx) const {
    // epsilon is used to avoid to be on the boundary of the mesh (1e-6mum=0.001nm)
    constexpr double     epsilon = 1.0e-6;
    const double         x_min   = get_bounding_box().get_x_min() + epsilon;
    const double         x_max   = get_bounding_box().get_x_max() - epsilon;
    std::size_t          Nx      = static_cast<int>((x_max - x_min) / dx);
    bbox                 x_line_box(x_min, x_max, y_const, y_const, z_const, z_const);
    std::vector<vector3> x_line_grid = x_line_box.generate_mesh_grid_3d(Nx, 1, 1);
    std::ofstream        file_export(filename);
    file_export << "X,";
    auto list_function = get_all_functions_names();
    for (auto &&dataset_name : list_function) {
        file_export << dataset_name << ",";
    }
    file_export << "BlankColumn" << std::endl;
    for (const auto &point : x_line_grid) {
        file_export << point.x() << ',';
        for (auto &&field_name : list_function) {
            file_export << interpolate_scalar_at_location(field_name, point) << ",";
        }
        file_export << "0\n";
    }
    file_export.close();
}

// Others

/**
 * @brief  This functions purpose is to change the regions of the mesh in a certain way.
 *
 * All the element on which the ElementData mundfab_data_name is equal to value_silicon will be transfer to the Silicon region.
 *
 * THIS FUNCTION IS FOR TEST ONLY.
 *
 * @param mundfab_data_name
 * @param value_silicon
 */
void mesh::mundfabisation(const std::string &mundfab_data_name, const double value_silicon, const std::string &silicon_region_name) {
    region *p_silicon_region = get_p_region(silicon_region_name);
    for (region &reg_bulk : m_ListRegionsBulk) {
        const std::string region_name = reg_bulk.get_name();
        LOG_DEBUG << "TRANSFERING ELEMENT FROM REGION : " << region_name;
        if (region_name == silicon_region_name) {
            continue;
        }
        std::vector<std::size_t> list_element_to_transfer;
        auto                     list_elements = reg_bulk.get_list_elements();
        for (auto &&reg_element : list_elements) {
            double value_mundfab = reg_element->get_scalar_data(mundfab_data_name);
            if (value_mundfab == value_silicon) {
                list_element_to_transfer.push_back(reg_element->get_index());
            }
        }
        LOG_DEBUG << "NUMBER OF ELEMENTS TO BE TRANSFERED : " << list_element_to_transfer.size();
        this->transfer_elements_to_other_region(list_element_to_transfer, &reg_bulk, p_silicon_region);

        auto *p_region_gas       = get_p_region("Gas_1");
        auto *p_region_silicon   = get_p_region("Silicon_1");
        auto *p_interface_region = get_p_interface_region_from_bulk_indices(p_region_gas->get_index(), p_region_silicon->get_index());
        remove_region(p_interface_region->get_index());
        compute_interface_region_betwwen_two_bulks(p_region_gas->get_index(), p_region_silicon->get_index(), list_element_to_transfer);
    }
}

}  // namespace mesh


} // namespace uepm