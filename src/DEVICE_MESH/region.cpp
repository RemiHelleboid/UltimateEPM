/**
 * @file region.cpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief
 * @version 0.1
 * @date 2021-07-10
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "region.hpp"

#include <algorithm>
#include <cassert>

namespace uepm {

namespace mesh {

/**
 * @brief Overloading of operator<< for RegionType enum class.
 *
 * @param os
 * @param type
 * @return std::ostream&
 */
std::ostream &operator<<(std::ostream &os, const RegionType &type) {
    if (type == RegionType::bulk) {
        os << "bulk";
    } else if (type == RegionType::interface) {
        os << "interface";
    } else if (type == RegionType::contact) {
        os << "contact";
    } else if (type == RegionType::unknown) {
        os << "unknown" << std::endl;
    } else {
        assert("Error: operator<< overloaded for RegionType does not support the type.");
    }
    return os;
}

void region::add_element(const std::shared_ptr<element> &sp_elem) {
    sp_elem->set_region_index(m_region_index);
    m_ListElements.insert({sp_elem->get_index(), sp_elem});
}

bool region::is_element_in_region(std::size_t element_index) { return m_ListElements.find(element_index) != m_ListElements.end(); }

std::shared_ptr<element> region::get_p_element(std::size_t index_element) { return m_ListElements.at(index_element); }

void region::remove_element(std::size_t element_index) { m_ListElements.erase(element_index); }

void region::remove_elements(const std::vector<std::size_t> &list_element_index) {
    for (auto &&index_element : list_element_index) {
        m_ListElements.erase(index_element);
    }
}

std::vector<std::shared_ptr<element>> region::get_list_elements() const {
    std::vector<std::shared_ptr<element>> list_sp_elements;
    list_sp_elements.reserve(m_ListElements.size());
    std::transform(m_ListElements.begin(), m_ListElements.end(), std::back_inserter(list_sp_elements), [](const auto &maps_element) {
        return maps_element.second;
    });
    return list_sp_elements;
}

std::vector<std::size_t> region::get_list_elements_index() const {
    std::vector<std::size_t> hey_elements_idx_list{};
    hey_elements_idx_list.reserve(m_ListElements.size());
    std::transform(m_ListElements.begin(), m_ListElements.end(), std::back_inserter(hey_elements_idx_list), [](const auto &maps_element) {
        return maps_element.first;
    });
    return hey_elements_idx_list;
}

void region::compute_unique_vertices() {
    m_UniqueVertices.clear();
    for (const auto &elem : m_ListElements) {
        auto list_vertices_index = (elem.second)->get_vertices_index();
        for (auto idx_vtx : list_vertices_index) {
            m_UniqueVertices.insert(idx_vtx);
        }
    }
}

std::vector<std::size_t> region::get_unique_vertices_as_vector() const {
    std::vector<std::size_t> v(m_UniqueVertices.begin(), m_UniqueVertices.end());
    return v;
}

std::vector<vertex *> region::get_list_all_p_vertices() const {
    std::vector<vertex *> list_vertices;
    for (const auto &elem : m_ListElements) {
        auto list_element_vertices = (elem.second)->get_vertices();
        std::copy(list_element_vertices.begin(), list_element_vertices.end(), std::back_inserter(list_vertices));
    }
    return list_vertices;
}

bool region::is_vertex_inside(std::size_t index_vertex) const {
    return (std::find(m_UniqueVertices.begin(), m_UniqueVertices.end(), index_vertex) != m_UniqueVertices.end());
}

bbox region::compute_bounding_box() const {
    auto                list_all_vertices = get_list_all_p_vertices();
    std::vector<double> X_coords(list_all_vertices.size());
    std::vector<double> Y_coords(list_all_vertices.size());
    std::vector<double> Z_coords(list_all_vertices.size());
    std::transform(list_all_vertices.begin(), list_all_vertices.end(), X_coords.begin(), [&](const auto &p_vtx) { return p_vtx->x(); });
    std::transform(list_all_vertices.begin(), list_all_vertices.end(), Y_coords.begin(), [&](const auto &p_vtx) { return p_vtx->y(); });
    std::transform(list_all_vertices.begin(), list_all_vertices.end(), Z_coords.begin(), [&](const auto &p_vtx) { return p_vtx->z(); });
    const double x_min = *std::min_element(X_coords.begin(), X_coords.end());
    const double x_max = *std::max_element(X_coords.begin(), X_coords.end());
    const double y_min = *std::min_element(Y_coords.begin(), Y_coords.end());
    const double y_max = *std::max_element(Y_coords.begin(), Y_coords.end());
    const double z_min = *std::min_element(Z_coords.begin(), Z_coords.end());
    const double z_max = *std::max_element(Z_coords.begin(), Z_coords.end());
    return bbox(x_min, x_max, y_min, y_max, z_min, z_max);
}

void region::print_info() const {
    std::cout << "Region name             : " << get_name() << std::endl;
    std::cout << "Region index            : " << get_index() << std::endl;
    std::cout << "Region type             : " << get_region_type() << std::endl;
    std::cout << "Number elements         : " << get_number_elements() << std::endl;
    std::cout << "Number unique vertices  : " << get_number_unique_vertices() << std::endl << std::endl;
}

}  // namespace mesh

}  // namespace uepm