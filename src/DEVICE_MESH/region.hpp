/**
 * @file region.hpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief Region class header.
 * @version 0.1
 * @date 2021-07-10
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <cassert>
#include <iostream>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "element.hpp"

namespace uepm {

namespace mesh {

enum class RegionType { bulk, interface, contact, unknown };
std::ostream &operator<<(std::ostream &os, const RegionType &);

class region {
 protected:
    int                                             m_dimension    = 0;
    unsigned int                                    m_region_index = 0;
    std::string                                     m_name;
    RegionType                                      m_region_type = RegionType::unknown;
    std::map<std::size_t, std::shared_ptr<element>> m_ListElements;
    std::set<unsigned int>                          m_UniqueVertices;

 public:
    // Constructors
    region() : m_dimension(-1), m_region_index(0), m_name("") { m_ListElements.clear(); }
    region(int dimension, const std::string &name, unsigned int region_index, const RegionType &type)
        : m_dimension(dimension),
          m_region_index(region_index),
          m_name(name),
          m_region_type(type),
          m_ListElements{} {}

    // Setters and Getters
    int                                   get_dimension() const { return m_dimension; }
    void                                  set_dimension(int dimension) { m_dimension = dimension; }
    const std::string &                   get_name() const { return m_name; }
    unsigned int                          get_index() const { return m_region_index; }
    std::size_t                           get_number_elements() const { return m_ListElements.size(); }
    bool                                  is_element_in_region(std::size_t element_index);
    std::optional<const element &>        get_element(std::size_t index_element) const;
    std::shared_ptr<element>              get_p_element(std::size_t index_element);
    std::vector<std::shared_ptr<element>> get_list_elements() const;
    RegionType                            get_region_type() const { return m_region_type; }
    std::vector<std::size_t>              get_list_elements_index() const;
    std::vector<vertex *>                 get_list_all_p_vertices() const;

    // Others
    const std::set<unsigned int> &get_unique_vertices() const { return m_UniqueVertices; }
    std::vector<std::size_t>      get_unique_vertices_as_vector() const;
    std::size_t                   get_number_unique_vertices() const { return m_UniqueVertices.size(); }
    void                          add_element(const std::shared_ptr<element> &sp_elem);
    void                          remove_element(std::size_t element_index);
    void                          remove_elements(const std::vector<std::size_t> &list_element_index);

    void compute_unique_vertices();
    bool is_vertex_inside(std::size_t index_vertex) const;
    bbox compute_bounding_box() const;

    void print_info() const;
};  // namespace mesh

class region_bulk : public region {
 private:
    std::string m_material;

 public:
    // Constructors
    region_bulk(int dimension, const std::string &name, unsigned int region_index, const std::string &material)
        : region(dimension, name, region_index, RegionType::bulk),
          m_material(material) {}
    // Setters and Getters
    std::string get_material() const { return m_material; }
    void        set_material(const std::string &material_name) { m_material = material_name; }
};

class region_interface : public region {
 private:
    int          m_index_region_bulk_0;
    unsigned int m_index_region_bulk_1;

 public:
    // Constructors
    region_interface(int dimension, const std::string &name, unsigned int region_index, int index_bulk_0, unsigned int index_bulk_1)
        : region(dimension, name, region_index, RegionType::interface),
          m_index_region_bulk_0(index_bulk_0),
          m_index_region_bulk_1(index_bulk_1) {}
    // Setters and Getters
    unsigned int        get_bulk_0() const { return m_index_region_bulk_0; }
    unsigned int        get_bulk_1() const { return m_index_region_bulk_1; }
    std::pair<int, int> get_bulk_region_indices() const { return {m_index_region_bulk_0, m_index_region_bulk_1}; }
};

class region_contact : public region {
 private:
    unsigned int m_index_region_bulk_0;

 public:
    // Constructors
    region_contact(int dimension, const std::string &name, unsigned int region_index, unsigned int index_bulk_0)
        : region(dimension, name, region_index, RegionType::contact),
          m_index_region_bulk_0(index_bulk_0) {}
    // Setters and Getters
    unsigned int get_bulk_0() const { return m_index_region_bulk_0; }
};

}  // namespace mesh

}  // namespace uepm