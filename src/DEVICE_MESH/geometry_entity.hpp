/**
 * @file geometry_entity.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-04-14
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "vector.hpp"

namespace uepm {

namespace mesh {

struct scalar_data {
    std::string m_name;
    double*     m_value;
    scalar_data(const std::string& name, double* value) : m_name(name), m_value(value) {}
};

struct vector_data {
    std::string m_name;
    vector3*    m_value;
    vector_data(const std::string& name, vector3* value) : m_name(name), m_value(value) {}
};

class geometry_entity {
 protected:
    std::size_t              m_index;
    std::vector<scalar_data> m_ListScalarData;
    std::vector<vector_data> m_ListVectorData;

 public:
    geometry_entity() : m_index{0} {}
    explicit geometry_entity(std::size_t index) : m_index{index} {}

    std::size_t get_index() const { return m_index; }
    void        set_index(std::size_t index) { m_index = index; }

    bool has_scalar_data(const std::string& query_name) const {
        return std::find_if(m_ListScalarData.begin(), m_ListScalarData.end(), [&query_name](const scalar_data& data) {
                   return data.m_name == query_name;
               }) != m_ListScalarData.end();
    }
    bool has_vector_data(const std::string& query_name) const {
        return std::find_if(m_ListVectorData.begin(), m_ListVectorData.end(), [&query_name](const vector_data& data) {
                   return data.m_name == query_name;
               }) != m_ListVectorData.end();
    }

    double  get_scalar_data(const std::string& fieldname) const;
    vector3 get_vector_data(const std::string& fieldname) const;
    void    add_scalar_data(const std::string& name, double* value);
    void    add_vector_data(const std::string& name, vector3* value);
    void    remove_scalar_data(const std::string& name);
    void    remove_vector_data(const std::string& name);
    void    remove_all_scalar_data() {
        m_ListScalarData.clear();
        m_ListScalarData.shrink_to_fit();
    }
    void remove_all_vector_data() {
        m_ListVectorData.clear();
        m_ListVectorData.shrink_to_fit();
    }
    void reset_all_data() {
        remove_all_scalar_data();
        remove_all_vector_data();
    }
};

}  // namespace mesh

}  // namespace uepm