/**
 * @file geometry_entity.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-04-14
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "geometry_entity.hpp"

namespace uepm {

namespace mesh {

void geometry_entity::add_scalar_data(const std::string &name, double *new_value) {
    const auto it_scalar_data = std::find_if(m_ListScalarData.begin(), m_ListScalarData.end(), [&](const auto &scalar_dataset) {
        return scalar_dataset.m_name == name;
    });
    if (it_scalar_data != m_ListScalarData.end()) {
        it_scalar_data->m_value = new_value;
    } else {
        m_ListScalarData.push_back(scalar_data{name, new_value});
    }
}

double geometry_entity::get_scalar_data(const std::string &fieldname) const {
    const auto it_data = std::find_if(m_ListScalarData.begin(), m_ListScalarData.end(), [&](const auto &scalar_dataset) {
        return scalar_dataset.m_name == fieldname;
    });
    if (it_data != m_ListScalarData.end()) {
        return *(it_data->m_value);
    }
    return std::numeric_limits<double>::quiet_NaN();
}

void geometry_entity::add_vector_data(const std::string &name, vector3 *new_value) {
    const auto it_vector_data = std::find_if(m_ListVectorData.begin(), m_ListVectorData.end(), [&](const auto &vector_dataset) {
        return vector_dataset.m_name == name;
    });
    if (it_vector_data != m_ListVectorData.end()) {
        it_vector_data->m_value = new_value;
    } else {
        m_ListVectorData.push_back(vector_data{name, new_value});
    }
}

vector3 geometry_entity::get_vector_data(const std::string &fieldname) const {
    const auto it_vector_data = std::find_if(m_ListVectorData.begin(), m_ListVectorData.end(), [&](const auto &vector_dataset) {
        return vector_dataset.m_name == fieldname;
    });
    if (it_vector_data != m_ListVectorData.end()) {
        return *(it_vector_data->m_value);
    }
    double NaN = std::numeric_limits<double>::quiet_NaN();
    return vector3(NaN, NaN, NaN);
}

void geometry_entity::remove_scalar_data(const std::string &name) {
    m_ListScalarData.erase(std::remove_if(m_ListScalarData.begin(),
                                          m_ListScalarData.end(),
                                          [&](const auto data_struct) { return data_struct.m_name == name; }),
                           m_ListScalarData.end());
}

void geometry_entity::remove_vector_data(const std::string &name) {
    m_ListVectorData.erase(std::remove_if(m_ListVectorData.begin(),
                                          m_ListVectorData.end(),
                                          [&](const auto data_struct) { return data_struct.m_name == name; }),
                           m_ListVectorData.end());
}

}  // namespace mesh

}  // namespace uepm