/**
 * @file dataset.hpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief Header file for implementation of dataset related class and struct.
 * @version 0.1
 * @date 2021-07-10
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "vector.hpp"

namespace uepm {

namespace mesh {

/**
 * @brief DataType : different types of "physical fields".
 *
 */
enum class DataType { scalar, vector, tensor, unknown };

std::ostream &operator<<(std::ostream &os, const DataType &);

/**
 * @brief DataLocationType: the kind of location on which the physical data is defined (vertex, cell).
 *
 */
enum class DataLocationType { vertex, cell, unknown };

/**
 * @brief Matching between GMSH and Armin location type (vertex, cell(element)) for datasets.
 *
 */
static const std::map<std::string, DataLocationType> msh_to_armin_data_location_type{{"NodeData", DataLocationType::vertex},
                                                                                     {"ElementData", DataLocationType::cell}};

/**
 * @brief Matching between STF and Armin location type (vertex, cell(element)) for datasets.
 *
 */
static const std::map<int, DataLocationType> STF_to_armin_data_location_type{{0, DataLocationType::vertex}, {3, DataLocationType::cell}};
static const std::map<DataLocationType, int> armin_to_STF_data_location_type{{DataLocationType::vertex, 0}, {DataLocationType::cell, 3}};

/**
 * @brief struct that contains multiple metedata attached to a dataset
 *
 */
struct STF_dataset_infos {
    std::map<std::string, int>              integer_attributes;
    std::map<std::string, std::vector<int>> integer_list_attributes;
    std::map<std::string, double>           double_attributes;
    std::map<std::string, std::string>      string_attributes;

    /**
     * @brief Construct a new STF dataset infos object, empty.
     *
     */
    STF_dataset_infos() {}

    /**
     * @brief add a integer metadata to the struct
     *
     * @param name name of the metadata field
     * @param value int value of the metadata
     */
    void add_integer_attribute(std::string name, int value) { integer_attributes.insert({name, value}); }

    /**
     * @brief add a list of integer metadata to the struct
     *
     * @param name name of the metadata field
     * @param value int value of the metadata
     */
    void add_integer_list_attribute(std::string name, std::vector<int> value) { integer_list_attributes.insert({name, value}); }

    /**
     * @brief add a double metadata to the struct
     *
     * @param name
     * @param value
     */
    void add_double_attribute(const std::string &name, double value) { double_attributes.insert({name, value}); }

    /**
     * @brief add a string metadata to the struct
     *
     * @param name
     * @param value
     */
    void add_string_attribute(const std::string &name, const std::string &value) { string_attributes.insert({name, value}); }

    /**
     * @brief Check if the metadata struct is empty
     *
     * @return true if there is no metadata of any kind
     * @return false if there is at least one metadata of any kind
     */
    bool is_empty() const { return (integer_attributes.empty() && double_attributes.empty() && string_attributes.empty()); }
};

static const STF_dataset_infos empty_dataset_info_struct{};

/**
 * @brief dataset class
 *
 * It represent a dataset for a "physical field", could be ElectricField, temperature, strain etc.
 *
 * @tparam T parameter template for the type of data (double, \ref mesh::vector3 "vector3", ...)
 */
template <typename T>
class dataset {
 protected:
    std::string                           m_name;
    unsigned int                          m_index;
    unsigned int                          m_index_region_validity;
    const std::shared_ptr<std::vector<T>> m_values;
    std::vector<std::size_t>              m_index_geometry_element;
    DataType                              m_data_type;
    DataLocationType                      m_data_location_type;  // vertex or cell
    int                                   m_dimension;           // scalar: 1  , 2d vector: 2, 3d vector: 3  etc...
    STF_dataset_infos                    m_dataset_metadata{};

 public:
    /**
     * @brief Construct a new dataset object : deleted.
     * The default constructor is deleted because an empty dataset should never be instanced.
     *
     */
    dataset() = delete;

    /**
     * @brief Construct a new dataset object with only basic member created (No data and no metadata).
     *
     * @param name name of the dataset (e.g. Temperature)
     * @param index Index of the dataset, strictly positive.
     * @param index_region_validity Index of the region on which the dataset is valid.
     * @param dimension Dimension of the dataset (scalar -> 1, vector -> 2 ou 3 etc.)
     */
    dataset(const std::string &name,
            unsigned int       index,
            unsigned int       index_region_validity,
            DataType           data_type,
            DataLocationType   data_location_type,
            int                dimension)
        : m_name(name),
          m_index(index),
          m_index_region_validity(index_region_validity),
          m_data_type(data_type),
          m_data_location_type(data_location_type),
          m_values(std::make_shared<std::vector<T>>(0)),
          m_dimension(dimension) {}

    /**
     * @brief Construct a new dataset object
     *
     * @param name name of the dataset (e.g. Temperature)
     * @param index Index of the dataset, strictly positive.
     * @param index_region_validity Index of the region on which the dataset is valid.
     * @param values std::vector<T> of values of the dataset
     * @param dimension Dimension of the dataset (scalar -> 1, vector -> 2 ou 3 etc.)
     * @param dataset_metadata Metadata for the dataset (can be empty)
     */
    dataset(const std::string &             name,
            unsigned int                    index,
            unsigned int                    index_region_validity,
            const std::vector<T> &          values,
            const std::vector<std::size_t> &index_geometry_elements,
            DataType                        data_type,
            DataLocationType                data_location_type,
            int                             dimension,
            const STF_dataset_infos &      dataset_metadata = empty_dataset_info_struct)
        : m_name(name),
          m_index(index),
          m_index_region_validity(index_region_validity),
          m_values(std::make_shared<std::vector<T>>(values)),
          m_index_geometry_element(index_geometry_elements),
          m_data_type(data_type),
          m_data_location_type(data_location_type),
          m_dimension(dimension),
          m_dataset_metadata(dataset_metadata) {}

    // Copy constructor
    dataset(const dataset &other) = default;
    // Move constructor
    dataset(dataset &&other) = default;

    dataset get_dataset_copy(const std::string& new_name, unsigned int new_index) const {
        std::vector<T> new_values(*m_values);
        return dataset(new_name,
                       new_index,
                       m_index_region_validity,
                       new_values,
                       m_index_geometry_element,
                       m_data_type,
                       m_data_location_type,
                       m_dimension);
    }

    /**
     * @brief Get the name of the dataset
     *
     * @return std::string
     */
    std::string get_name() const { return m_name; }

    /**
     * @brief Get the index of the dataset
     *
     * @return unsigned int
     */
    unsigned int get_index() const { return m_index; }

    /**
     * @brief Get the index of the dataset
     *
     * @return unsigned int
     */
    void set_index(std::size_t new_index) { m_index = new_index; }

    /**
     * @brief Get the index region validity of the dataset
     *
     * @return unsigned int
     */
    unsigned int get_index_region_validity() const { return m_index_region_validity; }

    std::size_t get_number_values() const { return m_values->size(); }

    /**
     * @brief Get the values vector of the dataset
     *
     * @return std::vector<T>
     */
    const std::vector<T> &get_values() const { return *m_values; }

    /**
     * @brief Get the p values object
     *
     * @return std::shared_ptr<std::vector<T>>
     */
    std::shared_ptr<std::vector<T>> get_p_values() const { return m_values; }

    const std::vector<std::size_t> &get_index_geometry_elements() const { return m_index_geometry_element; }

    /**
     * @brief Get the dataset metadata of the dataset
     *
     * @return const STF_dataset_infos&
     */
    const STF_dataset_infos &get_dataset_metadata() const { return m_dataset_metadata; }

    /**
     * @brief Get the dimension of dataset values
     *
     * @return int
     */
    int get_data_dimension() const { return m_dimension; }

    /**
     * @brief Get the data type object
     *
     * @return DataType
     */
    DataType get_data_type() const { return m_data_type; }

    /**
     * @brief Get the data location type.
     *
     * @return DataLocationType
     */
    DataLocationType get_data_location_type() const { return m_data_location_type; }

    void add(const T &value) {
        std::for_each(m_values->begin(), m_values->end(), [&value](T &v) { v += value; });
    }

    void add(const dataset<T> &other) {
        if (m_dimension != other.get_data_dimension()) {
            throw std::runtime_error("Dimension of dataset to add is not the same as the current one");
        }
        if (m_data_type != other.get_data_type()) {
            throw std::runtime_error("Data type of dataset to add is not the same as the current one");
        }
        if (m_data_location_type != other.get_data_location_type()) {
            throw std::runtime_error("Data location type of dataset to add is not the same as the current one");
        }
        if (m_index_region_validity != other.get_index_region_validity()) {
            throw std::runtime_error("Index region validity of dataset to add is not the same as the current one");
        }
        if (m_values->size() != other.get_values().size()) {
            throw std::runtime_error("Size of dataset to add is not the same as the current one");
        }
        const std::vector<T> &other_values = other.get_values();
        for (std::size_t i = 0; i < m_values->size(); ++i) {
            (*m_values)[i] += other_values[i];
        }

        // std::transform(m_values->begin(), m_values->end(), other_values.begin(), m_values->begin(), std::plus<T>());
    }

    void multiply(const T &value) {
        std::for_each(m_values->begin(), m_values->end(), [&value](T &v) { v *= value; });
    }

    void multiply(const dataset<T> &other) {
        if (m_dimension != other.get_data_dimension()) {
            throw std::runtime_error("Dimension of dataset to multiply is not the same as the current one");
        }
        if (m_data_type != other.get_data_type()) {
            throw std::runtime_error("Data type of dataset to multiply is not the same as the current one");
        }
        if (m_data_location_type != other.get_data_location_type()) {
            throw std::runtime_error("Data location type of dataset to multiply is not the same as the current one");
        }
        if (m_index_region_validity != other.get_index_region_validity()) {
            throw std::runtime_error("Index region validity of dataset to multiply is not the same as the current one");
        }
        if (m_values->size() != other.get_values().size()) {
            throw std::runtime_error("Size of dataset to multiply is not the same as the current one");
        }
        const std::vector<T> &other_values = other.get_values();
        std::transform(m_values->begin(), m_values->end(), other_values.begin(), m_values->begin(), std::multiplies<T>());
    }

    void apply_function(const std::function<T(T)> &function) {
        std::for_each(m_values->begin(), m_values->end(), [&function](T &v) { v = function(v); });
    }

    void apply_function(const std::function<T(T, T)> &function, const dataset<T> &other) {
        std::transform(m_values->begin(), m_values->end(), other.get_values().begin(), m_values->begin(), function);
    }

    /**
     * @brief Print information on the dataset on the standard output.
     *
     */
    void print_info() const;

    void dump_values_in_file(const std::string &filename) const;
};

}  // namespace mesh

}  // namespace uepm