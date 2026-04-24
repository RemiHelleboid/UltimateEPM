/**
 * @file function.hpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief Header file for function class implementation.
 * @version 0.1
 * @date 2021-07-10
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "dataset.hpp"
#include "vector.hpp"

namespace uepm {

namespace mesh {

template <typename T>
using sp_dataset = std::shared_ptr<dataset<T>>;

template <typename data_type>
class function {
 private:
    std::string                        m_name;
    DataType                           m_datatype = DataType::unknown;  //  scalar, vector, tensor ...
    DataLocationType                   m_data_location_type;            // vertex or cell
    std::vector<sp_dataset<data_type>> m_list_sp_datasets{};

 public:
    // Constructors
    function() = delete;
    function(const std::string &function_name, DataType datatype, DataLocationType data_location_type)
        : m_name(function_name),
          m_datatype(datatype),
          m_data_location_type(data_location_type) {}
    function(const std::string &                       function_name,
             DataType                                  datatype,
             DataLocationType                          data_location_type,
             const std::vector<sp_dataset<data_type>> &list_sp_datasets)
        : m_name(function_name),
          m_datatype(datatype),
          m_data_location_type(data_location_type),
          m_list_sp_datasets(list_sp_datasets) {}
    // Copy
    function(const function &other) = default;
    function &operator=(const function &other) = default;
    // Move
    function(function &&other) = default;
    function &operator=(function &&other) = default;

    function get_function_copy(const std::string &new_name) const {
        int      index_dataset = 0;
        function new_function(new_name, m_datatype, m_data_location_type);
        for (const auto &my_sp_dataset : m_list_sp_datasets) {
            index_dataset++;
            sp_dataset<data_type> new_dataset = std::make_shared<dataset<data_type>>(my_sp_dataset->get_dataset_copy(new_name, 0));
            // New values
            auto values = new_dataset->get_values();
            // Mean value
            auto mean_value = std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
            std::cout << "mean_value = " << mean_value << std::endl;
            new_function.add_dataset(new_dataset);
        }
        return new_function;
    }

    // Getters
    const std::string &                       get_name() const { return m_name; }
    DataType                                  get_datatype() const { return m_datatype; }
    DataLocationType                          get_location_type() const { return m_data_location_type; }
    const std::vector<sp_dataset<data_type>> &get_list_sp_datasets() const { return m_list_sp_datasets; }

    std::vector<unsigned int> get_list_valid_regions_index() const {
        std::vector<unsigned int> list_valid_region_index;
        for (const auto &sp_dataset : m_list_sp_datasets) {
            list_valid_region_index.push_back(sp_dataset->get_index_region_validity());
        }
        return list_valid_region_index;
    }

    std::optional<sp_dataset<data_type>> get_dataset_with_region_index(const unsigned int &region_index) const {
        auto it_sp_dataset =
            std::find_if(m_list_sp_datasets.begin(), m_list_sp_datasets.end(), [&region_index](const sp_dataset<data_type> &sp_dataset) {
                return sp_dataset->get_index_region_validity() == region_index;
            });
        if (it_sp_dataset == m_list_sp_datasets.end()) {
            return std::nullopt;
        }
        return *it_sp_dataset;
    }

    // Setters, modifiers
    void set_name(const std::string &new_name) { m_name = new_name; }
    void set_datatype(const DataType &new_datatype) { m_datatype = new_datatype; }

    void add_dataset(sp_dataset<data_type> sp_new_dataset) { m_list_sp_datasets.push_back(sp_new_dataset); }
    void remove_all_datasets() { m_list_sp_datasets.clear(); }

    void add(data_type value) {
        std::for_each(m_list_sp_datasets.begin(), m_list_sp_datasets.end(), [&value](auto &sp_dataset) { sp_dataset->add(value); });
    }

    void add(const function &function_to_add) {
        if (m_datatype != function_to_add.m_datatype) {
            throw std::runtime_error("Error: add() does not support adding two functions with different data types.");
        }
        if (m_data_location_type != function_to_add.m_data_location_type) {
            throw std::runtime_error("Error: add() does not support adding two functions with different data location types.");
        }
        for (auto &sp_dataset : m_list_sp_datasets) {
            std::cout << "Adding " << sp_dataset->get_name() << " to " << m_name << std::endl;
            auto opt_sp_dataset_to_add = function_to_add.get_dataset_with_region_index(sp_dataset->get_index_region_validity());
            if (opt_sp_dataset_to_add.has_value()) {
                sp_dataset->add(*(opt_sp_dataset_to_add.value()));
            } else {
                continue;
            }
        }
    }

    void multiply(const function &function_to_multiply) {
        // if (m_datatype != function_to_multiply.m_datatype) {
        //     throw std::runtime_error("Error: multiply() does not support multiplying two functions with different data types.");
        // }
        // if (m_data_location_type != function_to_multiply.m_data_location_type) {
        //     throw std::runtime_error("Error: multiply() does not support multiplying two functions with different data location types.");
        // }
        const std::string name_other         = function_to_multiply.get_name();
        std::size_t       number_of_datasets = m_list_sp_datasets.size();
        std::cout << "Number of datasets: " << number_of_datasets << std::endl;
        for (auto &my_sp_dataset : m_list_sp_datasets) {
            auto opt_sp_dataset_to_multiply =
                function_to_multiply.get_dataset_with_region_index(my_sp_dataset->get_index_region_validity());
            if (opt_sp_dataset_to_multiply.has_value()) {
                my_sp_dataset->multiply(*(opt_sp_dataset_to_multiply.value()));
            } else {
                std::cout << "No dataset to multiply with " << name_other << std::endl;
                continue;
            }
        }
    }

    void multiply(data_type value) {
        std::for_each(m_list_sp_datasets.begin(), m_list_sp_datasets.end(), [&value](auto &sp_dataset) { sp_dataset->multiply(value); });
    }

    void apply_function(const std::function<data_type(data_type)> &function_to_apply) {
        std::for_each(m_list_sp_datasets.begin(), m_list_sp_datasets.end(), [&function_to_apply](auto &sp_dataset) {
            sp_dataset->apply_function(function_to_apply);
        });
    }

    void print_function_infos() const {
        std::cout << "Function name: " << m_name << std::endl;
        // std::cout << "Function datatype: " << m_datatype << std::endl;
        // std::cout << "Function data location type: " << m_data_location_type << std::endl;
        std::cout << "Number of datasets: " << m_list_sp_datasets.size() << std::endl;
        std::cout << "Datasets: " << std::endl;
        for (const auto &my_sp_dataset : m_list_sp_datasets) {
            my_sp_dataset->print_info();
        }
        std::cout << "----------------------------------------------------------------------\n";
    }
};

}  //  namespace mesh

}  // namespace uepm