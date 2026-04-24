/**
 * @file export_vector_to_csv.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-11-10
 *
 * @copyright Copyright (c) 2021
 *
 */
#pragma once

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

namespace uepm {

namespace utils {

template <typename Tvalue>
void export_vector_to_csv(const std::string& filename, const std::string& header, const std::vector<Tvalue>& value_vector) {
    std::ofstream csv_file(filename);
    fmt::print(csv_file, "{}\n", header);
    fmt::print(csv_file, "{:.6e}\n", fmt::join(value_vector, "\n"));
    csv_file.close();
}

inline void export_multiple_vector_to_csv(const std::string&                      filename,
                                          const std::vector<std::string>&         header_columns,
                                          const std::vector<std::vector<double>>& value_vector_of_vector) {
    if (value_vector_of_vector.empty()) {
        fmt::print("Error: Empty value_vector_of_vector. Cannot export to '{}'\n", filename);
        return;
    }

    const std::size_t reference_vector_size = value_vector_of_vector[0].size();
    for (const auto& vector : value_vector_of_vector) {
        if (vector.size() != reference_vector_size) {
            fmt::print("Error: Mismatch between vector sizes in '{}'. Reference size: {}, Found: {}\n", filename, reference_vector_size,
                       vector.size());
            return;
        }
    }

    std::ofstream csv_file(filename);
    if (!csv_file.is_open()) {
        fmt::print("Error: Failed to open file '{}'\n", filename);
        return;
    }

    // Write header_columns
    fmt::print(csv_file, "{}\n", fmt::join(header_columns, ","));

    // Transpose the value_vector_of_vector for easier formatting
    std::vector<std::vector<double>> transposed_values; 
    for (std::size_t i = 0; i < reference_vector_size; ++i) {
        transposed_values.emplace_back(value_vector_of_vector.size());
        for (std::size_t j = 0; j < value_vector_of_vector.size(); ++j) {
            transposed_values[i][j] = value_vector_of_vector[j][i];
        }
    }

    // Write data
    for (const auto& row : transposed_values) {
        fmt::print(csv_file, "{}\n", fmt::join(row, ","));
    }

    csv_file.close();
}

}  // namespace utils

}  // namespace uepm