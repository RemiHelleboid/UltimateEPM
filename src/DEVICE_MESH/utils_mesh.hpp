/**
 * @file utils_mesh.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-04-28
 *
 * @copyright Copyright (c) 2022
 *
 */

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

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include "export_vector_to_csv.hpp"
#include "vector.hpp"

#pragma once

namespace uepm {


namespace utils {

inline std::vector<std::vector<double>> transform_vector_of_position_into_xyz_vectors(
    const std::vector<mesh::vector3> &vector_of_position) {
    std::vector<double> x_coord(vector_of_position.size());
    std::vector<double> y_coord(vector_of_position.size());
    std::vector<double> z_coord(vector_of_position.size());
    std::transform(vector_of_position.begin(), vector_of_position.end(), x_coord.begin(), [&](auto &&position) { return position.x(); });
    std::transform(vector_of_position.begin(), vector_of_position.end(), y_coord.begin(), [&](auto &&position) { return position.y(); });
    std::transform(vector_of_position.begin(), vector_of_position.end(), z_coord.begin(), [&](auto &&position) { return position.z(); });
    return std::vector<std::vector<double>>{x_coord, y_coord, z_coord};
}

inline void export_data_on_grid_to_csv(const std::string                      &filename,
                                       const std::vector<mesh::vector3>       &grid_points,
                                       const std::vector<std::string>         &header_columns,
                                       const std::vector<std::vector<double>> &value_vector_of_vector) {
    if (value_vector_of_vector.empty() || grid_points.empty() ||
        std::any_of(value_vector_of_vector.begin(), value_vector_of_vector.end(),
                    [](const auto &vector_value) { return vector_value.empty(); }) ||
        std::any_of(value_vector_of_vector.begin(), value_vector_of_vector.end(),
                    [&](const auto &vector_value) { return vector_value.size() != value_vector_of_vector[0].size(); })) {
        fmt::print("Error: Invalid data for exporting to '{}'\n", filename);
        return;
    }
    std::size_t reference_vector_size = value_vector_of_vector[0].size();

    std::ofstream csv_file(filename);
    if (!csv_file.is_open()) {
        fmt::print("Error: Failed to open file '{}'\n", filename);
        return;
    }

    // Write header_columns
    fmt::print(csv_file, "{},{},{},{}\n", "X", "Y", "Z", fmt::join(header_columns, ","));

    // Transpose the value_vector_of_vector for easier formatting
    std::vector<std::vector<double>> transposed_values(reference_vector_size, std::vector<double>(value_vector_of_vector.size()));
    for (std::size_t i = 0; i < value_vector_of_vector.size(); ++i) {
        for (std::size_t j = 0; j < reference_vector_size; ++j) {
            transposed_values[j][i] = value_vector_of_vector[i][j];
        }
    }

    // Write data
    for (std::size_t index_value = 0; index_value < value_vector_of_vector[0].size(); ++index_value) {
        fmt::print(csv_file, "{:.6e},{:.6e},{:.6e},", grid_points[index_value].x(), grid_points[index_value].y(),
                   grid_points[index_value].z());
        fmt::print(csv_file, "{:.6e}\n", fmt::join(transposed_values[index_value], ","));
    }

    csv_file.close();
}

inline void export_vector_postion_to_csv(const std::string                &filename,
                                         const std::string                &header,
                                         const std::vector<mesh::vector3> &value_vector) {
    std::vector<double> x_coord(value_vector.size());
    std::vector<double> y_coord(value_vector.size());
    std::vector<double> z_coord(value_vector.size());
    std::transform(value_vector.begin(), value_vector.end(), x_coord.begin(), [&](auto &&position) { return position.x(); });
    std::transform(value_vector.begin(), value_vector.end(), y_coord.begin(), [&](auto &&position) { return position.y(); });
    std::transform(value_vector.begin(), value_vector.end(), z_coord.begin(), [&](auto &&position) { return position.z(); });
    std::vector<std::string> header_columns = {header + "_X", header + "_Y", header + "_Z"};
    export_multiple_vector_to_csv(filename, header_columns, {x_coord, y_coord, z_coord});
}

}  // namespace utils

}  // namespace uepm