/**
 * @file grid_data.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-08-05
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <vector>

#include "bbox.hpp"
#include "mesh.hpp"
#include "quadrangular_interpolators.hpp"
#include "rapidcsv.h"
#include "vector.hpp"

namespace uepm {

namespace mesh {

class grid_data {
 private:
    int                 m_dimension = 2;
    std::size_t         m_number_point_x = 1;
    std::size_t         m_number_point_y = 1;
    std::size_t         m_number_point_z = 1;
    std::vector<double> m_x_coordinates;
    std::vector<double> m_y_coordinates;
    std::vector<double> m_z_coordinates;
    std::vector<double> m_data_values;

    bbox m_bbox;

 public:
    grid_data()  = default;
    ~grid_data() = default;

    void compute_bbox() {
        const double x_min = *std::min_element(m_x_coordinates.begin(), m_x_coordinates.end());
        const double x_max = *std::max_element(m_x_coordinates.begin(), m_x_coordinates.end());
        const double y_min = *std::min_element(m_y_coordinates.begin(), m_y_coordinates.end());
        const double y_max = *std::max_element(m_y_coordinates.begin(), m_y_coordinates.end());
        if (m_dimension == 2) {
            m_bbox = bbox(x_min, x_max, y_min, y_max, 0.0, 0.0);
        } else {
            const double z_min = *std::min_element(m_z_coordinates.begin(), m_z_coordinates.end());
            const double z_max = *std::max_element(m_z_coordinates.begin(), m_z_coordinates.end());
            m_bbox             = bbox(x_min, x_max, y_min, y_max, z_min, z_max);
        }
    }

    void set_values(const std::vector<double>& data_values) {
        std::size_t Nz = m_dimension == 3 ? m_number_point_z : 1;
        if (data_values.size() != m_number_point_x * m_number_point_y * Nz) {
            throw std::invalid_argument("Data values size is not equal to the number of points.");
        }
        m_data_values = data_values;
    }
    void set_data_value(const std::size_t i, const std::size_t j, const std::size_t k, const double value) {
        m_data_values[i + j * m_number_point_x + k * m_number_point_x * m_number_point_y] = value;
    }

    double get_data(const std::size_t i, const std::size_t j) const {
        if (i + j * m_number_point_x >= m_data_values.size()) {
            std::cout << "i = " << i << " j = " << j << " m_number_point_x = " << m_data_values.size() << std::endl;
            throw std::invalid_argument("Index out of range.");
        }
        return m_data_values[i + j * m_number_point_x];
    }

    double get_data(const std::size_t i, const std::size_t j, const std::size_t k) const {
        if (m_dimension == 2) {
            return m_data_values[i * m_number_point_y + j];
        } else {
            return m_data_values[i * m_number_point_y * m_number_point_z + j * m_number_point_z + k];
        }
    }

    /**
     * @brief Permute coordinates axis, e.g. (x,y,z) -> (y, x, z)
     *
     */
    void apply_coordinates_permutation(const std::string& permutation) {
        if (permutation == "xy" || permutation == "XY" || permutation == "yx" || permutation == "YX") {
            std::swap(m_x_coordinates, m_y_coordinates);
            std::swap(m_number_point_x, m_number_point_y);
        } else if (permutation == "xz" || permutation == "XZ" || permutation == "zx" || permutation == "ZX") {
            std::swap(m_x_coordinates, m_z_coordinates);
            std::swap(m_number_point_x, m_number_point_z);
        } else if (permutation == "yz" || permutation == "YZ" || permutation == "zy" || permutation == "ZY") {
            std::swap(m_y_coordinates, m_z_coordinates);
            std::swap(m_number_point_y, m_number_point_z);
        } else {
            throw std::invalid_argument("Invalid permutation. Permutation must be one of the following: xy, xz, yz, yx, zx, zy.");
        }
        compute_bbox();
    }

    // void reverse_coordinate_axis(const std::string axis) {
    //     if (axis == "x" || axis == "X") {
    //         std::reverse(m_x_coordinates.begin(), m_x_coordinates.end());
    //     } else if (axis == "y" || axis == "Y") {
    //         std::reverse(m_y_coordinates.begin(), m_y_coordinates.end());
    //     } else if (axis == "z" || axis == "Z") {
    //         std::reverse(m_z_coordinates.begin(), m_z_coordinates.end());
    //     } else {
    //         throw std::invalid_argument("Invalid axis. Axis must be one of the following: x, y, z.");
    //     }
    // }

    void generate_grid(int         dimension,
                       std::size_t number_point_x,
                       std::size_t number_point_y,
                       std::size_t number_point_z,
                       double      x_min,
                       double      x_max,
                       double      y_min,
                       double      y_max,
                       double      z_min,
                       double      z_max) {
        m_dimension      = dimension;
        m_number_point_x = number_point_x;
        m_number_point_y = number_point_y;
        m_number_point_z = 1;

        m_x_coordinates.resize(number_point_x);
        m_y_coordinates.resize(number_point_y);

        if (m_dimension != 2 && m_dimension != 3) {
            throw std::runtime_error("The dimension of the grid must be 2 or 3.");
        }

        if (number_point_x < 2) {
            throw std::invalid_argument("The number of points in the x direction must be greater than 1.");
        }
        if (number_point_y < 2) {
            throw std::invalid_argument("The number of points in the y direction must be greater than 1.");
        }
        if (number_point_z < 2 && dimension == 3) {
            throw std::invalid_argument("The number of points in the z direction must be greater than 1.");
        }

        double dx = (x_max - x_min) / (number_point_x - 1);
        double dy = (y_max - y_min) / (number_point_y - 1);

        for (std::size_t i = 0; i < number_point_x; ++i) {
            m_x_coordinates[i] = x_min + i * dx;
        }

        for (std::size_t i = 0; i < number_point_y; ++i) {
            m_y_coordinates[i] = y_min + i * dy;
        }

        if (m_dimension == 3) {
            m_number_point_z = number_point_z;
            m_z_coordinates.resize(number_point_z);
            double dz = (z_max - z_min) / (number_point_z - 1);
            for (std::size_t i = 0; i < number_point_z; ++i) {
                m_z_coordinates[i] = z_min + i * dz;
            }
        } else {
            m_z_coordinates.resize(1);
            m_z_coordinates[0] = 0.0;
        }
        compute_bbox();
    }

    void generate_grid_from_mesh(const mesh& Mesh, std::size_t number_point_x, std::size_t number_point_y, std::size_t number_point_z) {
        m_dimension = Mesh.get_dimension();
        bbox box    = Mesh.get_bounding_box();
        if (m_dimension == 2) {
            generate_grid(2, number_point_x, number_point_y, 1, box.get_x_min(), box.get_x_max(), box.get_y_min(), box.get_y_max(), 0, 0);
        } else {
            generate_grid(3,
                          number_point_x,
                          number_point_y,
                          number_point_z,
                          box.get_x_min(),
                          box.get_x_max(),
                          box.get_y_min(),
                          box.get_y_max(),
                          box.get_z_min(),
                          box.get_z_max());
        }
        compute_bbox();
    }

    /**
     * @brief Load the data from a csv file.
     *
     * The coordinates must be sorted in the right way which is:
     * (x0, y0, z0), (x0, y0, z1), ..., (x0, y1, z0), (x0, y1, z1),
     * ..., (x1, y0, z0), (x1, y0, z1), ..., (x1, y1, z0), (x1, y1, z1),
     * ...
     *
     *
     * @param filename
     * @param name_column_data
     */
    void load_scalar_from_csv(const std::string& filename,
                              const std::string& name_column_data,
                              const std::string& name_col_x = "X",
                              const std::string& name_col_y = "Y",
                              const std::string& name_col_z = "Z") {
        rapidcsv::Document  doc(filename);
        std::vector<double> x_values    = doc.GetColumn<double>(name_col_x);
        std::vector<double> y_values    = doc.GetColumn<double>(name_col_y);
        std::vector<double> z_values    = doc.GetColumn<double>(name_col_z);
        std::vector<double> data_values = doc.GetColumn<double>(name_column_data);

        std::cout << "Number values X: " << x_values.size() << std::endl;
        std::cout << "Number values Y: " << y_values.size() << std::endl;
        std::cout << "Number values Z: " << z_values.size() << std::endl;
        std::cout << "Number values data: " << data_values.size() << std::endl;

        // Check if the coordinates are sorted in the right way which is:
        // (x0, y0, z0), (x0, y0, z1), ..., (x0, y1, z0), (x0, y1, z1), ..., (x1, y0, z0), (x1, y0, z1), ..., (x1, y1, z0), (x1, y1,
        // z1),
        // ...

        // Get unique x, y, z values
        std::sort(x_values.begin(), x_values.end());
        std::sort(y_values.begin(), y_values.end());
        std::sort(z_values.begin(), z_values.end());
        auto last_x = std::unique(x_values.begin(), x_values.end());
        auto last_y = std::unique(y_values.begin(), y_values.end());
        auto last_z = std::unique(z_values.begin(), z_values.end());
        x_values.erase(last_x, x_values.end());
        y_values.erase(last_y, y_values.end());
        z_values.erase(last_z, z_values.end());
        // Check if sizes are consistant
        std::size_t number_point_x = x_values.size();
        std::size_t number_point_y = y_values.size();
        std::size_t number_point_z = z_values.size();
        // if (data_values.size() != number_point_x * number_point_y * (number_point_z == 1 ? 1 : number_point_z)) {
        //     throw std::runtime_error("The number of data values is not consistant with the number of points.");
        // }
        m_dimension      = (number_point_z == 1) ? 2 : 3;
        m_number_point_x = number_point_x;
        m_number_point_y = number_point_y;
        m_number_point_z = number_point_z;
        m_x_coordinates  = x_values;
        m_y_coordinates  = y_values;
        m_z_coordinates  = z_values;
        m_data_values    = data_values;

        compute_bbox();
        std::cout << "Bounding box: " << m_bbox << std::endl;
    }

    std::array<double, 3> get_coordinates(const std::size_t i, const std::size_t j, const std::size_t k) const {
        if (m_dimension == 2) {
            return {m_x_coordinates[i], m_y_coordinates[j], 0.0};
        } else {
            return {m_x_coordinates[i], m_y_coordinates[j], m_z_coordinates[k]};
        }
    }

    std::optional<std::array<std::size_t, 3>> get_index(const double x, const double y, const double z) const {
        if (!m_bbox.is_inside(vector3(x, y, z))) {
            // std::cout << "Point is outside of the bounding box." << std::endl;
            return std::nullopt;
        }

        std::size_t i = std::distance(m_x_coordinates.begin(), std::lower_bound(m_x_coordinates.begin(), m_x_coordinates.end(), x));
        std::size_t j = std::distance(m_y_coordinates.begin(), std::lower_bound(m_y_coordinates.begin(), m_y_coordinates.end(), y));
        std::size_t k = std::distance(m_z_coordinates.begin(), std::lower_bound(m_z_coordinates.begin(), m_z_coordinates.end(), z));

        if (i > 0) {
            i--;
        }
        if (j > 0) {
            j--;
        }
        if (m_dimension == 3 && k > 0) {
            k--;
        }

        double x_inf = m_x_coordinates[i];
        double x_sup = m_x_coordinates[i + 1];
        double y_inf = m_y_coordinates[j];
        double y_sup = m_y_coordinates[j + 1];

        if (x > x_sup || x < x_inf) {
            std::cout << "Error in x" << std::endl;
            std::cout << "X: " << x_inf << " " << x << " " << x_sup << std::endl;
        }
        if (y > y_sup || y < y_inf) {
            std::cout << "Error in y" << std::endl;
            std::cout << "Y: " << y_inf << " " << y << " " << y_sup << std::endl;
        }

        if (m_dimension == 3) {
            double z_inf = m_z_coordinates[k];
            double z_sup = m_z_coordinates[k + 1];
            if (m_dimension == 3 && (z > z_sup || z < z_inf)) {
                std::cout << "Error in z" << std::endl;
                std::cout << "Z: " << z_inf << " " << z << " " << z_sup << std::endl;
            }
        }

        return std::array<std::size_t, 3>{i, j, k};
    }

    std::optional<std::array<std::size_t, 3>> find_nearest_neighbor(const double x, const double y, const double z) {
        auto op_index = get_index(x, y, z);
        if (!op_index.has_value()) {
            return std::nullopt;
        }
        std::array<std::size_t, 3> index        = op_index.value();
        std::size_t                i            = index[0];
        std::size_t                j            = index[1];
        std::size_t                k            = index[2];
        double                     min_distance = std::numeric_limits<double>::max();
        for (std::size_t ii = 0; ii < 2; ii++) {
            for (std::size_t jj = 0; jj < 2; jj++) {
                for (std::size_t kk = 0; kk < 2; kk++) {
                    double distance = std::pow(m_x_coordinates[i + ii] - x, 2) + std::pow(m_y_coordinates[j + jj] - y, 2) +
                                      std::pow(m_z_coordinates[k + kk] - z, 2);
                    if (distance < min_distance) {
                        min_distance = distance;
                        i            = i + ii;
                        j            = j + jj;
                        k            = k + kk;
                    }
                }
            }
        }
        return std::array<std::size_t, 3>{i, j, k};
    }

    quadrangle_2d get_quadrangle_2d(const std::size_t i, const std::size_t j) const {
        if (m_dimension == 2) {
            vector3 p0(m_x_coordinates[i], m_y_coordinates[j], 0.0);
            vector3 p1(m_x_coordinates[i + 1], m_y_coordinates[j], 0.0);
            vector3 p2(m_x_coordinates[i + 1], m_y_coordinates[j + 1], 0.0);
            vector3 p3(m_x_coordinates[i], m_y_coordinates[j + 1], 0.0);
            quadrangle_2d         my_quad2d(p0, p1, p2, p3);
            my_quad2d.set_scalar_field_values(get_data(i, j), get_data(i + 1, j), get_data(i + 1, j + 1), get_data(i, j + 1));
            return my_quad2d;
        } else {
            throw std::runtime_error("The dimension is not 2.");
        }
    }

    quadrangle_3d get_quadrangle_3d(const std::size_t i, const std::size_t j, const std::size_t k) const {
        // Check if indices ar in bounded
        if (i >= m_number_point_x - 1 || j >= m_number_point_y - 1 || k >= m_number_point_z - 1) {
            throw std::out_of_range("The indices are out of bounds." + std::to_string(i) + " " + std::to_string(j) + " " +
                                    std::to_string(k));
        }

        if (m_dimension == 3) {
            vector3 p0(m_x_coordinates[i], m_y_coordinates[j], m_z_coordinates[k]);
            vector3 p1(m_x_coordinates[i + 1], m_y_coordinates[j], m_z_coordinates[k]);
            vector3 p2(m_x_coordinates[i + 1], m_y_coordinates[j + 1], m_z_coordinates[k]);
            vector3 p3(m_x_coordinates[i], m_y_coordinates[j + 1], m_z_coordinates[k]);
            vector3 p4(m_x_coordinates[i], m_y_coordinates[j], m_z_coordinates[k + 1]);
            vector3 p5(m_x_coordinates[i + 1], m_y_coordinates[j], m_z_coordinates[k + 1]);
            vector3 p6(m_x_coordinates[i + 1], m_y_coordinates[j + 1], m_z_coordinates[k + 1]);
            vector3 p7(m_x_coordinates[i], m_y_coordinates[j + 1], m_z_coordinates[k + 1]);
            quadrangle_3d         my_quad3d(p0, p1, p2, p3, p4, p5, p6, p7);
            my_quad3d.set_scalar_field_values(get_data(i, j, k),
                                              get_data(i + 1, j, k),
                                              get_data(i + 1, j + 1, k),
                                              get_data(i, j + 1, k),
                                              get_data(i, j, k + 1),
                                              get_data(i + 1, j, k + 1),
                                              get_data(i + 1, j + 1, k + 1),
                                              get_data(i, j + 1, k + 1));
            return my_quad3d;
        } else {
            throw std::runtime_error("The dimension is not 3.");
        }
    }

    double multi_linear_interpolate_scalar_at_location(const double x, const double y, const double z) const {
        auto op_index = get_index(x, y, z);
        if (!op_index.has_value()) {
            return 0.0;
        }
        std::array<std::size_t, 3> index = op_index.value();

        if (m_dimension == 2) {
            quadrangle_2d my_quad2d = get_quadrangle_2d(index[0], index[1]);
            return my_quad2d.bilinear_interpolation_scalar_at_location(vector3(x, y, 0.0));
        } else {
            quadrangle_3d my_quad3d = get_quadrangle_3d(index[0], index[1], index[2]);
            return my_quad3d.trilinear_interpolation_scalar_at_location(vector3(x, y, z));
        }
    }

    double multi_linear_interpolate_scalar_at_location(const vector3& location) const {
        return multi_linear_interpolate_scalar_at_location(location.x(), location.y(), location.z());
    }

    double nearest_neighbor_interpolate_scalar_at_location(const double x, const double y, const double z) const {
        auto op_index = get_index(x, y, z);
        if (!op_index.has_value()) {
            return 0.0;
        }
        std::array<std::size_t, 3> index = op_index.value();
        if (m_dimension == 2) {
            quadrangle_2d my_quad2d = get_quadrangle_2d(index[0], index[1]);
            return my_quad2d.nearest_neighbor_scalar_interpolation(vector3(x, y, 0.0));
        } else {
            quadrangle_3d my_quad3d = get_quadrangle_3d(index[0], index[1], index[2]);
            return my_quad3d.nearest_neighbor_scalar_interpolation(vector3(x, y, z));
        }
    }

    double nearest_neighbor_interpolate_scalar_at_location(const vector3& location) const {
        return nearest_neighbor_interpolate_scalar_at_location(location.x(), location.y(), location.z());
    }

    void print_grid_data_info() const {
        std::cout << "The dimension is " << m_dimension << std::endl;
        std::cout << "The number of points in x direction is " << m_number_point_x << std::endl;
        std::cout << "The number of points in y direction is " << m_number_point_y << std::endl;
        std::cout << "The number of points in z direction is " << m_number_point_z << std::endl;
        std::cout << "The number of data values is " << m_data_values.size() << std::endl;
    }
}; 

}  // namespace mesh

}  // namespace uepm