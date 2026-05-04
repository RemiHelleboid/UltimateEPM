/**
 * @file vtkWritter.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-11-08
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "vtkWriter.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "mesh.hpp"

namespace uepm {

namespace file {

constexpr int VTK_LINE       = 3;
constexpr int VTK_TRIANGLE   = 5;
constexpr int VTK_TETRA      = 10;
constexpr int VTK_HEXAHEDRON = 12;

int get_number_vertices_per_elements(int dimension) {
    if (dimension == 2) {
        return 3;
    } else if (dimension == 3) {
        return 4;
    } else {
        std::cerr << "Error: dimension not supported" << std::endl;
        throw std::invalid_argument("dimension not supported");
    }
}

void write_vtk_geometry(std::ofstream& file, const mesh::mesh& MyMesh) {
    file << "# vtk DataFile Version 2.1" << std::endl;
    file << "vtk output" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS " << MyMesh.get_nb_vertices() << " double" << std::endl;
    auto list_vertices = MyMesh.get_list_vertices();
    for (const auto& vertex : list_vertices) {
        file << vertex.x() << " " << vertex.y() << " " << vertex.z() << std::endl;
    }
    file << "\n";
    // Cells
    const int   dimension                 = MyMesh.get_dimension();
    auto        list_bulk_elements        = MyMesh.get_list_bulk_element();
    std::size_t nb_cells                  = list_bulk_elements.size();
    int         nb_nodes_per_bulk_element = get_number_vertices_per_elements(dimension);
    file << "CELLS " << nb_cells << " " << nb_cells * (nb_nodes_per_bulk_element + 1) << std::endl;
    for (const auto& element : list_bulk_elements) {
        file << nb_nodes_per_bulk_element << " ";
        for (const auto& vertex_index : element->get_vertices_index()) {
            file << vertex_index << " ";
        }
        file << std::endl;
    }
    // Cell types
    file << "CELL_TYPES " << nb_cells << std::endl;
    const int cell_type = (dimension == 2) ? VTK_TRIANGLE : VTK_TETRA;
    for (std::size_t i = 0; i < nb_cells; ++i) {
        file << cell_type << std::endl;
    }
    file << "\n";
}

void write_vtk_data(std::ofstream&                  file,
                    const mesh::mesh&               MyMesh,
                    const std::vector<std::string>& scalar_fields,
                    const std::vector<std::string>& vector_fields,
                    bool                            export_all_fields) {
    // Scalar fields on vertices
    // auto list_scalar_fields = (export_all_fields) ? MyMesh.get_scalar_functions_name() : scalar_fields;
    auto                            list_scalar_functions = MyMesh.get_list_scalar_functions();
    auto                            list_vector_functions = MyMesh.get_list_vector_functions();
    std::vector<mesh::sp_scalar_function> list_scalar_vertex_function;
    std::vector<mesh::sp_scalar_function> list_scalar_element_function;
    std::vector<mesh::sp_vector_function> list_vector_vertex_function;
    std::vector<mesh::sp_vector_function> list_vector_element_function;
    for (const auto& scalar_function : list_scalar_functions) {
        if (scalar_function->get_location_type() == mesh::DataLocationType::vertex) {
            list_scalar_vertex_function.push_back(scalar_function);
        } else if (scalar_function->get_location_type() == mesh::DataLocationType::cell) {
            list_scalar_element_function.push_back(scalar_function);
        }
    }
    for (const auto& vector_function : list_vector_functions) {
        if (vector_function->get_location_type() == mesh::DataLocationType::vertex) {
            list_vector_vertex_function.push_back(vector_function);
        } else if (vector_function->get_location_type() == mesh::DataLocationType::cell) {
            list_vector_element_function.push_back(vector_function);
        }
    }

    file << "POINT_DATA " << MyMesh.get_nb_vertices() << std::endl;
    auto list_vertices = MyMesh.get_list_vertices();
    for (const auto& my_function : list_scalar_vertex_function) {
        const std::string field_name = my_function->get_name();
        file << "SCALARS " << field_name << " double 1" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;
        for (const auto& vertex : list_vertices) {
            file << vertex.get_scalar_data(field_name) << std::endl;
        }
        file << "\n";
    }
    // Vector fields on vertices
    for (const auto& my_function : list_vector_vertex_function) {
        const std::string field_name = my_function->get_name();
        // std::cout << "Exporting vector field: " << field_name << std::endl;
        file << "VECTORS " << field_name << " double" << std::endl;
        for (const auto& vertex : list_vertices) {
            mesh::vector3 vector = vertex.get_vector_data(field_name);
            file << vector.x() << " " << vector.y() << " " << vector.z() << std::endl;
        }
        file << "\n";
    }
    // Scalar fields on elements
    // Debug file
    file << "CELL_DATA " << MyMesh.get_list_bulk_element().size() << std::endl;
    auto list_bulk_elements = MyMesh.get_list_bulk_element();
    for (const auto& my_function : list_scalar_element_function) {
        const std::string field_name = my_function->get_name();
        file << "SCALARS " << field_name << " double 1" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;
        for (const auto& element : list_bulk_elements) {
            file << element->get_scalar_data(field_name) << std::endl;
        }
        file << "\n";
    }
    // Vector fields on elements
    for (const auto& my_function : list_vector_element_function) {
        const std::string field_name = my_function->get_name();
        // std::cout << "Exporting vector field: " << field_name << std::endl;
        file << "VECTORS " << field_name << " double" << std::endl;
        for (const auto& element : list_bulk_elements) {
            mesh::vector3 vector = element->get_vector_data(field_name);
            file << vector.x() << " " << vector.y() << " " << vector.z() << std::endl;
        }
        file << "\n";
    }

    file.close();
}

void export_as_vtk(const mesh::mesh&               mesh,
                   const std::string&              filename,
                   const std::vector<std::string>& scalar_fields,
                   const std::vector<std::string>& vector_fields,
                   bool                            export_all_fields) {
    std::ofstream file;
    file.open(filename);
    write_vtk_geometry(file, mesh);
    if (!scalar_fields.empty() || !vector_fields.empty() || export_all_fields) {
        write_vtk_data(file, mesh, scalar_fields, vector_fields, export_all_fields);
    }
    file.close();
}

}  // namespace file

}