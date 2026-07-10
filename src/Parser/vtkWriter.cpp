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
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
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

void write_vtk_data(std::ofstream& file, const mesh::mesh& MyMesh) {
    auto                                  list_scalar_functions = MyMesh.get_list_scalar_functions();
    auto                                  list_vector_functions = MyMesh.get_list_vector_functions();
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

void export_as_vtk(const mesh::mesh& mesh, const std::string& filename) {
    std::ofstream file;
    file.open(filename);
    write_vtk_geometry(file, mesh);
    write_vtk_data(file, mesh);
    file.close();
}

void export_as_vtu(const mesh::mesh& mesh, const std::string& filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open VTU file: " + filename);
    }

    file << std::setprecision(std::numeric_limits<double>::max_digits10);

    const auto list_vertices      = mesh.get_list_vertices();
    const auto list_bulk_elements = mesh.get_list_bulk_element();

    const std::size_t number_points = list_vertices.size();
    const std::size_t number_cells  = list_bulk_elements.size();

    const int dimension = mesh.get_dimension();
    const int cell_type = (dimension == 2) ? VTK_TRIANGLE : VTK_TETRA;

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << number_points << "\" NumberOfCells=\"" << number_cells << "\">\n";

    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    file << "          ";
    for (const auto& vertex : list_vertices) {
        file << vertex.x() << ' ' << vertex.y() << ' ' << vertex.z() << ' ';
    }
    file << "\n";
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    file << "      <Cells>\n";

    file << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    file << "          ";
    for (const auto& element : list_bulk_elements) {
        for (const auto& vertex_index : element->get_vertices_index()) {
            file << static_cast<long long>(vertex_index) << ' ';
        }
    }
    file << "\n";
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    file << "          ";
    std::size_t offset = 0;
    for (const auto& element : list_bulk_elements) {
        offset += element->get_vertices_index().size();
        file << static_cast<long long>(offset) << ' ';
    }
    file << "\n";
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    file << "          ";
    for (std::size_t i = 0; i < number_cells; ++i) {
        file << cell_type << ' ';
    }
    file << "\n";
    file << "        </DataArray>\n";

    file << "      </Cells>\n";

    const auto list_scalar_functions = mesh.get_list_scalar_functions();
    const auto list_vector_functions = mesh.get_list_vector_functions();

    std::vector<mesh::sp_scalar_function> scalar_vertex_functions;
    std::vector<mesh::sp_scalar_function> scalar_cell_functions;
    std::vector<mesh::sp_vector_function> vector_vertex_functions;
    std::vector<mesh::sp_vector_function> vector_cell_functions;

    for (const auto& scalar_function : list_scalar_functions) {
        if (scalar_function->get_location_type() == mesh::DataLocationType::vertex) {
            scalar_vertex_functions.push_back(scalar_function);
        } else if (scalar_function->get_location_type() == mesh::DataLocationType::cell) {
            scalar_cell_functions.push_back(scalar_function);
        }
    }

    for (const auto& vector_function : list_vector_functions) {
        if (vector_function->get_location_type() == mesh::DataLocationType::vertex) {
            vector_vertex_functions.push_back(vector_function);
        } else if (vector_function->get_location_type() == mesh::DataLocationType::cell) {
            vector_cell_functions.push_back(vector_function);
        }
    }

    file << "      <PointData>\n";

    for (const auto& scalar_function : scalar_vertex_functions) {
        const std::string field_name = scalar_function->get_name();

        file << "        <DataArray type=\"Float64\" Name=\"" << field_name << "\" format=\"ascii\">\n";
        file << "          ";
        for (const auto& vertex : list_vertices) {
            file << vertex.get_scalar_data(field_name) << ' ';
        }
        file << "\n";
        file << "        </DataArray>\n";
    }

    for (const auto& vector_function : vector_vertex_functions) {
        const std::string field_name = vector_function->get_name();

        file << "        <DataArray type=\"Float64\" Name=\"" << field_name
             << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        file << "          ";
        for (const auto& vertex : list_vertices) {
            const auto vector = vertex.get_vector_data(field_name);
            file << vector.x() << ' ' << vector.y() << ' ' << vector.z() << ' ';
        }
        file << "\n";
        file << "        </DataArray>\n";
    }

    file << "      </PointData>\n";

    file << "      <CellData>\n";

    for (const auto& scalar_function : scalar_cell_functions) {
        const std::string field_name = scalar_function->get_name();

        file << "        <DataArray type=\"Float64\" Name=\"" << field_name << "\" format=\"ascii\">\n";
        file << "          ";
        for (const auto& element : list_bulk_elements) {
            file << element->get_scalar_data(field_name) << ' ';
        }
        file << "\n";
        file << "        </DataArray>\n";
    }

    for (const auto& vector_function : vector_cell_functions) {
        const std::string field_name = vector_function->get_name();

        file << "        <DataArray type=\"Float64\" Name=\"" << field_name
             << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        file << "          ";
        for (const auto& element : list_bulk_elements) {
            const auto vector = element->get_vector_data(field_name);
            file << vector.x() << ' ' << vector.y() << ' ' << vector.z() << ' ';
        }
        file << "\n";
        file << "        </DataArray>\n";
    }

    file << "      </CellData>\n";

    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";
}

}  // namespace file

}  // namespace uepm