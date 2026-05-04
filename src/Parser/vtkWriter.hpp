/**
 * @file vtkWritter.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-11-08
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "mesh.hpp"

namespace uepm {

namespace file {

void write_vtk_geometry(std::ofstream& file, const mesh::mesh& MyMesh);

void write_vtk_data(std::ofstream&                  file,
                    const mesh::mesh&               MyMesh,
                    const std::vector<std::string>& scalar_fields,
                    const std::vector<std::string>& vector_fields,
                    bool                            export_all_fields = false);

void export_as_vtk(const mesh::mesh&               mesh,
                   const std::string&              filename,
                   const std::vector<std::string>& scalar_fields     = {},
                   const std::vector<std::string>& vector_fields     = {},
                   const bool                      export_all_fields = false);

}  // namespace file

}  // namespace uepm