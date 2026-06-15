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

void write_vtk_data(std::ofstream& file, const mesh::mesh& MyMesh);

void export_as_vtk(const mesh::mesh& mesh, const std::string& filename);

void export_as_vtu(const mesh::mesh& mesh, const std::string& filename);

}  // namespace file

}  // namespace uepm