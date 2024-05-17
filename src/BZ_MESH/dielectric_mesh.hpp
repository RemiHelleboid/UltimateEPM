/**
 * @file dielectric_mesh.hpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-05-17
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Material.h"
#include "bz_mesh.hpp"

namespace bz_mesh {

class DielectricMesh : public MeshBZ {
 protected:
    int m_nb_bands = 0;

 public:
    DielectricMesh(const EmpiricalPseudopotential::Material& material) : MeshBZ(material) {}
    void read_dielectric_file(const std::string& filename);
};

}  // namespace bz_mesh