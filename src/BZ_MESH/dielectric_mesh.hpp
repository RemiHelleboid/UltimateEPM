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

typedef std::complex<double> complex_d;

class DielectricMesh : public MeshBZ {
 protected:
    /**
     * @brief Store the energies of the dielectric function.
     *
     */
    std::vector<double> m_energies;

    /**
     * @brief Store the dielectric function for each k-point of the mesh.
     * m_dielectric_function[idx_node][idx_energy]  is the dielectric function at the k-point idx_node and energy idx_energy (from
     * m_energies).
     *
     */
    std::vector<std::vector<complex_d>> m_dielectric_function;

 public:
    DielectricMesh(const EmpiricalPseudopotential::Material& material) : MeshBZ(material) {}

    /**
     * @brief Read the dielectric function from a .msf file (created by epsilon.epm).
     *
     * @param filename Path to the file containing the dielectric function.
     */
    void read_dielectric_file(const std::string& filename);
};

}  // namespace bz_mesh