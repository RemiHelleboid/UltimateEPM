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

namespace uepm::mesh_bz {

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
    DielectricMesh() = default;
    DielectricMesh(const uepm::pseudopotential::Material& material) : MeshBZ(material) {}

    /**
     * @brief Read the dielectric function from a .msf file (created by epsilon.epm).
     *
     * @param filename Path to the file containing the dielectric function.
     */
    void read_dielectric_file(const std::string& filename);

    /**
     * @brief Find the closest energy in the list of energies.
     * Return the index (idx) of the stored energy directly below the given energy and the fraction (t) of the distance between the two
     * closest energies. The dielectric function at the given energy can be interpolated as: m_dielectric_function[idx_node][idx] * (1 - t)
     * + m_dielectric_function[idx_node][idx + 1] * t
     *
     *
     * @param energy
     * @return std::pair<std::size_t, double>
     */
    std::pair<std::size_t, double> find_closest_energy(double energy) const;

    /**
     * @brief Interpolate the dielectric function at a given k-point and energy.
     *
     * @param k
     * @param energy
     * @return complex_d
     */
    complex_d interpolate_dielectric_function(const vector3& k, double energy) const;
};

}  // namespace uepm::mesh_bz