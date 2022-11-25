/**
 * @file DielectricFunction.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-11-24
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <atomic>
#include <vector>

#include "Material.h"
#include "SymmetryPoints.h"
#include "Vector3D.h"

namespace EmpiricalPseudopotential {

class DielectricFunction {
 protected:
    std::vector<Vector3D<int>>    m_basisVectors;
    std::vector<Vector3D<double>> m_kpoints;
    const Material&               m_material;
    const int                     m_nb_bands;

 public:
    DielectricFunction() = default;
    DielectricFunction(const Material& material, const std::vector<Vector3D<int>>& basisVectors, const int nb_bands);

    /**
     * @brief Randomly generate a list of k-points in the irreducible wedge of the first Brillouin zone.
     *
     * @param nb_points
     */
    void generate_k_points_random(std::size_t nb_points);

    /**
     * @brief Generate a list of k-points in the irreducible wedge of the first Brillouin zone, using the
     * Monkhorst-Pack algorithm.
     *
     * @param nb_points
     */
    void generate_k_points_grid(std::size_t Nx, std::size_t Ny, std::size_t Nz);

    /** Get the list of k-points.
     * @return const std::vector<Vector3D<double>>&
     */
    const std::vector<Vector3D<double>>& get_kpoints() const { return m_kpoints; }

    /**
     * @brief Compute the dielectric function at a given frequency and q-vector.
     *
     * @param omega
     * @return Eigen::Matrix3cd
     */
    std::vector<double> compute_dielectric_function(const Vector3D<double>&    q_vect,
                                                    const std::vector<double>& list_energies,
                                                    double                     eta_smearing=1e-2,
                                                    int                        nb_threads=1) const;

    void export_kpoints(const std::string& filename) const;
};

}  // namespace EmpiricalPseudopotential