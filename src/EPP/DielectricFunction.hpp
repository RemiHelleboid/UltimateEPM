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

    std::vector<Vector3D<double>> m_qpoints;
    std::vector<double>           m_energies;

    std::vector<Eigen::VectorXd>  m_eigenvalues_k;
    std::vector<Eigen::MatrixXcd> m_eigenvectors_k;

    /**
     * @brief m_dielectric_function[idx_q][idx_energy] is the dielectric function
     * q = m_qpoints[idx_q]
     * energy = m_energies[idx_energy]
     *
     */
    std::vector<std::vector<double>> m_dielectric_function;

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
    void generate_k_points_grid(std::size_t Nx, std::size_t Ny, std::size_t Nz, double shift, bool irreducible_wedge);

    /** Get the list of k-points.
     * @return const std::vector<Vector3D<double>>&
     */
    const std::vector<Vector3D<double>>& get_kpoints() const { return m_kpoints; }

    /**
     * @brief Set the list of q-points for which the dielectric function will be computed.
     *
     * @param kpoints
     */
    void set_qpoints(const std::vector<Vector3D<double>>& qpoints) { m_qpoints = qpoints; }

    /**
     * @brief Set the list of energies for which the dielectric function will be computed.
     *
     * @param energies
     */
    void set_energies(const std::vector<double>& energies) { m_energies = energies; }

    /**
     * @brief Compute the dielectric function at a given frequency and q-vector.
     *
     * @param omega
     * @return Eigen::Matrix3cd
     */
    void compute_dielectric_function(double eta_smearing = 1e-2);

    /**
     * @brief Get the dielectric function result.
     *
     * @return const std::vector<std::vector<double>>&
     */
    const std::vector<std::vector<double>>& get_dielectric_function() const { return m_dielectric_function; }

    void export_kpoints(const std::string& filename) const;

    void export_dielectric_function_at_q(const std::string& filename, std::size_t idx_q, bool name_auto) const;
};

}  // namespace EmpiricalPseudopotential