/**
 * @file bz_dielectric.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-12-20
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <Eigen/Dense>

#include "Material.h"
#include "bz_mesh.hpp"

namespace bz_mesh {

class BZ_States : public MeshBZ {
 private:
    int m_nb_bands = 0;

    std::vector<Vector3D<int>> m_basisVectors;

    std::vector<Eigen::VectorXd>  m_eigenvalues_k;
    std::vector<Eigen::VectorXd>  m_eigenvalues_k_plus_q;
    std::vector<Eigen::MatrixXcd> m_eigenvectors_k;
    std::vector<Eigen::MatrixXcd> m_eigenvectors_k_plus_q;

 public:
    BZ_States(const EmpiricalPseudopotential::Material& material) : MeshBZ(material) {}

    void set_nb_bands(int nb_bands) { m_nb_bands = nb_bands; }
    void set_basis_vectors(const std::vector<Vector3D<int>>& basis_vectors) { m_basisVectors = basis_vectors; }
    void compute_eigenstates(int nb_threads = 1);
    void compute_dielectric_function(const std::vector<double>& energies, double eta_smearing, int nb_threads = 1);

    double compute_direct_impact_ionization_matrix_element(int                     idx_n1,
                                                           int                     idx_n1_prime,
                                                           int                     idx_n2,
                                                           int                     idx_n2_prime,
                                                           const Vector3D<double>& k1,
                                                           const Vector3D<double>& k2,
                                                           const Vector3D<double>& k1_prime,
                                                           const Vector3D<double>& k2_prime) const;

    void export_full_eigenstates() const;
};

}  // namespace bz_mesh