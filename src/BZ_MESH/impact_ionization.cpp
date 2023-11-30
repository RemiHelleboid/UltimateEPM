/**
 * @file bz_states.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-12-20
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "impact_ionization.hpp"

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>

#include "Hamiltonian.h"
#include "bz_mesh.hpp"
#include "bz_states.hpp"
#include "omp.h"

namespace bz_mesh {

void ImpactIonization::compute_impact_ionization_rate() {
    this->compute_eigenstates();
    auto vertices = this->get_list_vertices();

}

double ImpactIonization::compute_direct_impact_ionization_matrix_element(int                     idx_n1,
                                                                         int                     idx_n1_prime,
                                                                         int                     idx_n2,
                                                                         int                     idx_n2_prime,
                                                                         const Vector3D<double>& k1,
                                                                         const Vector3D<double>& k1_prime,
                                                                         const Vector3D<double>& k2,
                                                                         const Vector3D<double>& k2_prime) const {
    // Momentum conservation
    // Vector3D<double> q = k1 + k1_prime - k2 - k2_prime;
    // if (q.norm() > 1e-12) {
    //     std::cout << "Momentum conservation is not satisfied" << std::endl;
    //     return 0.0;
    // }

    const auto& basis_vectors    = get_basis_vectors();
    std::size_t nb_basis_vectors = basis_vectors.size();
    double      matrix_element   = 0.0;

    Eigen::VectorXd  m_eigenvalues_k1;
    Eigen::VectorXd  m_eigenvalues_k1_prime;
    Eigen::VectorXd  m_eigenvalues_k2;
    Eigen::VectorXd  m_eigenvalues_k2_prime;
    Eigen::MatrixXcd m_eigenvectors_k1;
    Eigen::MatrixXcd m_eigenvectors_k1_prime;
    Eigen::MatrixXcd m_eigenvectors_k2;
    Eigen::MatrixXcd m_eigenvectors_k2_prime;
    // Get the eigenvalues and eigenvectors for the four k-points (TODO)

    double w_a = (m_eigenvalues_k1_prime[idx_n1_prime] - m_eigenvalues_k1[idx_n1]) / EmpiricalPseudopotential::Constants::h_bar;

    // Quadruple sum over the basis vectors
    for (std::size_t idx_g1 = 0; idx_g1 < nb_basis_vectors; ++idx_g1) {
        for (std::size_t idx_g2 = 0; idx_g2 < nb_basis_vectors; ++idx_g2) {
            for (std::size_t idx_g3 = 0; idx_g3 < nb_basis_vectors; ++idx_g3) {
                for (std::size_t idx_g4 = 0; idx_g4 < nb_basis_vectors; ++idx_g4) {
                    auto g1 = basis_vectors[idx_g1];
                    auto g2 = basis_vectors[idx_g2];
                    auto g3 = basis_vectors[idx_g3];
                    auto g4 = basis_vectors[idx_g4];

                    auto A1_G1 = m_eigenvectors_k1(idx_n1, idx_g1);
                }
            }
        }
    }

    return matrix_element;
}

}  // namespace bz_mesh