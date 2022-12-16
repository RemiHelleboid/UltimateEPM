#include "Hamiltonian.h"

#define _USE_MATH_DEFINES 1
#include <math.h>

#include <complex>
#include <iostream>
#include <vector>

#include "Constants.hpp"
#include "NonLocalFunctional.hpp"
#include "Vector3D.h"

namespace EmpiricalPseudopotential {

Hamiltonian::Hamiltonian(const Material& material, const std::vector<Vector3D<int>>& basisVectors)
    : m_material(material),
      m_basisVectors(basisVectors) {
    const unsigned int basisSize   = static_cast<unsigned int>(basisVectors.size());
    m_constant_non_diagonal_matrix = Eigen::MatrixXcd::Zero(basisSize, basisSize);
    SetConstantNonDiagonalMatrix();
    matrix.resize(basisSize, basisSize);
}

/**
 * @brief Set the constant part of the Hamiltonian matrix.
 * The non-diagonal part of the matrix does not depend on the k-point, except for the non-local correction.
 *
 */
void Hamiltonian::SetConstantNonDiagonalMatrix() {
    const std::size_t     basisSize       = static_cast<unsigned int>(m_basisVectors.size());
    const double          latticeConstant = m_material.get_lattice_constant_meter();
    const Pseudopotential pseudopotential = m_material.get_pseudopotential();

    constexpr double one_eight = 1.0 / 8.0;
    Vector3D<double> tau{one_eight * latticeConstant, one_eight * latticeConstant, one_eight * latticeConstant};
    for (unsigned int i = 0; i < basisSize; ++i) {
        for (unsigned int j = 0; j < basisSize; ++j) {
            m_constant_non_diagonal_matrix(i, j) = pseudopotential.GetValue(m_basisVectors[i] - m_basisVectors[j], tau, latticeConstant);
        }
    }
}

void Hamiltonian::SetMatrix(const Vector3D<double>& k, bool add_non_local_correction) {
    const unsigned int basisSize       = static_cast<unsigned int>(m_basisVectors.size());
    const double       latticeConstant = m_material.get_lattice_constant_meter();
    constexpr double   one_eight       = 1.0 / 8.0;
    Vector3D<double>   tau{one_eight * latticeConstant, one_eight * latticeConstant, one_eight * latticeConstant};
    const double       fourier_factor = 2.0 * M_PI / latticeConstant;
    matrix = m_constant_non_diagonal_matrix;
    // diagonal elements
    const double diag_factor = pow(Constants::h_bar, 2) / (2.0 * Constants::m0 * Constants::q);
    for (unsigned int i = 0; i < basisSize; ++i) {
        Vector3D<double> real_k_vector = (k + m_basisVectors[i]);
        const double     KG2           = fourier_factor * fourier_factor * diag_factor * (real_k_vector * real_k_vector);
        matrix(i, i)                   = std::complex<double>(KG2, 0);
        // Non local correction
        if (add_non_local_correction) {
            std::complex<double> nl_correction = m_material.compute_pseudopotential_non_local_correction(real_k_vector, real_k_vector, tau);
            matrix(i, i) += nl_correction;
        }
    }

    // const Pseudopotential pseudopotential = m_material.get_pseudopotential();
    if (add_non_local_correction) {
        for (unsigned int i = 0; i < basisSize; ++i) {
            for (unsigned int j = 0; j < basisSize; ++j) {
                std::cout << "NL x \n";
                Vector3D<double> k_vector_i = (k + m_basisVectors[i]);
                Vector3D<double> k_vector_j = (k + m_basisVectors[j]);
                matrix(i, j) += m_material.compute_pseudopotential_non_local_correction(k_vector_i, k_vector_j, tau);
            }
        }
    }
} 

void Hamiltonian::Diagonalize(bool keep_eigenvectors) {
    solver.compute(matrix, keep_eigenvectors ? Eigen::ComputeEigenvectors : Eigen::EigenvaluesOnly);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Eigenvalue decomposition failed");
    }
}

}  // namespace EmpiricalPseudopotential