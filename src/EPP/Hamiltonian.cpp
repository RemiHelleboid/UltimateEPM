#include "Hamiltonian.h"

#include <math.h>

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "Hamiltonian.h"
#include "NonLocalFunctional.hpp"
#include "SpinOrbitFunctional.hpp"
#include "Vector3D.h"
#include "physical_constants.hpp"

namespace uepm::pseudopotential {

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
    const Eigen::Index basisSize       = static_cast<Eigen::Index>(m_basisVectors.size());
    const double       latticeConstant = m_material.get_lattice_constant_meter();
    const auto&        pseudopotential = m_material.get_pseudopotential();

    constexpr double one_eight = 1.0 / 8.0;
    Vector3D<double> tau{one_eight * latticeConstant, one_eight * latticeConstant, one_eight * latticeConstant};

    for (Eigen::Index i = 0; i < basisSize; ++i) {
        for (Eigen::Index j = 0; j < basisSize; ++j) {
            m_constant_non_diagonal_matrix(i, j) =
                pseudopotential.GetValue(m_basisVectors[static_cast<std::size_t>(i)] - m_basisVectors[static_cast<std::size_t>(j)],
                                         tau,
                                         latticeConstant);
        }
    }
}

void Hamiltonian::SetMatrix(const Vector3D<double>& k, bool add_non_local_correction, bool enable_soc) {
    const Eigen::Index basisSize       = static_cast<Eigen::Index>(m_basisVectors.size());
    const double       latticeConstant = m_material.get_lattice_constant_meter();
    constexpr double   one_eight       = 1.0 / 8.0;
    Vector3D<double>   tau{one_eight * latticeConstant, one_eight * latticeConstant, one_eight * latticeConstant};

    const double two_pi_over_a = 2.0 * M_PI / latticeConstant;

    // Start from the k-independent local potential (includes diagonal local V)
    matrix = m_constant_non_diagonal_matrix;

    // Precompute k + G for all basis vectors
    std::vector<Vector3D<double>> k_plus_G(static_cast<std::size_t>(basisSize));
    for (Eigen::Index i = 0; i < basisSize; ++i) {
        k_plus_G[static_cast<std::size_t>(i)] = k + m_basisVectors[static_cast<std::size_t>(i)];
    }

    // Diagonal kinetic energy (add on top of local potential; do not overwrite)
    const double diag_factor = (uepm::constants::h_bar * uepm::constants::h_bar) / (2.0 * uepm::constants::m_e * uepm::constants::q_e);

    for (Eigen::Index i = 0; i < basisSize; ++i) {
        const Vector3D<double>& real_k_vector = k_plus_G[static_cast<std::size_t>(i)];
        const double            KG2           = two_pi_over_a * two_pi_over_a * diag_factor * (real_k_vector * real_k_vector);
        matrix(i, i) += std::complex<double>(KG2, 0.0);
    }

    // Non-local correction: build full matrix once (avoid double-adding diagonal)
    if (add_non_local_correction) {
        for (Eigen::Index i = 0; i < basisSize; ++i) {
            const auto& ki = k_plus_G[static_cast<std::size_t>(i)];
            for (Eigen::Index j = 0; j < basisSize; ++j) {
                const auto& kj = k_plus_G[static_cast<std::size_t>(j)];
                auto nl_correction = m_material.compute_pseudopotential_non_local_correction(ki, kj, tau);
                // std::cout << "Non-local correction (" << i << "," << j << "): " << nl_correction << std::endl;
                matrix(i, j) += nl_correction;
            }
        }
    }

    if (enable_soc) {
        // Spin-orbit correction
        SpinOrbitParameters SpinParams = m_material.get_spin_orbit_parameters();
        SpinOrbitCorrection soc_correction(m_material, SpinParams);

        const Eigen::Index N = matrix.rows();

        Eigen::MatrixXcd UpUpMatrix     = matrix;                        // start from current scalar H
        Eigen::MatrixXcd DownDownMatrix = matrix;                        // same for ↓↓ block
        Eigen::MatrixXcd UpDownMatrix   = Eigen::MatrixXcd::Zero(N, N);  // ↑↓
        Eigen::MatrixXcd DownUpMatrix   = Eigen::MatrixXcd::Zero(N, N);  // ↓↑

        for (Eigen::Index i = 0; i < N; ++i) {
            const auto& ki = k_plus_G[static_cast<std::size_t>(i)];
            for (Eigen::Index j = 0; j < N; ++j) {
                const auto&                               kj = k_plus_G[static_cast<std::size_t>(j)];
                Eigen::Matrix<std::complex<double>, 2, 2> soc_contribution =
                    soc_correction.compute_soc_contribution(ki,
                                                            kj,
                                                            m_basisVectors[static_cast<std::size_t>(i)],
                                                            m_basisVectors[static_cast<std::size_t>(j)],
                                                            tau);
                UpUpMatrix(i, j) += soc_contribution(0, 0);
                UpDownMatrix(i, j) += soc_contribution(1, 0);
                DownUpMatrix(i, j) += soc_contribution(0, 1);
                DownDownMatrix(i, j) += soc_contribution(1, 1);
            }
        }

        matrix.resize(2 * N, 2 * N);
        matrix << UpUpMatrix, UpDownMatrix, DownUpMatrix, DownDownMatrix;
    }
}

void Hamiltonian::Diagonalize(bool keep_eigenvectors) {
    // Enforce Hermiticity to mitigate numerical noise prior to diagonalization
    // matrix = 0.5 * (matrix + matrix.adjoint());

    solver.compute(matrix, keep_eigenvectors ? Eigen::ComputeEigenvectors : Eigen::EigenvaluesOnly);
    if (solver.info() != Eigen::Success) {
        std::cout << matrix << std::endl;
        throw std::runtime_error("Eigenvalue decomposition failed");
    }
}

}  // namespace uepm::pseudopotential
