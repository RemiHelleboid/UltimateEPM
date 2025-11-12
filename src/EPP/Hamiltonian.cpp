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

    // Diagonal kinetic energy 
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
                const auto& kj            = k_plus_G[static_cast<std::size_t>(j)];
                auto        nl_correction = m_material.compute_pseudopotential_non_local_correction(ki, kj, tau);
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
                // DBEUG
                // soc_contribution *= 1e-10;
                UpUpMatrix(i, j) += soc_contribution(0, 0);
                UpDownMatrix(i, j) += soc_contribution(0, 1);
                DownUpMatrix(i, j) += soc_contribution(1, 0);
                DownDownMatrix(i, j) += soc_contribution(1, 1);
            }
        }

        matrix.resize(2 * N, 2 * N);
        matrix << UpUpMatrix, UpDownMatrix, DownUpMatrix, DownDownMatrix;
    }
}

void Hamiltonian::Diagonalize(bool keep_eigenvectors) {
    // Enforce Hermiticity to mitigate numerical noise prior to diagonalization
    matrix = 0.5 * (matrix + matrix.adjoint());

    solver.compute(matrix, keep_eigenvectors ? Eigen::ComputeEigenvectors : Eigen::EigenvaluesOnly);
    if (solver.info() != Eigen::Success) {
        std::cout << matrix << std::endl;
        throw std::runtime_error("Eigenvalue decomposition failed");
    }
}

/**
 * @brief Compute the gradient of the Hamiltonian at a specific k-point and energy level using Hellmann-Feynman theorem.
 * 
 * WARING : IT ONLY TAKES INTO ACCOUNT THE LOCAL PART OF THE EPM, IF NL OR SOC ARE ENABLED, THE RESULT IS "WRONG".
 * WARNINAG : IT IS VERY IMPORTANT TO DISENTEGLE THE BANDS BEFORE CALLING THIS FUNCTION, OTHERWISE THE RESULT WILL BE LOCCALLY WRONG .
 * 
 * @param k_point 
 * @param level_index 
 * @return Vector3D<double> 
 */
Vector3D<double> Hamiltonian::compute_gradient_at_level(const Vector3D<double>& k_point, unsigned int level_index) const {
    // Uses Hellmann–Feynman on the kinetic term; local V is k-independent.
    using std::size_t;

    const Eigen::Index Nbasis = static_cast<Eigen::Index>(m_basisVectors.size());
    if (solver.eigenvectors().size() == 0) {
        throw std::runtime_error("compute_gradient_at_level: eigenvectors not available. Call Diagonalize(true) first.");
    }
    if (level_index >= static_cast<unsigned int>(solver.eigenvectors().cols())) {
        throw std::out_of_range("compute_gradient_at_level: level_index out of range");
    }

    const double two_pi_over_a = 2.0 * M_PI / m_material.get_lattice_constant_meter();

    // ħ²/(2 m_e q_e) converts J to eV
    constexpr double diag_factor = (uepm::constants::h_bar * uepm::constants::h_bar) / (2.0 * uepm::constants::m_e * uepm::constants::q_e);

    // For fractional-k derivative (dimensionless): pref_frac = 2 * diag_factor * (2π/a)^2
    const double pref_frac = 2.0 * diag_factor * (two_pi_over_a * two_pi_over_a);

    // Build k+G (fractional components as elsewhere in SetMatrix)
    std::vector<Vector3D<double>> k_plus_G(static_cast<size_t>(Nbasis));
    for (Eigen::Index i = 0; i < Nbasis; ++i) {
        k_plus_G[static_cast<size_t>(i)] = k_point + m_basisVectors[static_cast<size_t>(i)];
    }

    // Eigenvector
    const Eigen::VectorXcd ev       = solver.eigenvectors().col(static_cast<Eigen::Index>(level_index));
    const bool             has_spin = (ev.size() == 2 * Nbasis);

    // Hellmann-Feynman sum
    Vector3D<double> sum_kG{0.0, 0.0, 0.0};
    if (!has_spin) {
        for (Eigen::Index i = 0; i < Nbasis; ++i) {
            sum_kG += k_plus_G[static_cast<size_t>(i)] * std::norm(ev(i));
        }
    } else {
        for (Eigen::Index i = 0; i < Nbasis; ++i) {
            const double w = std::norm(ev(i)) + std::norm(ev(i + Nbasis));
            sum_kG += k_plus_G[static_cast<size_t>(i)] * w;
        }
    }

    Vector3D<double> grad_frac     = sum_kG * pref_frac;
    const double     frac_to_phys  = 1.0 / two_pi_over_a;
    Vector3D<double> grad_phys_eVm = grad_frac * frac_to_phys;

    return grad_phys_eVm;
}

}  // namespace uepm::pseudopotential
