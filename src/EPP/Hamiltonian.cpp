#include "Hamiltonian.h"

#define _USE_MATH_DEFINES 1
#include <math.h>

#include <complex>
#include <iostream>
#include <vector>

#include "Constants.hpp"
#include "Vector3D.h"

namespace EmpiricalPseudopotential {

Hamiltonian::Hamiltonian(const Material& material, const std::vector<Vector3D<int>>& basisVectors)
    : m_material(material),
      m_basisVectors(basisVectors) {
    const unsigned int basisSize = static_cast<unsigned int>(basisVectors.size());
    matrix.resize(basisSize, basisSize);
}

void Hamiltonian::SetMatrix(const Vector3D<double>& k, bool add_non_local_correction) {
    const unsigned int basisSize       = static_cast<unsigned int>(m_basisVectors.size());
    const double       latticeConstant = m_material.get_lattice_constant_meter();
    constexpr double   one_eight       = 1.0 / 8.0;
    Vector3D<double>   tau{one_eight * latticeConstant, one_eight * latticeConstant, one_eight * latticeConstant};
    const double       fourier_factor = 2.0 * M_PI / latticeConstant;
    for (unsigned int i = 0; i < basisSize; ++i) {
        for (unsigned int j = 0; j < basisSize; ++j) {
            // matrix(i, j) = 0.0;
            // std::cout << "k + vi = " << (k + m_basisVectors[j]) << std::endl;
            Vector3D<double> real_k_vector_i = (k + m_basisVectors[i]);
            Vector3D<double> real_k_vector_j = (k + m_basisVectors[j]);
            matrix(i, j) = m_material.m_pseudopotential.GetValue(m_basisVectors[i] - m_basisVectors[j], tau, latticeConstant);
            // std::cout << "Local: " << matrix(i, j) << std::endl;
            // Non local correction
            if (add_non_local_correction) {
                // std::cout << "real_k_vector_i: " << real_k_vector_i << std::endl;
                std::complex<double> nl_correction =
                    m_material.compute_pseudopotential_non_local_correction(real_k_vector_i, real_k_vector_j, tau);
                // std::cout << "nl_correction UNDIAG = " << real_k_vector_i << "   ->   " << nl_correction << std::endl;

                matrix(i, j) += nl_correction;
            }
        }
    }

    const double diag_factor = pow(Constants::h_bar, 2) / (2.0 * Constants::m0 * Constants::q);
    for (unsigned int i = 0; i < basisSize; ++i) {
        // diagonal elements
        Vector3D<double> real_k_vector = (k + m_basisVectors[i]);
        const double     KG2           = fourier_factor * fourier_factor * diag_factor * real_k_vector * real_k_vector;
        matrix(i, i)                   = std::complex<double>(KG2, 0);
        // Non local correction
        if (add_non_local_correction) {
            std::complex<double> nl_correction = m_material.compute_pseudopotential_non_local_correction(real_k_vector, real_k_vector, tau);
            // std::cout << "local diag         = " << real_k_vector << "   ->   " << KG2 << std::endl;
            // std::cout << "nl_correction diag = " << real_k_vector << "   ->   " << nl_correction << std::endl;
            matrix(i, i) += nl_correction;
        }
    }
}

void Hamiltonian::Diagonalize() {
    solver.compute(matrix, false);
    assert(solver.info() == Eigen::ComputationInfo::Success);
}

}  // namespace EmpiricalPseudopotential