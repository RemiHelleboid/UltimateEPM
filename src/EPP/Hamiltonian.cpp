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
    const unsigned int basisSize = static_cast<unsigned int>(basisVectors.size());
    matrix.resize(basisSize, basisSize);
}

void Hamiltonian::SetMatrix(const Vector3D<double>& k, bool add_non_local_correction) {
    const unsigned int basisSize       = static_cast<unsigned int>(m_basisVectors.size());
    const double       latticeConstant = m_material.get_lattice_constant_meter();
    constexpr double   one_eight       = 1.0 / 8.0;
    Vector3D<double>   tau{one_eight * latticeConstant, one_eight * latticeConstant, one_eight * latticeConstant};
    const double       fourier_factor = 2.0 * M_PI / latticeConstant;

    // NonLocalFunctor non_local_functor(m_material.get_non_local_parameters(), m_material, tau);

    const Pseudopotential pseudopotential = m_material.get_pseudopotential();
    for (unsigned int i = 0; i < basisSize; ++i) {
        for (unsigned int j = 0; j < basisSize; ++j) {
            Vector3D<double> k_vector_i = (k + m_basisVectors[i]);
            Vector3D<double> k_vector_j = (k + m_basisVectors[j]);
            matrix(i, j)                = pseudopotential.GetValue(m_basisVectors[i] - m_basisVectors[j], tau, latticeConstant);
            if (add_non_local_correction) {
                matrix(i, j) += m_material.compute_pseudopotential_non_local_correction(k_vector_i, k_vector_j, tau);
            }
        }
    }

    // diagonal elements
    const double diag_factor = pow(Constants::h_bar, 2) / (2.0 * Constants::m0 * Constants::q);
    for (unsigned int i = 0; i < basisSize; ++i) {
        Vector3D<double> real_k_vector = (k + m_basisVectors[i]);
        const double     KG2           = fourier_factor * fourier_factor * diag_factor * real_k_vector * real_k_vector;
        matrix(i, i)                   = std::complex<double>(KG2, 0);
        // Non local correction
        if (add_non_local_correction) {
            std::complex<double> nl_correction = m_material.compute_pseudopotential_non_local_correction(real_k_vector, real_k_vector, tau);
            matrix(i, i) += nl_correction;
        }
    }
}

void Hamiltonian::Diagonalize() {
    solver.compute(matrix, false);
    assert(solver.info() == Eigen::ComputationInfo::Success);
}

}  // namespace EmpiricalPseudopotential