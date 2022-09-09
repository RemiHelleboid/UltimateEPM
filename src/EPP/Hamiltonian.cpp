#include "Hamiltonian.h"

#define _USE_MATH_DEFINES 1
#include <math.h>

#include <complex>
#include <iostream>
#include <vector>

#include "Constants.hpp"

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
    Vector3D           tau{one_eight * latticeConstant, one_eight * latticeConstant, one_eight * latticeConstant};
    double             factor = 2 * M_PI / latticeConstant;
    for (unsigned int i = 0; i < basisSize; ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            matrix(i, j) = 1.0 * m_material.m_pseudopotential.GetValue(m_basisVectors[i] - m_basisVectors[j], tau, latticeConstant);
            // Non local correction
            if (add_non_local_correction) {
                std::complex<double> nl_correction =
                    m_material.compute_pseudopotential_non_local_correction(k + m_basisVectors[i], k + m_basisVectors[j], tau);
                matrix(i, j) += nl_correction;
            }
        }
    }

    for (unsigned int i = 0; i < basisSize; ++i) {
        // diagonal elements
        const double factor      = 2 * M_PI / latticeConstant;
        const double diag_factor = pow(Constants::h_bar, 2) / (2.0 * Constants::m0 * Constants::q);
        const double KG2         = diag_factor * factor * factor * (k + m_basisVectors[i]) * (k + m_basisVectors[i]);
        matrix(i, i)             = std::complex<double>(KG2);
        // Non local correction
        if (add_non_local_correction) {
            std::complex<double> nl_correction =
                m_material.compute_pseudopotential_non_local_correction(k + m_basisVectors[i], k + m_basisVectors[i], tau);
            std::cout << "nl_correction diag = " << k + m_basisVectors[i] << "   ->   " << nl_correction << std::endl;
            matrix(i, i) += nl_correction;
        }
    }
}

void Hamiltonian::Diagonalize() {
    solver.compute(matrix, Eigen::DecompositionOptions::EigenvaluesOnly);
    assert(solver.info() == Eigen::ComputationInfo::Success);
}

}  // namespace EmpiricalPseudopotential