#include "Hamiltonian.h"

#define _USE_MATH_DEFINES 1
#include <math.h>

namespace EmpiricalPseudopotential {

Hamiltonian::Hamiltonian(const Material& material, const std::vector<Vector3D<int>>& basisVectors)
    : m_material(material),
      m_basisVectors(basisVectors) {
    const unsigned int basisSize = static_cast<unsigned int>(basisVectors.size());
    matrix.resize(basisSize, basisSize);
}

void Hamiltonian::SetMatrix(const Vector3D<double>& k, bool add_non_local_correction) {
    const unsigned int basisSize = static_cast<unsigned int>(m_basisVectors.size());
    constexpr double   one_eight = 1.0 / 8.0;
    const Vector3D     tau{one_eight, one_eight, one_eight};
    for (unsigned int i = 0; i < basisSize; ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            // only the lower triangular of matrix is set because the diagonalization method only needs that
            matrix(i, j) = m_material.m_pseudopotential.GetValue(m_basisVectors[i] - m_basisVectors[j]);
            if (add_non_local_correction) {
                std::complex<double> nl_correction =
                    m_material.compute_pseudopotential_non_local_correction(k + m_basisVectors[i], k + m_basisVectors[j], tau);
                // std::cout << "nl_correction = " << nl_correction << std::endl;
                matrix(i, j) += nl_correction;
            }
        }
    }
    for (unsigned int i = 0; i < basisSize; ++i) {
        // diagonal elements
        // this is actually with 2 * M_PI, but I optimized it with the /2. from the kinetic energy term
        const Vector3D<double> KG        = M_PI / m_material.m_lattice_constant * (k + m_basisVectors[i]);
        constexpr double       const_two = 2.0;
        matrix(i, i) = std::complex<double>(const_two * KG * KG);  // 2* comes from the above optimization, instead of a /2
        if (add_non_local_correction) {
            std::complex<double> nl_correction =
                m_material.compute_pseudopotential_non_local_correction(k + m_basisVectors[i], k + m_basisVectors[i], tau);
            std::cout << "nl_correction diag = " << k + m_basisVectors[i] << "   ->   " << nl_correction << std::endl;
            matrix(i, i) += nl_correction;
        }
    }
    std::cout << "Hamiltonian matrix set" << std::endl;
}

void Hamiltonian::Diagonalize() {
    solver.compute(matrix, Eigen::DecompositionOptions::EigenvaluesOnly);
    assert(solver.info() == Eigen::ComputationInfo::Success);
}

}  // namespace EmpiricalPseudopotential