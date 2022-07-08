#include "Hamiltonian.h"

#define _USE_MATH_DEFINES 1
#include <math.h>

namespace EmpiricalPseudopotential {

Hamiltonian::Hamiltonian(const Material& material, const std::vector<Vector3D<int>>& basisVectors)
    : m_material(material),
      m_basisVectors(basisVectors) {
    const unsigned int basisSize = static_cast<unsigned int>(basisVectors.size());
    // std::cout << "Basis size: " << basisSize << std::endl;
    matrix.resize(basisSize, basisSize);
}

void Hamiltonian::SetMatrix(const Vector3D<double>& k) {
    const unsigned int basisSize = static_cast<unsigned int>(m_basisVectors.size());

    for (unsigned int i = 0; i < basisSize; ++i)
        for (unsigned int j = 0; j < i;++j)  // only the lower triangular of matrix is set because the diagonalization method only needs that
            // off diagonal elements
            matrix(i, j) = m_material.pseudopotential.GetValue(m_basisVectors[i] - m_basisVectors[j]);

    for (unsigned int i = 0; i < basisSize; ++i) {
        // diagonal elements
        // this is actually with 2 * M_PI, but I optimized it with the /2. from the kinetic energy term
        const Vector3D<double> KG = M_PI / m_material.m_a * (k + m_basisVectors[i]);

        matrix(i, i) = std::complex<double>(2. * KG * KG);  // 2* comes from the above optimization, instead of a /2
    }
}

void Hamiltonian::Diagonalize() {
    solver.compute(matrix, Eigen::DecompositionOptions::EigenvaluesOnly);

    assert(solver.info() == Eigen::ComputationInfo::Success);
}

}  // namespace EmpiricalPseudopotential