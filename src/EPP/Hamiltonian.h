#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <complex>

#include "Material.h"

namespace EmpiricalPseudopotential {

class Hamiltonian {
 public:
    Hamiltonian(const Material& material, const std::vector<Vector3D<int>>& basisVectors);

    void SetMatrix(const Vector3D<double>& k, bool add_non_local_correction = false);
    void Diagonalize();

    const Eigen::VectorXd& eigenvalues() const { return solver.eigenvalues(); }
    bool                   check_matrix_is_symmetric() const { return matrix.isApprox(matrix.transpose()); }

 protected:
    const Material&                   m_material;
    const std::vector<Vector3D<int>>& m_basisVectors;

    Eigen::MatrixXcd                                matrix;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver;
};

}  // namespace EmpiricalPseudopotential
