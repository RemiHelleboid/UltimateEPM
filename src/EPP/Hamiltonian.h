#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <complex>

#include "Material.h"

namespace EmpiricalPseudopotential {

class Hamiltonian {
 public:
    Hamiltonian(const Material& material, const std::vector<Vector3D<int>>& basisVectors);
    void SetConstantNonDiagonalMatrix();

    void SetMatrix(const Vector3D<double>& k, bool add_non_local_correction = false);
    void Diagonalize(bool keep_eigenvectors = false);

    const Eigen::VectorXd&  get_eigenvalues() const { return solver.eigenvalues(); }
    const Eigen::MatrixXcd& get_eigenvectors() const { return solver.eigenvectors(); }

    const Eigen::VectorXd& eigenvalues() const { return solver.eigenvalues(); }
    bool                   check_matrix_is_symmetric() const { return matrix.isApprox(matrix.transpose()); }

 protected:
    const Material&                   m_material;
    const std::vector<Vector3D<int>>& m_basisVectors;

    Eigen::MatrixXcd                                m_constant_non_diagonal_matrix;
    Eigen::MatrixXcd                                matrix;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver;
};

}  // namespace EmpiricalPseudopotential
