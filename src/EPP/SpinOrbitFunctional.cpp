/**
 * @file SpinOrbitDunctional.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-09-30
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "SpinOrbitFunctional.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace EmpiricalPseudopotential {

double SpinOrbitCorrection::compute_B2_cation(const Vector3D<double>& K) const {
    double kappa = K.Length() * (Constants::bohr_radius / m_soc_parameters.m_radial_extent_cation);
    double B2    = 1.0 / std::pow((1.0 + kappa * kappa), 3.0);
    return B2;
}

double SpinOrbitCorrection::compute_B2_anion(const Vector3D<double>& K) const {
    double kappa = K.Length() * (Constants::bohr_radius / m_soc_parameters.m_radial_extent_anion);
    double B2    = 1.0 / std::pow((1.0 + kappa * kappa), 3.0);
    return B2;
}

double SpinOrbitCorrection::compute_B3_cation(const Vector3D<double>& K) const {
    double kappa = K.Length() * (Constants::bohr_radius / m_soc_parameters.m_radial_extent_cation);
    double B3    = (5 - kappa * kappa) / (5.0 * pow((1.0 + kappa * kappa), 4.0));
    return B3;
}

double SpinOrbitCorrection::compute_B3_anion(const Vector3D<double>& K) const {
    double kappa = K.Length() * (Constants::bohr_radius / m_soc_parameters.m_radial_extent_anion);
    double B3    = (5 - kappa * kappa) / (5.0 * pow((1.0 + kappa * kappa), 4.0));
    return B3;
}

double SpinOrbitCorrection::compute_B4_cation(const Vector3D<double>& K) const {
    double kappa = K.Length() * (Constants::bohr_radius / m_soc_parameters.m_radial_extent_cation);
    double B4    = (5.0 - 3.0 * kappa * kappa) / (5.0 * pow((1.0 + kappa * kappa), 5.0));
    return B4;
}

double SpinOrbitCorrection::compute_B4_anion(const Vector3D<double>& K) const {
    double kappa = K.Length() * (Constants::bohr_radius / m_soc_parameters.m_radial_extent_anion);
    double B4    = (5.0 - 3.0 * kappa * kappa) / (5.0 * pow((1.0 + kappa * kappa), 5.0));
    return B4;
}

double SpinOrbitCorrection::compute_lambda_1(const Vector3D<double>& K, const Vector3D<double>& Kp) const {
    double lambda_1 = m_soc_parameters.m_mu * compute_B2_cation(K) * compute_B2_cation(Kp);
    return lambda_1;
}

double SpinOrbitCorrection::compute_lambda_2(const Vector3D<double>& K, const Vector3D<double>& Kp) const {
    double lambda_2 = m_soc_parameters.m_alpha * m_soc_parameters.m_mu * compute_B2_anion(K) * compute_B2_anion(Kp);
    return lambda_2;
}

double SpinOrbitCorrection::compute_lambda_sym(const Vector3D<double>& K, const Vector3D<double>& Kp) const {
    double lambda_1   = compute_lambda_1(K, Kp);
    double lambda_2   = compute_lambda_2(K, Kp);
    double lambda_sym = (lambda_1 + lambda_2) / 2;
    return lambda_sym;
}

double SpinOrbitCorrection::compute_lambda_antisym(const Vector3D<double>& K, const Vector3D<double>& Kp) const {
    double lambda_1       = compute_lambda_1(K, Kp);
    double lambda_2       = compute_lambda_2(K, Kp);
    double lambda_antisym = (lambda_1 - lambda_2) / 2;
    return lambda_antisym;
}

Eigen::Matrix<std::complex<double>, 2, 2> SpinOrbitCorrection::compute_pauli_state_dot_product(Vector3D<double> myVect) {
    using namespace std::complex_literals;
    std::complex<double>  a00 = myVect.Z;
    std::complex<double>  a01 = myVect.X - myVect.Y * 1i;
    std::complex<double>  a10 = myVect.X + myVect.Y * 1i;
    std::complex<double>  a11 = -myVect.Z;
    Eigen::Matrix<std::complex<double>, 2, 2> res_matrix;
    res_matrix << a00, a10, a01, a11;
    return res_matrix;
}

}  // namespace EmpiricalPseudopotential