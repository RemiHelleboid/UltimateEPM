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

namespace uepm::pseudopotential {

double SpinOrbitCorrection::compute_B3_cation(const Vector3D<double>& K) const {
    const double kappa = K.Length() * uepm::constants::bohr_radius * m_soc_parameters.m_radial_extent_cation;  // extent in Bohr radii
    return (5.0 - kappa * kappa) / (5.0 * std::pow(1.0 + kappa * kappa, 4.0));
}

double SpinOrbitCorrection::compute_B3_anion(const Vector3D<double>& K) const {
    const double kappa = K.Length() * uepm::constants::bohr_radius * m_soc_parameters.m_radial_extent_anion;
    return (5.0 - kappa * kappa) / (5.0 * std::pow(1.0 + kappa * kappa, 4.0));
}

double SpinOrbitCorrection::compute_B2_cation(const Vector3D<double>& K) const {
    const double kappa = K.Length() * uepm::constants::bohr_radius * m_soc_parameters.m_radial_extent_cation;  // extent in a0
    return 1.0 / std::pow(1.0 + kappa * kappa, 3.0);
}

double SpinOrbitCorrection::compute_B2_anion(const Vector3D<double>& K) const {
    const double kappa = K.Length() * uepm::constants::bohr_radius * m_soc_parameters.m_radial_extent_anion;
    return 1.0 / std::pow(1.0 + kappa * kappa, 3.0);
}

double SpinOrbitCorrection::compute_B4_cation(const Vector3D<double>& K) const {
    const double kappa = K.Length() * uepm::constants::bohr_radius * m_soc_parameters.m_radial_extent_cation;
    return (5.0 - 3.0 * kappa * kappa) / (5.0 * std::pow(1.0 + kappa * kappa, 5.0));
}

double SpinOrbitCorrection::compute_B4_anion(const Vector3D<double>& K) const {
    const double kappa = K.Length() * uepm::constants::bohr_radius * m_soc_parameters.m_radial_extent_anion;
    return (5.0 - 3.0 * kappa * kappa) / (5.0 * std::pow(1.0 + kappa * kappa, 5.0));
}

double SpinOrbitCorrection::compute_lambda_1(const Vector3D<double>& K, const Vector3D<double>& Kp) const {
    double lambda_1 = m_soc_parameters.m_mu * compute_B3_cation(K) * compute_B3_cation(Kp);
    return lambda_1;
}

double SpinOrbitCorrection::compute_lambda_2(const Vector3D<double>& K, const Vector3D<double>& Kp) const {
    double lambda_2 = m_soc_parameters.m_alpha * m_soc_parameters.m_mu * compute_B3_anion(K) * compute_B3_anion(Kp);
    return lambda_2;
}

double SpinOrbitCorrection::compute_lambda_sym(const Vector3D<double>& K, const Vector3D<double>& Kp) const {
    double lambda_1   = compute_lambda_1(K, Kp);
    double lambda_2   = compute_lambda_2(K, Kp);
    double lambda_sym = (lambda_1 + lambda_2) / 2.0;
    return lambda_sym;
}

double SpinOrbitCorrection::compute_lambda_antisym(const Vector3D<double>& K, const Vector3D<double>& Kp) const {
    double lambda_1       = compute_lambda_1(K, Kp);
    double lambda_2       = compute_lambda_2(K, Kp);
    double lambda_antisym = (lambda_1 - lambda_2) / 2.0;
    return lambda_antisym;
}

Eigen::Matrix<std::complex<double>, 2, 2> SpinOrbitCorrection::compute_pauli_state_dot_product(const Vector3D<double>& myVect) {
    using namespace std::complex_literals;
    std::complex<double>                      a00 = myVect.Z;
    std::complex<double>                      a01 = myVect.X - myVect.Y * 1i;
    std::complex<double>                      a10 = myVect.X + myVect.Y * 1i;
    std::complex<double>                      a11 = -myVect.Z;
    Eigen::Matrix<std::complex<double>, 2, 2> res_matrix;
    res_matrix << a00, a01, a10, a11;
    return res_matrix;
}

Eigen::Matrix<std::complex<double>, 2, 2> SpinOrbitCorrection::compute_soc_contribution(const Vector3D<double>& K,
                                                                                        const Vector3D<double>& Kp,
                                                                                        const Vector3D<double>& G,
                                                                                        const Vector3D<double>& Gp,
                                                                                        const Vector3D<double>& tau) const {
    using namespace std::complex_literals;

    const double a    = m_material.get_lattice_constant_meter();
    const double kfac = (2.0 * M_PI) / a;

    // Physical K for λ (B3 expects 1/m)
    const Vector3D<double> Kphys  = kfac * K;
    const Vector3D<double> Kpphys = kfac * Kp;

    // Elemental Si: antisymmetric = 0
    const double lambda_sym     = compute_lambda_sym(Kphys, Kpphys);
    const double lambda_antisym = 0.0;

    // (K × K')·σ using dimensionless K, then scale by kfac^2 once
    Eigen::Matrix<std::complex<double>, 2, 2> res_matrix = compute_pauli_state_dot_product(cross_product(K, Kp));

    const double               Gtau  = kfac * (tau * (G - Gp));
    const std::complex<double> phase = (-1i * lambda_sym * std::cos(Gtau)) + (lambda_antisym * std::sin(Gtau));

    res_matrix *= (phase * kfac * kfac);  // <-- multiply, don't divide

    return res_matrix;
}
}  // namespace uepm::pseudopotential