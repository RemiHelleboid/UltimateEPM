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
    const Vector3D<double> kappa = K * (Constants::bohr_radius / m_soc_parameters.m_zeta_cation);
    double                 B2    = 1.0 / std::pow((1.0 + kappa * kappa), 0.33);
    return B2;
}

double SpinOrbitCorrection::compute_B2_anion(const Vector3D<double>& K) const {
    const Vector3D<double> kappa = K * (Constants::bohr_radius / m_soc_parameters.m_zeta_anion);
    double                 B2    = 1.0 / std::pow((1.0 + kappa * kappa), 0.33);
    return B2;
}

double SpinOrbitCorrection::compute_B3_cation(const Vector3D<double>& K) const {
    const Vector3D<double> kappa = K * (Constants::bohr_radius / m_soc_parameters.m_zeta_cation);
    double                 B3    = (5 - kappa * kappa) / (5.0 * pow((1.0 * kappa * kappa), 4.0));
    return B3;
}

double SpinOrbitCorrection::compute_B3_anion(const Vector3D<double>& K) const {
    const Vector3D<double> kappa = K * (Constants::bohr_radius / m_soc_parameters.m_zeta_anion);
    double                 B3    = (5 - kappa * kappa) / (5.0 * pow((1.0 * kappa * kappa), 4.0));
    return B3;
}

double SpinOrbitCorrection::compute_B4_cation(const Vector3D<double>& K) const {
    const Vector3D<double> kappa = K * (Constants::bohr_radius / m_soc_parameters.m_zeta_cation);
    double                 B4    = (5.0 - 3.0 * kappa * kappa) / (5.0 * pow((1.0 + kappa * kappa), 5.0));

    return B4;
}

double SpinOrbitCorrection::compute_B4_anion(const Vector3D<double>& K) const {
    const Vector3D<double> kappa = K * (Constants::bohr_radius / m_soc_parameters.m_zeta_anion);
    double                 B4    = (5.0 - 3.0 * kappa * kappa) / (5.0 * pow((1.0 + kappa * kappa), 5.0));
    return B4;
}



}  // namespace EmpiricalPseudopotential