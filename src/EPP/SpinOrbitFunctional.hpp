/**
 * @file SpinOrbitFunctions.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-09-30
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include <string>
#include <system_error>
#include <vector>

#include "Material.h"
#include "SpinOrbitParameters.hpp"
#include "Vector3D.h"
#include "physical_constants.hpp"

namespace uepm::pseudopotential {

/**
 * @brief Functor class to compute the spin-orbit Hamiltonian.
 * The name of the functions follow the notation of the paper:
 * Pötz, W. & Vogl, P. Theory of optical-phonon deformation
 * potentials in tetrahedral semiconductors. Phys. Rev. B 24, 2025–2037 (1981).
 *
 */
class SpinOrbitCorrection {
 protected:
    SpinOrbitParameters m_soc_parameters;
    Material            m_material;

 public:
    SpinOrbitCorrection() = delete;
    SpinOrbitCorrection(const Material& material, const SpinOrbitParameters& SpinParams)
        : m_material(material),
          m_soc_parameters(SpinParams){};

    double compute_B2_cation(const Vector3D<double>& K) const;
    double compute_B2_anion(const Vector3D<double>& K) const;
    double compute_B3_cation(const Vector3D<double>& K) const;
    double compute_B3_anion(const Vector3D<double>& K) const;
    double compute_B4_cation(const Vector3D<double>& K) const;
    double compute_B4_anion(const Vector3D<double>& K) const;

    double compute_lambda_1(const Vector3D<double>& K, const Vector3D<double>& Kp) const;
    double compute_lambda_2(const Vector3D<double>& K, const Vector3D<double>& Kp) const;

    double compute_lambda_sym(const Vector3D<double>& K, const Vector3D<double>& Kp) const;
    double compute_lambda_antisym(const Vector3D<double>& K, const Vector3D<double>& Kp) const;

    Eigen::Matrix<std::complex<double>, 2, 2> compute_soc_contribution(const Vector3D<double>& K,
                                                                       const Vector3D<double>& Kp,
                                                                       const Vector3D<double>& G,
                                                                       const Vector3D<double>& Gp,
                                                                       const Vector3D<double>& tau) const;

    static Eigen::Matrix<std::complex<double>, 2, 2> compute_pauli_state_dot_product(const Vector3D<double>& a);
};

}  // namespace uepm::pseudopotential
