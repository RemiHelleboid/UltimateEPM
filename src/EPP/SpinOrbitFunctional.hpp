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

#include <iostream>
#include <string>
#include <system_error>
#include <vector>

#include "Constants.hpp"
#include "Material.h"
#include "SpinOrbitParameters.hpp"
#include "Vector3D.h"

namespace EmpiricalPseudopotential {

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
    SpinOrbitCorrection(const SpinOrbitParameters& parameters, const Material& material)
        : m_soc_parameters(parameters),
          m_material(material) {}

    double compute_B2_cation(const Vector3D<double>& K) const;
    double compute_B2_anion(const Vector3D<double>& K) const;
    double compute_B3_cation(const Vector3D<double>& K) const;
    double compute_B3_anion(const Vector3D<double>& K) const;
    double compute_B4_cation(const Vector3D<double>& K) const;
    double compute_B4_anion(const Vector3D<double>& K) const;

    double compute_lambda_1(const Vector3D<double>& K) const;
    double compute_lambda_2(const Vector3D<double>& K) const;
};

}  // namespace EmpiricalPseudopotential
