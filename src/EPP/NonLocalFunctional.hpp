// /**
//  * @file NonLocalFunctional.hpp
//  * @author remzerrr (remi.helleboid@gmail.com)
//  * @brief
//  * @version 0.1
//  * @date 2022-10-01
//  *
//  * @copyright Copyright (c) 2022
//  *
//  */

// #pragma once

// #include <algorithm>
// #include <cmath>
// #include <complex>
// #include <iostream>
// #include <string>
// #include <system_error>
// #include <vector>

// #include "Constants.hpp"
// #include "Material.h"
// #include "NonLocalParameters.hpp"
// #include "Vector3D.h"

// namespace EmpiricalPseudopotential {

// /**
//  * @brief Functor class to compute the non-local correction to the Hamiltonian.
//  * The notation of the functions follow the paper:
//  * Pötz, W. & Vogl, P. Theory of optical-phonon deformation potentials in tetrahedral semiconductors.
//  * Phys. Rev. B 24, 2025–2037 (1981).
//  *
//  */
// class NonLocalFunctor {
//  private:
//     NonLocalParameters      m_non_local_parameters;
//     Material                m_material;
//     Vector3D<double>        m_tau;
//     const double            m_cinetic_factor;
//     const double            m_fourrier_factor;
//     const double            m_V0_pref_factor;
//     const double            m_V2_square_well_pref_factor;
//     const double            m_V2_gaussian_well_pref_factor;
//     static constexpr double m_epsilon = 1e-10;

//  public:
//     NonLocalFunctor() = delete;
//     NonLocalFunctor(const NonLocalParameters& non_local_parameters, const Material& material, const Vector3D<double>& tau);

//     double F_0_Gamma(double atomic_radii) const;
//     double F_0_diag(double norm_K, double atomic_radii) const;
//     double F_0_coupling(double norm_K1, double norm_K2, double atomic_radii) const;

//     double F_2_Gamma_square_potential() const { return 0.0; }
//     double F_2_diag_square_potential(double norm_K, double atomic_radii) const;
//     double F_2_coupling_square_potential(double norm_K1, double norm_K2, double atomic_radii) const;
//     double F_2_gaussian(double norm_K1, double norm_K2, double atomic_radii) const;

//     double compute_anion_non_local_correction(const Vector3D<double>& K1_normalized, const Vector3D<double>& K2_normalized2) const;
//     double compute_cation_non_local_correction(const Vector3D<double>& K1_normalized, const Vector3D<double>& K2_normalized) const;

//     std::complex<double> operator()(const Vector3D<double>& K1_normalized, const Vector3D<double>& K2_normalized) const;
// };

// }  // namespace EmpiricalPseudopotential