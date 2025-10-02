// /**
//  * @file NonLocalFunctional.cpp
//  * @author remzerrr (remi.helleboid@gmail.com)
//  * @brief
//  * @version 0.1
//  * @date 2022-10-01
//  *
//  * @copyright Copyright (c) 2022
//  *
//  */

// #include "NonLocalFunctional.hpp"

// #include "Vector3D.h"
// #include "bessel_func.hpp"

// namespace EmpiricalPseudopotential {

// /**
//  * @brief Constructor of the NonLocalFunctional functor.
//  *
//  * @param non_local_parameters
//  * @param material
//  * @param tau
//  */
// NonLocalFunctor::NonLocalFunctor(const NonLocalParameters& non_local_parameters, const Material& material, const Vector3D<double>& tau)
//     : m_non_local_parameters(non_local_parameters),
//       m_material(material),
//       m_tau(tau),
//       m_cinetic_factor{(Constants::h_bar * Constants::h_bar) / (2.0 * Constants::m_e * Constants::q_e)},
//       m_fourrier_factor{2.0 * M_PI / m_material.get_lattice_constant_meter()},
//       m_V0_pref_factor(2.0 * M_PI / m_material.get_lattice_constant_meter()),
//       m_V2_square_well_pref_factor{m_V0_pref_factor},
//       m_V2_gaussian_well_pref_factor{pow(M_PI, 1.5) * (pow(m_non_local_parameters.m_R2_anion, 3.0) / m_material.get_atomic_volume())} {}


// /**
//  * @brief Compute the so called F_l function, which is used in the non-local pseudopotential correction.
//  * (See Chelikowsky, J. R. & Cohen, M. L. Nonlocal pseudopotential calculations for the electronic structure of eleven diamond
//  * and zinc-blende semiconductors. Phys. Rev. B 14, 556–582 (1976).)
//  *
//  * For the values of function, see: Bloomfield, J. K., Face, S. H. P. & Moss, Z. Indefinite Integrals of Spherical Bessel
//  * Functions. Preprint at http://arxiv.org/abs/1703.06428 (2017). Equations 49 and 59.
//  *
//  * @param K1
//  * @param K2
//  * @param atomic_radii
//  * @param l
//  * @return double
//  */
// double NonLocalFunctorF_l_function(const Vector3D<double>& K1, const Vector3D<double>& K2, double atomic_radii, int l) {
//     // This epsilon is used to avoid division by zero in the case of K1 == K2.
//     // The value is quite big, but lower values lead to numerical instabilities (noisy bands).
//     // Reason: K1 and K2 are of the order of  2PI / a_0 ~ 1e10 !
//     constexpr double EPSILON = 1.0e-4;
//     const double     norm_K1 = K1.Length();
//     const double     norm_K2 = K2.Length();
//     if (fabs(norm_K1 - norm_K2) > EPSILON) {
//         const double pre_factor = pow(atomic_radii, 2.0) / (norm_K1 * norm_K1 - norm_K2 * norm_K2);
//         const double F = norm_K1 * generalized_bessel(l + 1, norm_K1 * atomic_radii) * generalized_bessel(l, norm_K2 * atomic_radii) -
//                          norm_K2 * generalized_bessel(l + 1, norm_K2 * atomic_radii) * generalized_bessel(l, norm_K1 * atomic_radii);
//         return pre_factor * F;
//     } else if (norm_K1 > EPSILON) {
//         const double pre_factor = pow(atomic_radii, 3.0) / (2.0);
//         const double F          = pow(generalized_bessel(l, norm_K1 * atomic_radii), 2.0) -
//                          generalized_bessel(l - 1, norm_K1 * atomic_radii) * generalized_bessel(l + 1, norm_K1 * atomic_radii);
//         return pre_factor * F;
//     } else {
//         return (l==0) ? pow(atomic_radii, 3.0) / (3.0) : 0.0;
//     }
// }

// double F_2_function_gaussian(const Vector3D<double>& K1, const Vector3D<double>& K2, double atomic_radii) {
//     const double norm_K1    = K1.Length();
//     const double norm_K2    = K2.Length();
//     const double bessel_arg = 0.5 * (atomic_radii * atomic_radii) * norm_K1 * norm_K2;
//     return bessel_2nd_order_first_kind(bessel_arg) * exp(-0.25 * (norm_K1 * norm_K1 + norm_K2 * norm_K2) * atomic_radii * atomic_radii);
// }

// /**
//  * @brief Compute the non local correction to the EPM Hamiltonian.
//  * It follows: Chelikowsky, J. R. & Cohen, M. L. Nonlocal pseudopotential calculations for the electronic structure of eleven diamond
//  * and zinc-blende semiconductors. Phys. Rev. B 14, 556–582 (1976).
//  * See also: Pötz, W. & Vogl, P. Theory of optical-phonon deformation
//  * potentials in tetrahedral semiconductors. Phys. Rev. B 24, 2025–2037 (1981)
//  *
//  * K1 = (k + G)
//  * K2 = (k + G')
//  * tau = 1/8 * a * (1, 1, 1)
//  *
//  * @warning This function aims to be as close as possible to the original implementation of the authors.
//  * It might not be the most efficient way, even though the compiler may optimize it for us.
//  *
//  * @warning Only square well pseudopotential are supported for now on.
//  *
//  * @param K1
//  * @param K2
//  * @param tau
//  * @return std::complex<double>
//  */
// std::complex<double> Material::compute_pseudopotential_non_local_correction(const Vector3D<double>& K1_normalized,
//                                                                             const Vector3D<double>& K2_normalized,
//                                                                             const Vector3D<double>& tau) const {
//     const double           diag_factor       = pow(Constants::h_bar, 2) / (2.0 * Constants::m_e * Constants::q_e);
//     const double           fourier_factor    = 2.0 * M_PI / get_lattice_constant_meter();
//     const Vector3D<double> G_diff_normalized = (K1_normalized - K2_normalized);
//     const Vector3D<double> K1                = K1_normalized * fourier_factor;
//     const Vector3D<double> K2                = K2_normalized * fourier_factor;
//     const double           norm_K1           = K1.Length();
//     const double           norm_K2           = K2.Length();
//     const double           cos_angle_K1_K2   = compte_cos_angle(K1, K2);
//     const double           V_pre_factor      = 4.0 * M_PI / get_atomic_volume();
//     const double           legendre_0        = 1.0;
//     const double           legendre_2        = 0.5 * (3 * cos_angle_K1_K2 * cos_angle_K1_K2 - 1);

//     // First atomic species: anion
//     double V_anion = 0;
//     // l = 0
//     const double A_0_anion = m_non_local_parameters.m_alpha_0_anion + diag_factor * m_non_local_parameters.m_beta_0_anion *
//                                                                           (norm_K1 * norm_K2 - pow(this->get_fermi_momentum(), 2.0));
//     const double F_0_anion = (m_non_local_parameters.m_R0_anion == 0.0) ? 0.0 : F_l_function(K1, K2, m_non_local_parameters.m_R0_anion, 0);
//     V_anion += V_pre_factor * A_0_anion * (2 * 0 + 1) * 1.0 * F_0_anion;
//     // l = 2
//     double V_anion_2 = 0.0;
//     if (m_non_local_parameters.m_A2_anion != 0) {
//         const double A_2_anion = m_non_local_parameters.m_A2_anion;
//         double       F_2_anion = 0.0;
//         if (m_non_local_parameters.m_well_type == non_local_well_type::square) {
//             F_2_anion = F_l_function(K1, K2, m_non_local_parameters.m_R2_anion, 2);
//             V_anion_2 = V_pre_factor * A_2_anion * (2 * 2 + 1) * legendre_2 * F_2_anion;
//         } else {
//             F_2_anion = F_2_function_gaussian(K1, K2, m_non_local_parameters.m_R2_anion);
//             V_anion_2 = 5.0 * pow(M_PI, 1.5) * (pow(m_non_local_parameters.m_R2_anion, 3.0) / get_atomic_volume()) * A_2_anion *
//                         legendre_2 * F_2_anion;
//         }
//         V_anion += V_anion_2;
//     }

//     // Second atomic species: cation
//     double V_cation = 0;
//     // l = 0
//     const double F_0_cation = (m_non_local_parameters.m_R0_cation == 0) ? 0.0 : F_l_function(K1, K2, m_non_local_parameters.m_R0_cation, 0);
//     const double A_0_cation = m_non_local_parameters.m_alpha_0_cation + m_non_local_parameters.m_beta_0_cation * diag_factor *
//                                                                             (norm_K1 * norm_K2 - pow(this->get_fermi_momentum(), 2.0));
//     V_cation += V_pre_factor * A_0_cation * (2 * 0 + 1) * legendre_0 * F_0_cation;
//     // l = 2
//     double V_cation_2 = 0.0;
//     if (m_non_local_parameters.m_A2_cation != 0.0) {
//         const double A_2_cation = m_non_local_parameters.m_A2_cation;
//         double       F_2_cation = 0.0;
//         if (m_non_local_parameters.m_well_type == non_local_well_type::square) {
//             F_2_cation = F_l_function(K1, K2, m_non_local_parameters.m_R2_cation, 2);
//             V_cation_2 = V_pre_factor * A_2_cation * (2 * 2 + 1) * legendre_2 * F_2_cation;
//         } else {
//             F_2_cation = F_2_function_gaussian(K1, K2, m_non_local_parameters.m_R2_cation);
//             V_cation_2 = 5.0 * pow(M_PI, 1.5) * (pow(m_non_local_parameters.m_R2_cation, 3.0) / get_atomic_volume()) * A_2_cation *
//                          legendre_2 * F_2_cation;
//         }
//         V_cation += V_cation_2;
//     }

//     const double V_symmetric     = 1.0 * (V_anion + V_cation) / 2.0;
//     const double V_antisymmetric = 1.0 * (V_anion - V_cation) / 2.0;

//     constexpr double const_two        = 2.0;
//     const double     lattice_constant = this->get_lattice_constant_meter();
//     const double     Gtau             = (tau / lattice_constant) * (G_diff_normalized);

//     return std::complex<double>(cos(const_two * M_PI * Gtau) * V_symmetric, sin(const_two * M_PI * Gtau) * V_antisymmetric);
// }

// // /**
// //  * @brief Return the value of the F_l function for l = 0, at K1 == K2 == 0.
// //  *
// //  * @param norm_K
// //  * @param atomic_radii
// //  * @return double
// //  */
// // double NonLocalFunctor::F_0_Gamma(double atomic_radii) const { return pow(atomic_radii, 3.0) / 3.0; }

// // /**
// //  * @brief Return the value of the F_l function for l = 0, at K1 == K2 (diagonal term).
// //  *
// //  * @param norm_K1
// //  * @param atomic_radii
// //  * @return double
// //  */
// // double NonLocalFunctor::F_0_diag(double norm_K1, double atomic_radii) const {
// //     constexpr int l          = 0;
// //     const double  pre_factor = pow(atomic_radii, 3.0) / 2.0;
// //     const double  F          = pow(generalized_bessel(l, norm_K1 * atomic_radii), 2.0) -
// //                      generalized_bessel(l - 1, norm_K1 * atomic_radii) * generalized_bessel(l + 1, norm_K1 * atomic_radii);
// //     return m_V0_pref_factor * pre_factor * F;
// // }

// // /**
// //  * @brief Return the value of the F_l function for l = 0, at K1 != K2 (coupling term).
// //  *
// //  * @param norm_K1
// //  * @param norm_K2
// //  * @param atomic_radii
// //  * @return double
// //  */
// // double NonLocalFunctor::F_0_coupling(double norm_K1, double norm_K2, double atomic_radii) const {
// //     constexpr int l          = 0;
// //     const double  pre_factor = pow(atomic_radii, 2.0) / (norm_K1 * norm_K1 - norm_K2 * norm_K2);
// //     const double  F          = norm_K1 * generalized_bessel(l + 1, norm_K1 * atomic_radii) * generalized_bessel(l, norm_K2 * atomic_radii) -
// //                      norm_K2 * generalized_bessel(l + 1, norm_K2 * atomic_radii) * generalized_bessel(l, norm_K1 * atomic_radii);
// //     return m_V0_pref_factor * pre_factor * F;
// // }

// // /**
// //  * @brief Return the value of the F_l function for l = 2, at K1 == K2 == 0 for a square potential.
// //  *
// //  * @param norm_K
// //  * @param atomic_radii
// //  * @return double
// //  */
// // double NonLocalFunctor::F_2_diag_square_potential(double norm_K, double atomic_radii) const {
// //     constexpr int l          = 2;
// //     const double  pre_factor = pow(atomic_radii, 3.0) / 2.0;
// //     const double  F          = pow(generalized_bessel(l, norm_K * atomic_radii), 2.0) -
// //                      generalized_bessel(l - 1, norm_K * atomic_radii) * generalized_bessel(l + 1, norm_K * atomic_radii);
// //     return m_V2_square_well_pref_factor * pre_factor * F;
// // }

// // /**
// //  * @brief Return the value of the F_l function for l = 2, at K1 =! K2 for a square potential.
// //  *
// //  * @param norm_K1
// //  * @param norm_K2
// //  * @param atomic_radii
// //  * @return double
// //  */
// // double NonLocalFunctor::F_2_coupling_square_potential(double norm_K1, double norm_K2, double atomic_radii) const {
// //     constexpr int l          = 2;
// //     const double  pre_factor = pow(atomic_radii, 2.0) / (norm_K1 * norm_K1 - norm_K2 * norm_K2);
// //     const double  F          = norm_K1 * generalized_bessel(l + 1, norm_K1 * atomic_radii) * generalized_bessel(l, norm_K2 * atomic_radii) -
// //                      norm_K2 * generalized_bessel(l + 1, norm_K2 * atomic_radii) * generalized_bessel(l, norm_K1 * atomic_radii);
// //     return m_V2_square_well_pref_factor * pre_factor * F;
// // }

// // /**
// //  * @brief Return the value of the F_l function for l = 2, at any K1, K2 for a gaussian potential.
// //  *
// //  * @param norm_K
// //  * @param atomic_radii
// //  * @return double
// //  */
// // double NonLocalFunctor::F_2_gaussian(double norm_K1, double norm_K2, double atomic_radii) const {
// //     const double bessel_arg = 0.5 * (atomic_radii * atomic_radii) * norm_K1 * norm_K2;
// //     return m_V2_gaussian_well_pref_factor * bessel_2nd_order_first_kind(bessel_arg) *
// //            exp(-0.25 * (norm_K1 * norm_K1 + norm_K2 * norm_K2) * atomic_radii * atomic_radii);
// // }

// // double NonLocalFunctor::compute_anion_non_local_correction(const Vector3D<double>& K1_normalized,
// //                                                            const Vector3D<double>& K2_normalized) const {
// //     if (m_non_local_parameters.m_R0_anion == 0.0) {
// //         return 0.0;
// //     }
// //     const double cos_angle_K1_K2   = compte_cos_angle(K1_normalized, K2_normalized);
// //     const double legendre_0        = 1.0;
// //     const double legendre_2        = 0.5 * (3 * cos_angle_K1_K2 * cos_angle_K1_K2 - 1);
// //     const double norm_K1_normlized = K1_normalized.Length();
// //     const double norm_K2_normlized = K2_normalized.Length();
// //     const double norm_K1           = norm_K1_normlized * m_fourrier_factor;
// //     const double norm_K2           = norm_K2_normlized * m_fourrier_factor;
// //     const double A_0_anion         = m_non_local_parameters.m_alpha_0_anion + m_cinetic_factor * m_non_local_parameters.m_beta_0_anion *
// //                                                                           (norm_K1 * norm_K2 - pow(m_material.get_fermi_momentum(), 2.0));
// //     double F_0 = 0.0;
// //     double F_2 = 0.0;
// //     if (fabs(norm_K1_normlized - norm_K2_normlized) > m_epsilon) {
// //         F_0 = F_0_coupling(norm_K1, norm_K2, m_non_local_parameters.m_R0_anion);
// //         F_2 = (m_non_local_parameters.m_well_type == non_local_well_type::square)
// //                   ? F_2_coupling_square_potential(norm_K1, norm_K2, m_non_local_parameters.m_R0_anion)
// //                   : F_2_gaussian(norm_K1, norm_K2, m_non_local_parameters.m_R0_anion);
// //     } else if (norm_K1_normlized > m_epsilon) {
// //         F_0 = F_0_diag(norm_K1, m_non_local_parameters.m_R0_anion);
// //         F_2 = (m_non_local_parameters.m_well_type == non_local_well_type::square)
// //                   ? F_2_diag_square_potential(norm_K1, m_non_local_parameters.m_R0_anion)
// //                   : F_2_gaussian(norm_K1, norm_K2, m_non_local_parameters.m_R0_anion);
// //     } else if (norm_K1 < m_epsilon) {
// //         F_0 = F_0_Gamma(m_non_local_parameters.m_R0_anion);
// //         F_2 = (m_non_local_parameters.m_well_type == non_local_well_type::square)
// //                   ? F_2_Gamma_square_potential()
// //                   : F_2_gaussian(norm_K1, norm_K2, m_non_local_parameters.m_R0_anion);
// //     }
// //     return A_0_anion * legendre_0 * F_0 + m_non_local_parameters.m_A2_anion * legendre_2 * F_2;
// // }

// // double NonLocalFunctor::compute_cation_non_local_correction(const Vector3D<double>& K1_normalized,
// //                                                             const Vector3D<double>& K2_normalized) const {
// //     if (m_non_local_parameters.m_R0_cation == 0.0) {
// //         return 0.0;
// //     }
// //     const double cos_angle_K1_K2   = compte_cos_angle(K1_normalized, K2_normalized);
// //     const double legendre_0        = 1.0;
// //     const double legendre_2        = 0.5 * (3 * cos_angle_K1_K2 * cos_angle_K1_K2 - 1);
// //     const double norm_K1_normlized = K1_normalized.Length();
// //     const double norm_K2_normlized = K2_normalized.Length();
// //     const double norm_K1           = norm_K1_normlized * m_fourrier_factor;
// //     const double norm_K2           = norm_K2_normlized * m_fourrier_factor;
// //     const double A_0_cation        = m_non_local_parameters.m_alpha_0_cation + m_cinetic_factor * m_non_local_parameters.m_beta_0_cation *
// //                                                                             (norm_K1 * norm_K2 - pow(m_material.get_fermi_momentum(), 2.0));
// //     double F_0 = 0.0;
// //     double F_2 = 0.0;
// //     if (fabs(norm_K1_normlized - norm_K2_normlized) > m_epsilon) {
// //         F_0 = F_0_coupling(norm_K1, norm_K2, m_non_local_parameters.m_R0_cation);
// //         F_2 = (m_non_local_parameters.m_well_type == non_local_well_type::square)
// //                   ? F_2_coupling_square_potential(norm_K1, norm_K2, m_non_local_parameters.m_R0_cation)
// //                   : F_2_gaussian(norm_K1, norm_K2, m_non_local_parameters.m_R0_cation);
// //     } else if (norm_K1_normlized > m_epsilon) {
// //         F_0 = F_0_diag(norm_K1, m_non_local_parameters.m_R0_cation);
// //         F_2 = (m_non_local_parameters.m_well_type == non_local_well_type::square)
// //                   ? F_2_diag_square_potential(norm_K1, m_non_local_parameters.m_R0_cation)
// //                   : F_2_gaussian(norm_K1, norm_K2, m_non_local_parameters.m_R0_cation);
// //     } else if (norm_K1 < m_epsilon) {
// //         F_0 = F_0_Gamma(m_non_local_parameters.m_R0_cation);
// //         F_2 = (m_non_local_parameters.m_well_type == non_local_well_type::square)
// //                   ? F_2_Gamma_square_potential()
// //                   : F_2_gaussian(norm_K1, norm_K2, m_non_local_parameters.m_R0_cation);
// //     }
// //     return A_0_cation * legendre_0 * F_0 + m_non_local_parameters.m_A2_cation * legendre_2 * F_2;
// // }

// // std::complex<double> NonLocalFunctor::operator()(const Vector3D<double>& K1_normalized, const Vector3D<double>& K2_normalized) const {
// //     const Vector3D<double> G_diff_normalized           = K1_normalized - K2_normalized;
// //     const double           anion_non_local_correction  = compute_anion_non_local_correction(K1_normalized, K2_normalized);
// //     const double           cation_non_local_correction = compute_cation_non_local_correction(K1_normalized, K2_normalized);
// //     const double           V_symmetric                 = anion_non_local_correction + cation_non_local_correction;
// //     const double           V_antisymmetric             = anion_non_local_correction - cation_non_local_correction;
// //     constexpr double       const_two                   = 2.0;
// //     const double           lattice_constant            = m_material.get_lattice_constant_meter();
// //     const double           Gtau                        = (m_tau / lattice_constant) * (G_diff_normalized);

// //     return std::complex<double>(cos(const_two * M_PI * Gtau) * V_symmetric, sin(const_two * M_PI * Gtau) * V_antisymmetric);
// // }

// }  // namespace EmpiricalPseudopotential