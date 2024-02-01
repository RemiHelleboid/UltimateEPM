/**
 * @file bz_states.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-12-20
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "impact_ionization.hpp"

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>

#include "Constants.hpp"
#include "Hamiltonian.h"
#include "bz_mesh.hpp"
#include "bz_states.hpp"
#include "omp.h"

namespace bz_mesh {

double ImpactIonization::compute_impact_ionization_rate(int idx_n1, const Vector3D<double>& k1) {
    double rate                = 0.0;
    int    nb_conduction_bands = 4;
    int    nb_valence_bands    = 4;
    // Pseudo code

    return rate;
}

double ImpactIonization::compute_direct_impact_ionization_matrix_element(int idx_n1,
                                                                         int idx_n1_prime,
                                                                         int idx_n2,
                                                                         int idx_n2_prime,
                                                                         int idx_k1,
                                                                         int idx_k1_prime,
                                                                         int idx_k2_prime) const {
    vector3 k1       = m_list_vertices[idx_k1].get_position();
    vector3 k1_prime = m_list_vertices[idx_k1_prime].get_position();
    vector3 k2_prime = m_list_vertices[idx_k2_prime].get_position();
    vector3 k2       = k1_prime + k2_prime - k1;

    std::size_t index_k2       = get_nearest_k_index(k2);
    std::complex<double>      Ma = 0.0;

    double e_charge       = EmpiricalPseudopotential::Constants::q;
    double eps_zero       = EmpiricalPseudopotential::Constants::eps_zero;
    double volume         = compute_mesh_volume();
    double fourier_factor = 2.0 * M_PI / volume;

    auto        basis_vectors      = get_basis_vectors();
    std::size_t size_basis_vectors = basis_vectors.size();
    for (std::size_t idx_g1 = 0; idx_g1 < size_basis_vectors; ++idx_g1) {
        for (std::size_t idx_g2 = 0; idx_g2 < size_basis_vectors; ++idx_g2) {
            for (std::size_t idx_g1_prime = 0; idx_g1_prime < size_basis_vectors; ++idx_g1_prime) {
                for (std::size_t idx_g2_prime = 0; idx_g2_prime < size_basis_vectors; ++idx_g2_prime) {
                    vector3 g1             = {basis_vectors[idx_g1].X * fourier_factor,
                                              basis_vectors[idx_g1].Y * fourier_factor,
                                              basis_vectors[idx_g1].Z * fourier_factor};
                    vector3 g1_prime       = {basis_vectors[idx_g1_prime].X * fourier_factor,
                                              basis_vectors[idx_g1_prime].Y * fourier_factor,
                                              basis_vectors[idx_g1_prime].Z * fourier_factor};
                    vector3 q_a            = k1_prime + g1_prime - k1 - g1;.0000000000000000000000000000000000000000000000000000000000000000000000
                    auto    state_k1       = m_eigenvectors_k[idx_k1](idx_n1, idx_g1);
                    auto    state_k2       = m_eigenvectors_k[index_k2](idx_n2, idx_g2);
                    auto    state_k1_prime = m_eigenvectors_k[idx_k1_prime](idx_n1_prime, idx_g1_prime);
                    auto    state_k2_prime = m_eigenvectors_k[idx_k2_prime](idx_n2_prime, idx_g2_prime);

                    double               omega_a = m_eigenvalues_k[idx_k1](idx_n1) - m_eigenvalues_k[idx_k1_prime](idx_n1_prime);
                    double               omega_b = m_eigenvalues_k[idx_k2_prime](idx_n2_prime) - m_eigenvalues_k[idx_k1](idx_n1);
                    std::complex<double> micro_matrix_element_a = std::conj(state_k1_prime) * std::conj(state_k2_prime) * state_k1 * state_k2;
                    std::complex<double> dielectric_func      = get_dielectric_function(q_a, omega_a);
                    std::complex<double> pre_factor = e_charge * e_charge / (eps_zero * dielectric_func * q_a.norm() * q_a.norm() * volume);

                    Ma += pre_factor * micro_matrix_element;
                }
            }
        }
    }
}

}  // namespace bz_mesh