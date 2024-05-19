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

#include "BandStructure.h"
#include "Constants.hpp"
#include "Hamiltonian.h"
#include "Material.h"
#include "Options.h"
#include "bz_mesh.hpp"
#include "bz_states.hpp"
#include "omp.h"

namespace bz_mesh {

ImpactIonization::ImpactIonization(const EmpiricalPseudopotential::Material& material, const std::string& initial_mesh_path) {
    std::filesystem::path path(initial_mesh_path);
    if (!std::filesystem::exists(path)) {
        std::cerr << "Error: file " << initial_mesh_path << " does not exist." << std::endl;
        exit(1);
    }
    std::ifstream file(initial_mesh_path);
    if (!file.is_open()) {
        std::cerr << "Error: could not open file " << initial_mesh_path << std::endl;
        exit(1);
    }
    m_material          = material;
    m_initial_mesh_path = initial_mesh_path;
}

void ImpactIonization::read_dielectric_file(const std::string& filename) {
    bool normalize_by_fourier_factor = false;
    m_dielectric_mesh.read_mesh_geometry_from_msh_file(filename, normalize_by_fourier_factor);
    m_dielectric_mesh.build_search_tree();
    m_dielectric_mesh.read_dielectric_file(filename);
}

void ImpactIonization::compute_eigenstates() {
    int                nb_bands_to_use = 4;
    bz_mesh::BZ_States my_bz_mesh(m_material);
    my_bz_mesh.set_nb_bands(nb_bands_to_use);
    EmpiricalPseudopotential::BandStructure band_structure{};
    int                                     nb_threads           = 1;
    int                                     nb_nearest_neighbors = 10;
    bool                                    nonlocal_epm         = false;
    bool                                    enable_soc           = false;
    band_structure.Initialize(m_material, nb_bands_to_use, {}, nb_nearest_neighbors, nonlocal_epm, enable_soc);
    auto basis = band_structure.get_basis_vectors();
    my_bz_mesh.set_basis_vectors(basis);

    my_bz_mesh.read_mesh_geometry_from_msh_file(m_initial_mesh_path);

    const vector3 b1 = {-1.0, 1.0, 1.0};
    const vector3 b2 = {1.0, -1.0, 1.0};
    const vector3 b3 = {1.0, 1.0, -1.0};

    // std::cout << "k: " << k << std::endl;

    // test k + G
    std::vector<int> list_n_k = {0, 1, -1, 2, -2, 3, -3, 4, -4};
    // std::vector<int> list_n_k = {0, 1, -1, 2, -2};

    std::cout << "Max radius G0 BZ: " << m_max_radius_G0_BZ << std::endl;
    for (auto&& n_k_x : list_n_k) {
        for (auto&& n_k_y : list_n_k) {
            for (auto&& n_k_z : list_n_k) {
                vector3 G_BZ = n_k_x * b1 + n_k_y * b2 + n_k_z * b3;
                if (G_BZ.norm() > m_max_radius_G0_BZ) {
                    continue;
                }
                std::cout << "G_BZ: " << G_BZ << std::endl;
                auto ptr_BZ_states = std::make_unique<BZ_States>(m_material);
                ptr_BZ_states->set_nb_bands(nb_bands_to_use);
                ptr_BZ_states->set_basis_vectors(basis);
                ptr_BZ_states->read_mesh_geometry_from_msh_file(m_initial_mesh_path);
                ptr_BZ_states->shift_bz_center(G_BZ);
                ptr_BZ_states->compute_eigenstates(nb_threads);
                m_list_BZ_states.push_back(std::move(ptr_BZ_states));
            }
        }
    }
    std::cout << "Number of BZ states: " << m_list_BZ_states.size() << std::endl;
}

double ImpactIonization::compute_impact_ionization_rate(int idx_n1, const Vector3D<double>& k1) {
    double rate                = 0.0;
    int    nb_conduction_bands = 4;
    int    nb_valence_bands    = 4;
    // Pseudo code

    return rate;
}

// /**
//  * @brief Compute the direct impact ionization matrix element (Ma) for a given set of indices.
//  * k2 is not an input parameter, it is computed from k1, k1_prime and k2_prime, using the conservation of momentum.
//  * k2 = k1_prime + k2_prime - k1
//  *
//  *
//  * @param idx_n1
//  * @param idx_n1_prime
//  * @param idx_n2
//  * @param idx_n2_prime
//  * @param idx_k1
//  * @param idx_k1_prime
//  * @param idx_k2_prime
//  * @return double
//  */
// std::array<complex_d, 2> ImpactIonization::compute_direct_indirect_impact_ionization_matrix_element(int idx_n1,
//                                                                                                     int idx_n1_prime,
//                                                                                                     int idx_n2,
//                                                                                                     int idx_n2_prime,
//                                                                                                     int idx_k1,
//                                                                                                     int idx_k1_prime,
//                                                                                                     int idx_k2,
//                                                                                                     int idx_k2_prime) const {
//     vector3   k1       = m_list_vertices[idx_k1].get_position();
//     vector3   k2       = m_list_vertices[idx_k2].get_position();
//     vector3   k1_prime = m_list_vertices[idx_k1_prime].get_position();
//     vector3   k2_prime = m_list_vertices[idx_k2_prime].get_position();
//     complex_d Ma       = 0.0;
//     complex_d Mb       = 0.0;

//     double e_charge       = EmpiricalPseudopotential::Constants::q;
//     double eps_zero       = EmpiricalPseudopotential::Constants::eps_zero;
//     double volume         = compute_mesh_volume();
//     double fourier_factor = 2.0 * M_PI / volume;

//     auto        basis_vectors      = get_basis_vectors();
//     std::size_t size_basis_vectors = basis_vectors.size();
//     for (std::size_t idx_g1 = 0; idx_g1 < size_basis_vectors; ++idx_g1) {
//         for (std::size_t idx_g2 = 0; idx_g2 < size_basis_vectors; ++idx_g2) {
//             for (std::size_t idx_g1_prime = 0; idx_g1_prime < size_basis_vectors; ++idx_g1_prime) {
//                 for (std::size_t idx_g2_prime = 0; idx_g2_prime < size_basis_vectors; ++idx_g2_prime) {
//                     vector3 g1             = {basis_vectors[idx_g1].X * fourier_factor,
//                                               basis_vectors[idx_g1].Y * fourier_factor,
//                                               basis_vectors[idx_g1].Z * fourier_factor};
//                     vector3 g1_prime       = {basis_vectors[idx_g1_prime].X * fourier_factor,
//                                               basis_vectors[idx_g1_prime].Y * fourier_factor,
//                                               basis_vectors[idx_g1_prime].Z * fourier_factor};
//                     vector3 g2_prime       = {basis_vectors[idx_g2_prime].X * fourier_factor,
//                                               basis_vectors[idx_g2_prime].Y * fourier_factor,
//                                               basis_vectors[idx_g2_prime].Z * fourier_factor};
//                     vector3 q_a            = k1_prime + g1_prime - k1 - g1;
//                     vector3 q_b            = k2_prime + g2_prime - k1 - g1;
//                     auto    state_k1       = m_eigenvectors_k[idx_k1](idx_n1, idx_g1);
//                     auto    state_k2       = m_eigenvectors_k[idx_k2](idx_n2, idx_g2);
//                     auto    state_k1_prime = m_eigenvectors_k[idx_k1_prime](idx_n1_prime, idx_g1_prime);
//                     auto    state_k2_prime = m_eigenvectors_k[idx_k2_prime](idx_n2_prime, idx_g2_prime);

//                     double    omega_a           = m_eigenvalues_k[idx_k1](idx_n1) - m_eigenvalues_k[idx_k1_prime](idx_n1_prime);
//                     double    omega_b           = m_eigenvalues_k[idx_k2_prime](idx_n2_prime) - m_eigenvalues_k[idx_k1](idx_n1);
//                     complex_d dielectric_func_a = get_dielectric_function(q_a, omega_a);
//                     complex_d pre_factor_a      = e_charge * e_charge / (eps_zero * dielectric_func_a * q_a.norm() * q_a.norm() *
//                     volume); complex_d dielectric_func_b = get_dielectric_function(q_b, omega_b); complex_d core_overlap      =
//                     std::conj(state_k1_prime) * std::conj(state_k2_prime) * state_k1 * state_k2; complex_d micro_matrix_element_a =
//                     core_overlap * pre_factor_a; complex_d pre_factor_b = e_charge * e_charge / (eps_zero * dielectric_func_b *
//                     q_b.norm() * q_b.norm() * volume); complex_d micro_matrix_element_b = core_overlap * pre_factor_b;

//                     Ma += micro_matrix_element_a;
//                     Mb += micro_matrix_element_b;
//                 }
//             }
//         }
//     }

//     return {Ma, Mb};
// }

}  // namespace bz_mesh