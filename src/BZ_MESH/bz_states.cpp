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

#include "bz_states.hpp"

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>

#include "Hamiltonian.h"
#include "bz_mesh.hpp"
#include "omp.h"

namespace bz_mesh {

void BZ_States::compute_eigenstates(int nb_threads) {
    double     normalization_factor = 2.0 * M_PI / m_material.get_lattice_constant_meter();
    const bool m_nonlocal_epm       = false;
    const bool keep_eigenvectors    = true;
    m_eigenvalues_k.resize(m_list_vertices.size());
    m_eigenvectors_k.resize(m_list_vertices.size());
    std::vector<EmpiricalPseudopotential::Hamiltonian> hamiltonian_per_thread;
    for (int i = 0; i < nb_threads; i++) {
        hamiltonian_per_thread.push_back(EmpiricalPseudopotential::Hamiltonian(m_material, m_basisVectors));
    }
#pragma omp parallel for schedule(dynamic) num_threads(nb_threads)
    for (std::size_t idx_k = 0; idx_k < m_list_vertices.size(); ++idx_k) {
        std::cout << "\rComputing eigenstates for k = " << idx_k << "/" << m_list_vertices.size() << std::flush;
        auto k_point    = Vector3D<double>(m_list_vertices[idx_k].get_position().x(),
                                        m_list_vertices[idx_k].get_position().y(),
                                        m_list_vertices[idx_k].get_position().z());
        k_point         = k_point * 1.0 / normalization_factor;
        auto idx_thread = omp_get_thread_num();
        hamiltonian_per_thread[idx_thread].SetMatrix(k_point, m_nonlocal_epm);
        hamiltonian_per_thread[idx_thread].Diagonalize(keep_eigenvectors);
        m_eigenvalues_k[idx_k]  = hamiltonian_per_thread[idx_thread].eigenvalues();
        m_eigenvectors_k[idx_k] = hamiltonian_per_thread[idx_thread].get_eigenvectors();
        auto nb_rows            = m_eigenvectors_k[idx_k].rows();
        m_eigenvectors_k[idx_k].conservativeResize(nb_rows, m_nb_bands);
    }
}

// double BZ_States::compute_direct_impact_ionization_matrix_element(int                     idx_n1,
//                                                                   int                     idx_n1_prime,
//                                                                   int                     idx_n2,
//                                                                   int                     idx_n2_prime,
//                                                                   const Vector3D<double>& k1,
//                                                                   const Vector3D<double>& k2,
//                                                                   const Vector3D<double>& k1_prime,
//                                                                   const Vector3D<double>& k2_prime) const {
//     double     matrix_element       = 0.0;
//     double     normalization_factor = 2.0 * M_PI / m_material.get_lattice_constant_meter();
//     const bool m_nonlocal_epm       = false;
//     const bool keep_eigenvectors    = true;
// }

void BZ_States::export_full_eigenstates() const {
    std::filesystem::remove_all("eigenstates");
    std::filesystem::create_directory("eigenstates");
    std::filesystem::create_directory("eigenstates/eigenvectors");
    std::filesystem::create_directory("eigenstates/eigenvalues");
    for (std::size_t idx_k = 0; idx_k < m_list_vertices.size(); ++idx_k) {
        std::ofstream eigenvalues_file("eigenstates/eigenvalues/eigenvalues_" + std::to_string(idx_k) + ".txt");
        eigenvalues_file << m_eigenvalues_k[idx_k].transpose() << std::endl;
        eigenvalues_file.close();

        std::ofstream eigenvectors_file("eigenstates/eigenvectors/eigenvectors_" + std::to_string(idx_k) + ".txt");
        eigenvectors_file << m_eigenvectors_k[idx_k] << std::endl;
        eigenvectors_file.close();
    }
}

}  // namespace bz_mesh