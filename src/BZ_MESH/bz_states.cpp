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
        if (omp_get_thread_num() == 0) {
            std::cout << "\rComputing eigenstates for k = " << idx_k << "/" << m_list_vertices.size() << std::flush;
        }
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
        // std::cout << "\r Nb rows = " << m_eigenvectors_k[idx_k].rows() << " Nb cols = " << m_eigenvectors_k[idx_k].cols() << std::endl;
    }
    std::cout << std::endl;
}

void BZ_States::compute_shifted_eigenstates(const Vector3D<double>& q_shift, int nb_threads) {
    m_q_shift                       = q_shift;
    double     normalization_factor = 2.0 * M_PI / m_material.get_lattice_constant_meter();
    const bool m_nonlocal_epm       = false;
    const bool keep_eigenvectors    = true;
    m_eigenvalues_k_plus_q.resize(m_list_vertices.size());
    m_eigenvectors_k_plus_q.resize(m_list_vertices.size());
    std::vector<EmpiricalPseudopotential::Hamiltonian> hamiltonian_per_thread;
    for (int i = 0; i < nb_threads; i++) {
        hamiltonian_per_thread.push_back(EmpiricalPseudopotential::Hamiltonian(m_material, m_basisVectors));
    }
    for (std::size_t idx_k = 0; idx_k < m_list_vertices.size(); ++idx_k) {
        std::cout << "\rComputing eigenstates for k+q = " << idx_k << "/" << m_list_vertices.size() << std::flush;
        auto k_point = Vector3D<double>(m_list_vertices[idx_k].get_position().x(),
                                        m_list_vertices[idx_k].get_position().y(),
                                        m_list_vertices[idx_k].get_position().z());
        k_point      = k_point * 1.0 / normalization_factor;
        k_point += q_shift;
        auto idx_thread = omp_get_thread_num();
        hamiltonian_per_thread[idx_thread].SetMatrix(k_point, m_nonlocal_epm);
        hamiltonian_per_thread[idx_thread].Diagonalize(keep_eigenvectors);
        m_eigenvalues_k_plus_q[idx_k]  = hamiltonian_per_thread[idx_thread].eigenvalues();
        m_eigenvectors_k_plus_q[idx_k] = hamiltonian_per_thread[idx_thread].get_eigenvectors();
        auto nb_rows                   = m_eigenvectors_k[idx_k].rows();
        m_eigenvectors_k_plus_q[idx_k].conservativeResize(nb_rows, m_nb_bands);
    }
    std::cout << std::endl;
}

/**
 * @brief Compute the dielectric function for a given list of energies and a given smearing.
 * The integration is performed by summing the contribution of each tetrahedron to the dielectric function.
 *
 * @param energies
 * @param eta_smearing
 * @param nb_threads
 */
void BZ_States::compute_dielectric_function(const std::vector<double>& list_energies, double eta_smearing, int nb_threads) {
    m_list_energies                         = list_energies;
    const int   index_first_conduction_band = 4;
    std::size_t nb_tetra                    = m_list_tetrahedra.size();
    m_dielectric_function_real.resize(list_energies.size());
    constexpr double one_fourth = 1.0 / 4.0;

    std::vector<double> dielectric_function_real_at_energies(list_energies.size(), 0.0);
    double              total_volume = 0.0;
    for (std::size_t idx_tetra = 0; idx_tetra < nb_tetra; ++idx_tetra) {
        std::cout << "\rComputing dielectric function for tetrahedron " << idx_tetra << "/" << nb_tetra << std::flush;
        std::array<std::size_t, 4>    list_idx_vertices = m_list_tetrahedra[idx_tetra].get_list_indices_vertices();
        const std::array<Vertex*, 4>& list_vertices     = m_list_tetrahedra[idx_tetra].get_list_vertices();
        double                        volume_tetra      = std::fabs(m_list_tetrahedra[idx_tetra].compute_signed_volume());
        total_volume += volume_tetra;
        // std::cout << "Volume tetra = " << volume_tetra << std::endl;
        std::vector<double> sum_dielectric_function_real_tetra_at_energies(list_energies.size(), 0.0);
        // Loop over the vertices of the tetrahedron
        for (std::size_t idx_vertex = 0; idx_vertex < 4; ++idx_vertex) {
            std::size_t index_k = list_idx_vertices[idx_vertex];
            for (int idx_conduction_band = index_first_conduction_band; idx_conduction_band < m_nb_bands; ++idx_conduction_band) {
                for (int idx_valence_band = 0; idx_valence_band < index_first_conduction_band; ++idx_valence_band) {
                    double overlap_integral = pow(std::fabs(m_eigenvectors_k_plus_q[index_k]
                                                                .col(idx_conduction_band)
                                                                .adjoint()
                                                                .dot(m_eigenvectors_k[index_k].col(idx_valence_band))),
                                                  2);
                    double delta_energy = m_eigenvalues_k_plus_q[index_k][idx_conduction_band] - m_eigenvalues_k[index_k][idx_valence_band];
                    for (std::size_t index_energy = 0; index_energy < list_energies.size(); ++index_energy) {
                        double energy = list_energies[index_energy];
                        double factor_1 =
                            (delta_energy - energy) / ((delta_energy - energy) * (delta_energy - energy) + eta_smearing * eta_smearing);
                        double factor_2 =
                            (delta_energy + energy) / ((delta_energy + energy) * (delta_energy + energy) + eta_smearing * eta_smearing);
                        double total_factor = factor_1 + factor_2;
                        sum_dielectric_function_real_tetra_at_energies[index_energy] += overlap_integral * total_factor;
                    }
                }
            }
        }
        for (std::size_t index_energy = 0; index_energy < list_energies.size(); ++index_energy) {
            sum_dielectric_function_real_tetra_at_energies[index_energy] *= volume_tetra * one_fourth;
            dielectric_function_real_at_energies[index_energy] += sum_dielectric_function_real_tetra_at_energies[index_energy];
        }
    }
    std::cout << "\n";
    std::cout << "Total volume: " << total_volume << std::endl;

    double q_squared  = m_q_shift.Length() * m_q_shift.Length();
    double pre_factor = 2.0 * M_PI / q_squared;
    for (std::size_t index_energy = 0; index_energy < list_energies.size(); ++index_energy) {
        m_dielectric_function_real[index_energy] = 1.0 + pre_factor * dielectric_function_real_at_energies[index_energy] / total_volume;
    }
    std::cout << "EPS[0] = " << m_dielectric_function_real[0] << std::endl;
}

// Export the dielectric function to a file in the format (energy, dielectric function) (csv format).
void BZ_States::export_dielectric_function(const std::string& prefix) const {
    std::ofstream dielectric_function_file(prefix + "_dielectric_function.csv");
    dielectric_function_file << "Energy (eV), Dielectric function" << std::endl;
    for (std::size_t index_energy = 0; index_energy < m_dielectric_function_real.size(); ++index_energy) {
        dielectric_function_file << m_list_energies[index_energy] << ", " << m_dielectric_function_real[index_energy] << std::endl;
    }
    dielectric_function_file.close();
}

void BZ_States::export_full_eigenstates() const {
    std::filesystem::remove_all("eigenstates");
    std::filesystem::create_directory("eigenstates");
    std::filesystem::create_directory("eigenstates/eigenvectors");
    std::filesystem::create_directory("eigenstates/eigenvalues");
    for (std::size_t idx_k = 0; idx_k < m_list_vertices.size(); ++idx_k) {
        std::ofstream eigenvalues_file("eigenstates/eigenvalues/eigenvalues_" + std::to_string(idx_k) + ".txt");
        eigenvalues_file << m_eigenvalues_k[idx_k].transpose() << std::endl;
        eigenvalues_file.close();
        std::ofstream shiftedeigenvalues_file("eigenstates/eigenvalues/shiftedeigenvalues_" + std::to_string(idx_k) + ".txt");
        shiftedeigenvalues_file << m_eigenvalues_k_plus_q[idx_k].transpose() << std::endl;
        shiftedeigenvalues_file.close();

        std::ofstream eigenvectors_file("eigenstates/eigenvectors/eigenvectors_" + std::to_string(idx_k) + ".txt");
        eigenvectors_file << m_eigenvectors_k[idx_k] << std::endl;
        eigenvectors_file.close();
        std::ofstream shiftedeigenvectors_file("eigenstates/eigenvectors/shiftedeigenvectors_" + std::to_string(idx_k) + ".txt");
        shiftedeigenvectors_file << m_eigenvectors_k_plus_q[idx_k] << std::endl;
        shiftedeigenvectors_file.close();
    }
}

// void BZ_States::populate_vtx_dielectric_function(const std::vector<double>& energies, double eta_smearing);

}  // namespace bz_mesh