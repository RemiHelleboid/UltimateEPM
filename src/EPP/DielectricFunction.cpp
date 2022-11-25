/**
 * @file DielectricFunction.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-11-24
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "DielectricFunction.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "Hamiltonian.h"
#include "Material.h"
#include "omp.h"

namespace EmpiricalPseudopotential {

bool is_in_irreducible_wedge(const Vector3D<double>& k) {
    return (k.Z >= 0.0) && (k.Y >= k.Z) && (k.X >= k.Y) && (k.X <= 1.0) && (k.X + k.Y + k.Z <= 3.0 / 2.0);
}

bool is_in_first_BZ(const Vector3D<double>& k) {
    return (k.X >= 0.0) && (k.X <= 1.0) && (k.Y >= 0.0) && (k.Y <= 1.0) && (k.Z >= 0.0) && (k.Z <= 1.0);
}

DielectricFunction::DielectricFunction(const Material& material, const std::vector<Vector3D<int>>& basisVectors, const int nb_bands)
    : m_basisVectors(basisVectors),
      m_material(material),
      m_nb_bands(nb_bands) {}

void DielectricFunction::generate_k_points_random(std::size_t nb_points) {
    std::random_device               rd;
    std::mt19937                     gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    while (m_kpoints.size() < nb_points) {
        Vector3D<double> k(dis(gen), dis(gen), dis(gen));
        if (is_in_irreducible_wedge(k)) {
            m_kpoints.push_back(k);
        }
    }
}

void DielectricFunction::generate_k_points_grid(std::size_t Nx, std::size_t Ny, std::size_t Nz) {
    m_kpoints.clear();
    for (std::size_t i = 0; i < Nx; ++i) {
        for (std::size_t j = 0; j < Ny; ++j) {
            for (std::size_t k = 0; k < Nz; ++k) {
                Vector3D<double> k_vect(i / static_cast<double>(Nx), j / static_cast<double>(Ny), k / static_cast<double>(Nz));
                if (is_in_irreducible_wedge(k_vect)) {
                    m_kpoints.push_back(k_vect);
                }
            }
        }
    }
}

std::vector<double> DielectricFunction::compute_dielectric_function(const Vector3D<double>&    q_vect,
                                                                    const std::vector<double>& list_energies,
                                                                    double                     eta_smearing,
                                                                    int                        nb_threads) const {
    const int                     index_first_conduction_band = 4;
    double                        q_squared                   = pow(q_vect.Length(), 2);
    std::vector<Vector3D<double>> k_plus_q_vects(m_kpoints.size());
    std::transform(m_kpoints.begin(), m_kpoints.end(), k_plus_q_vects.begin(), [&q_vect](const Vector3D<double>& k) { return k + q_vect; });
    // Hamiltonian         hamiltonian_k(m_material, m_basisVectors);
    // Hamiltonian         hamiltonian_k_plus_q(m_material, m_basisVectors);
    std::vector<double>      iter_dielectric_function(m_kpoints.size());
    bool                     keep_eigenvectors = true;
    std::vector<Hamiltonian> hamiltonian_k_per_thread;
    std::vector<Hamiltonian> hamiltonian_k_plus_q_per_thread;

    for (int i = 0; i < nb_threads; i++) {
        hamiltonian_k_per_thread.push_back(Hamiltonian(m_material, m_basisVectors));
        hamiltonian_k_plus_q_per_thread.push_back(Hamiltonian(m_material, m_basisVectors));
    }

    std::vector<double> list_total_sum(list_energies.size());
#pragma omp parallel for schedule(dynamic) num_threads(nb_threads)
    for (std::size_t index_k = 0; index_k < m_kpoints.size(); ++index_k) {
        auto k_vect        = m_kpoints[index_k];
        auto k_plus_q_vect = k_plus_q_vects[index_k];
        int  thread_id     = omp_get_thread_num();
        hamiltonian_k_per_thread[thread_id].SetMatrix(k_vect);
        hamiltonian_k_plus_q_per_thread[thread_id].SetMatrix(k_plus_q_vect);
        hamiltonian_k_per_thread[thread_id].Diagonalize(keep_eigenvectors);
        hamiltonian_k_plus_q_per_thread[thread_id].Diagonalize(keep_eigenvectors);
        const auto&         eigenvalues_k         = hamiltonian_k_per_thread[thread_id].eigenvalues();
        const auto&         eigenvalues_k_plus_q  = hamiltonian_k_plus_q_per_thread[thread_id].eigenvalues();
        const auto&         eigenvectors_k        = hamiltonian_k_per_thread[thread_id].get_eigenvectors();
        const auto&         eigenvectors_k_plus_q = hamiltonian_k_plus_q_per_thread[thread_id].get_eigenvectors();
        std::vector<double> list_k_sum(list_energies.size());
        for (int idx_conduction_band = index_first_conduction_band; idx_conduction_band < m_nb_bands; ++idx_conduction_band) {
            for (int idx_valence_band = 0; idx_valence_band < index_first_conduction_band; ++idx_valence_band) {
                double overlap_integral =
                    pow(std::abs(eigenvectors_k_plus_q.col(idx_conduction_band).dot(eigenvectors_k.col(idx_valence_band))), 2);
                double delta_energy = eigenvalues_k_plus_q[idx_conduction_band] - eigenvalues_k[idx_valence_band];
                for (std::size_t index_energy = 0; index_energy < list_energies.size(); ++index_energy) {
                    double energy = list_energies[index_energy];
                    double factor_1 =
                        (delta_energy - energy) / ((delta_energy - energy) * (delta_energy - energy) + eta_smearing * eta_smearing);
                    double factor_2 =
                        (delta_energy + energy) / ((delta_energy + energy) * (delta_energy + energy) + eta_smearing * eta_smearing);
                    double total_factor = factor_1 + factor_2;
                    list_k_sum[index_energy] += overlap_integral * total_factor;
                }
            }
        }
        for (std::size_t index_energy = 0; index_energy < list_energies.size(); ++index_energy) {
            list_total_sum[index_energy] += list_k_sum[index_energy];
        }
        if (thread_id == 0) {
            std::cout << "\r"
                      << "Computing dielectric function: " << index_k << "/" << m_kpoints.size() << std::flush;
        }
    }
    for (std::size_t index_energy = 0; index_energy < list_energies.size(); ++index_energy) {
        list_total_sum[index_energy] /= m_kpoints.size();
    }
    std::vector<double> list_epsilon(list_energies.size());
    for (std::size_t index_energy = 0; index_energy < list_energies.size(); ++index_energy) {
        list_epsilon[index_energy] = 1.0 + (4.0 * M_PI / q_squared) * list_total_sum[index_energy];
    }
    return list_epsilon;
}

void DielectricFunction::export_kpoints(const std::string& filename) const {
    std::ofstream file(filename);
    file << "X,Y,Z" << std::endl;
    for (const auto& k : m_kpoints) {
        file << k.X << "," << k.Y << "," << k.Z << std::endl;
    }
}

}  // namespace EmpiricalPseudopotential