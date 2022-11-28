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
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "Hamiltonian.h"
#include "Material.h"

namespace EmpiricalPseudopotential {

bool is_in_irreducible_wedge(const Vector3D<double>& k) {
    return (k.Z >= 0.0) && (k.Y >= k.Z) && (k.X >= k.Y) && (k.X <= 1.0) && (k.X + k.Y + k.Z <= 3.0 / 2.0);
}

bool is_in_first_BZ(const Vector3D<double>& k, bool one_eighth = false) {
    bool cond_1      = fabs(k.X) <= 1.0 && fabs(k.Y) <= 1.0 && fabs(k.Z) <= 1.0;
    bool cond_2      = fabs(k.X) + fabs(k.Y) + fabs(k.Z) <= 3.0 / 2.0;
    bool cond_eighth = (k.X >= 0.0 && k.Y >= 0.0 && k.Z >= 0.0);
    return cond_1 && cond_2 && (one_eighth ? cond_eighth : true);
}

DielectricFunction::DielectricFunction(const Material& material, const std::vector<Vector3D<int>>& basisVectors, const int nb_bands)
    : m_basisVectors(basisVectors),
      m_material(material),
      m_nb_bands(nb_bands) {}

void DielectricFunction::generate_k_points_random(std::size_t nb_points) {
    std::random_device               rd;
    std::mt19937                     gen(rd());
    std::uniform_real_distribution<> dis(-1 - .0, 1.0);
    while (m_kpoints.size() < nb_points) {
        Vector3D<double> k(dis(gen), dis(gen), dis(gen));
        if (is_in_first_BZ(k)) {
            m_kpoints.push_back(k);
        }
    }
}

void DielectricFunction::generate_k_points_grid(std::size_t Nx, std::size_t Ny, std::size_t Nz, double shift, bool irreducible_wedge) {
    m_kpoints.clear();
    double min = -1.0 - shift;
    double max = 1.0 + shift;
    for (std::size_t i = 0; i < Nx; ++i) {
        for (std::size_t j = 0; j < Ny; ++j) {
            for (std::size_t k = 0; k < Nz; ++k) {
                Vector3D<double> k_vect(min + (max - min) * i / (Nx - 1),
                                        min + (max - min) * j / (Ny - 1),
                                        min + (max - min) * k / (Nz - 1));
                if (is_in_first_BZ(k_vect) && (!irreducible_wedge || is_in_irreducible_wedge(k_vect))) {
                    m_kpoints.push_back(k_vect);
                }
            }
        }
    }
}

/**
 * @brief Compute the energy and wavevector dependent dielectric function.
 * The formula used is the one from the paper "
 * 
 * @param eta_smearing 
 */
void DielectricFunction::compute_dielectric_function(double eta_smearing) {
    const int           index_first_conduction_band = 4;
    std::vector<double> iter_dielectric_function(m_kpoints.size());
    bool                keep_eigenvectors = true;

    std::size_t nb_kpoints = m_kpoints.size();
    m_eigenvalues_k.resize(nb_kpoints);
    m_eigenvectors_k.resize(nb_kpoints);
    auto        start = std::chrono::high_resolution_clock::now();
    Hamiltonian hamiltonian_k(m_material, m_basisVectors);
    Hamiltonian hamiltonian_k_plus_q(m_material, m_basisVectors);
    for (std::size_t index_q; index_q < m_qpoints.size(); ++index_q) {
        std::cout << "Computing dielectric function for q = " << m_qpoints[index_q] << std::endl;
        Vector3D<double>              q_vect    = m_qpoints[index_q];
        double                        q_squared = pow(q_vect.Length(), 2);
        std::vector<Vector3D<double>> k_plus_q_vects(m_kpoints.size());
        std::transform(m_kpoints.begin(), m_kpoints.end(), k_plus_q_vects.begin(), [&q_vect](const Vector3D<double>& k) {
            return k + q_vect;
        });
        std::vector<double> list_total_sum(m_energies.size());
        for (std::size_t index_k = 0; index_k < m_kpoints.size(); ++index_k) {
            // std::cout << "\rIndex k = " << index_k << std::flush;
            if (index_q == 0) {
                auto k_vect = m_kpoints[index_k];
                hamiltonian_k.SetMatrix(k_vect);
                hamiltonian_k.Diagonalize(keep_eigenvectors);

                m_eigenvalues_k[index_k]  = hamiltonian_k.eigenvalues();
                m_eigenvectors_k[index_k] = hamiltonian_k.get_eigenvectors();
                // Keep only firsts columns
                m_eigenvectors_k[index_k].conservativeResize(m_eigenvectors_k[index_k].size(), m_nb_bands);
            }
            auto k_plus_q_vect = k_plus_q_vects[index_k];
            hamiltonian_k_plus_q.SetMatrix(k_plus_q_vect);
            hamiltonian_k_plus_q.Diagonalize(keep_eigenvectors);
            const auto&         eigenvalues_k_plus_q  = hamiltonian_k_plus_q.eigenvalues();
            const auto&         eigenvectors_k_plus_q = hamiltonian_k_plus_q.get_eigenvectors();
            std::vector<double> list_k_sum(m_energies.size());
            for (int idx_conduction_band = index_first_conduction_band; idx_conduction_band < m_nb_bands; ++idx_conduction_band) {
                for (int idx_valence_band = 0; idx_valence_band < index_first_conduction_band; ++idx_valence_band) {
                    const auto& eigen_vect_k = m_eigenvectors_k[index_k].col(idx_valence_band);
                    double      overlap_integral =
                        pow(std::abs(eigenvectors_k_plus_q.col(idx_conduction_band).dot(m_eigenvectors_k[index_k].col(idx_valence_band))),
                            2);
                    double delta_energy = eigenvalues_k_plus_q[idx_conduction_band] - m_eigenvalues_k[index_k][idx_valence_band];
                    for (std::size_t index_energy = 0; index_energy < m_energies.size(); ++index_energy) {
                        double energy = m_energies[index_energy];
                        double factor_1 =
                            (delta_energy - energy) / ((delta_energy - energy) * (delta_energy - energy) + eta_smearing * eta_smearing);
                        double factor_2 =
                            (delta_energy + energy) / ((delta_energy + energy) * (delta_energy + energy) + eta_smearing * eta_smearing);
                        double total_factor = factor_1 + factor_2;
                        list_k_sum[index_energy] += overlap_integral * total_factor;
                    }
                }
            }
            if (m_qpoints.size()<=1) {
                // if there is only one q point in the list, we don't keep the eigenvectors, to save memory.
                m_eigenvectors_k[index_k].resize(1, 1);
            }
            for (std::size_t index_energy = 0; index_energy < m_energies.size(); ++index_energy) {
                list_total_sum[index_energy] += list_k_sum[index_energy];
            }
            // if (thread_id == 0) {
            //     std::cout << "\r"
            //               << "Computing dielectric function: " << index_k << "/" << m_kpoints.size() << std::flush;
            // }
        }
        for (std::size_t index_energy = 0; index_energy < m_energies.size(); ++index_energy) {
            list_total_sum[index_energy] /= m_kpoints.size();
        }
        std::vector<double> list_epsilon(m_energies.size());
        for (std::size_t index_energy = 0; index_energy < m_energies.size(); ++index_energy) {
            list_epsilon[index_energy] = 1.0 + (4.0 * M_PI / q_squared) * list_total_sum[index_energy];
        }
        auto end = std::chrono::high_resolution_clock::now();
        // std::cout << std::endl;
        std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0 << " s"
                  << std::endl;
        start = std::chrono::high_resolution_clock::now();
        m_dielectric_function.push_back(list_epsilon);
        this->export_dielectric_function_at_q("", index_q, true);
    }
}

void DielectricFunction::export_dielectric_function_at_q(const std::string& filename, std::size_t idx_q, bool name_auto) const {
    std::string outname;
    if (name_auto) {
        outname = "Q100/dielectric_function_q_" + std::to_string(m_qpoints[idx_q].X) + "_" + std::to_string(m_qpoints[idx_q].Y) + "_" +
                  std::to_string(m_qpoints[idx_q].Z) + ".csv";
    } else {
        outname = filename;
    }
    std::ofstream outfile(outname);
    outfile << "Energy (eV),Epsilon" << std::endl;
    for (std::size_t index_energy = 0; index_energy < m_energies.size(); ++index_energy) {
        outfile << m_energies[index_energy] << "," << m_dielectric_function[idx_q][index_energy] << std::endl;
    }
    outfile.close();
}

void DielectricFunction::export_kpoints(const std::string& filename) const {
    std::ofstream file(filename);
    file << "X,Y,Z" << std::endl;
    for (const auto& k : m_kpoints) {
        file << k.X << "," << k.Y << "," << k.Z << std::endl;
    }
}

}  // namespace EmpiricalPseudopotential