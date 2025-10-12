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

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "Hamiltonian.h"
#include "Material.h"

namespace uepm::pseudopotential {

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
    std::uniform_real_distribution<> dis(-1, 1.0);
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
    for (std::size_t i = 0; i < Nx - 1; ++i) {
        for (std::size_t j = 0; j < Ny - 1; ++j) {
            for (std::size_t k = 0; k < Nz - 1; ++k) {
                Vector3D<double> k_vect(min + (max - min) * i / static_cast<double>(Nx - 1),
                                        min + (max - min) * j / static_cast<double>(Ny - 1),
                                        min + (max - min) * k / static_cast<double>(Nz - 1));
                if (is_in_first_BZ(k_vect) && (!irreducible_wedge || is_in_irreducible_wedge(k_vect))) {
                    m_kpoints.push_back(k_vect);
                }
            }
        }
    }
}

/**
 * @brief Compute the energy and wave vector dependent dielectric function.
 * The formula used is the one from the paper "
 *
 * @param eta_smearing
 */
void DielectricFunction::compute_dielectric_function(double eta_smearing, int mpi_rank) {
    const int           index_first_conduction_band = 4;
    std::vector<double> iter_dielectric_function(m_kpoints.size());
    const bool          keep_eigenvectors = true;

    std::size_t nb_kpoints = m_kpoints.size();
    m_eigenvalues_k.resize(nb_kpoints);
    m_eigenvectors_k.resize(nb_kpoints);
    auto        start = std::chrono::high_resolution_clock::now();
    Hamiltonian hamiltonian_k(m_material, m_basisVectors);
    Hamiltonian hamiltonian_k_plus_q(m_material, m_basisVectors);
    for (std::size_t index_q = 0; index_q < m_qpoints.size(); ++index_q) {
        Vector3D<double>              q_vect = m_qpoints[index_q];
        if (q_vect.Length() <= 1e-15) {
            q_vect = Vector3D<double>(1e-15, 1e-15, 1e-15);
        }
        std::vector<Vector3D<double>> k_plus_q_vects(m_kpoints.size());
        std::transform(m_kpoints.begin(), m_kpoints.end(), k_plus_q_vects.begin(), [&q_vect](const Vector3D<double>& k) {
            return k + q_vect;
        });
        std::vector<double> list_total_sum(m_energies.size());
        for (std::size_t index_k = m_offset_k_index; index_k < m_offset_k_index + m_nb_kpoints; ++index_k) {
            if (index_q == 0) {
                auto k_vect = m_kpoints[index_k];
                hamiltonian_k.SetMatrix(k_vect, m_nonlocal_epm);
                hamiltonian_k.Diagonalize(keep_eigenvectors);
                m_eigenvalues_k[index_k]  = hamiltonian_k.eigenvalues();
                m_eigenvectors_k[index_k] = hamiltonian_k.get_eigenvectors();
                // Keep only firsts columns
                auto nb_rows = m_eigenvectors_k[index_k].rows();
                m_eigenvectors_k[index_k].conservativeResize(nb_rows, m_nb_bands);
            }
            auto k_plus_q_vect = k_plus_q_vects[index_k];
            hamiltonian_k_plus_q.SetMatrix(k_plus_q_vect, m_nonlocal_epm);
            hamiltonian_k_plus_q.Diagonalize(keep_eigenvectors);
            const auto&         eigenvalues_k_plus_q  = hamiltonian_k_plus_q.eigenvalues();
            const auto&         eigenvectors_k_plus_q = hamiltonian_k_plus_q.get_eigenvectors();
            std::vector<double> list_k_sum(m_energies.size());
            for (int idx_conduction_band = index_first_conduction_band; idx_conduction_band < m_nb_bands; ++idx_conduction_band) {
                for (int idx_valence_band = 0; idx_valence_band < index_first_conduction_band; ++idx_valence_band) {
                    double overlap_integral = pow(
                        std::abs(
                            eigenvectors_k_plus_q.col(idx_conduction_band).adjoint().dot(m_eigenvectors_k[index_k].col(idx_valence_band))),
                        2);
                    double delta_energy = (eigenvalues_k_plus_q[idx_conduction_band]) - m_eigenvalues_k[index_k][idx_valence_band];
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
            if (m_qpoints.size() <= 1) {
                // if there is only one q point in the list, we don't keep the eigenvectors, to save memory.
                m_eigenvectors_k[index_k].resize(1, 1);
                m_eigenvalues_k[index_k].resize(1);
            }
            for (std::size_t index_energy = 0; index_energy < m_energies.size(); ++index_energy) {
                list_total_sum[index_energy] += list_k_sum[index_energy];
            }
        }
        std::vector<double> list_epsilon(m_energies.size());
        for (std::size_t index_energy = 0; index_energy < m_energies.size(); ++index_energy) {
            list_epsilon[index_energy] = list_total_sum[index_energy];
        }
        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
        if (mpi_rank == 0) {
            std::cout << "Computed dielectric function for q = " << m_qpoints[index_q] << " -> " << index_q + 1 << "/" << m_qpoints.size()
                      << " in " << elapsed << " s" << std::endl;
        }
        start = std::chrono::high_resolution_clock::now();
        m_dielectric_function_real.push_back(list_epsilon);
    }
}

DielectricFunction DielectricFunction::merge_results(DielectricFunction                                  RootDielectricFunction,
                                                     const std::vector<std::vector<std::vector<double>>>& dielectric_function_results,
                                                     std::vector<int>                                    nb_kpoints_per_instance) {
    std::vector<std::vector<double>> total_dielectric_function;
    if (dielectric_function_results.size() == 0) {
        throw std::runtime_error("No results to merge");
    }
    if (dielectric_function_results.size() != nb_kpoints_per_instance.size()) {
        throw std::runtime_error("Number of results and number of k-points per instance do not match");
    }
    std::size_t total_number_kpoints = 0;
    for (std::size_t index_instance = 0; index_instance < nb_kpoints_per_instance.size(); ++index_instance) {
        total_number_kpoints += nb_kpoints_per_instance[index_instance];
    }
    std::cout << "Number total kpoint to the merge : " << total_number_kpoints << std::endl;
    // Add-up the k-contributions.
    for (std::size_t index_instance = 0; index_instance < dielectric_function_results.size(); ++index_instance) {
        // Add the contributions to epsilon for each q-point and energy.
        for (std::size_t index_q = 0; index_q < dielectric_function_results[index_instance].size(); ++index_q) {
            if (index_instance == 0) {
                total_dielectric_function.push_back(dielectric_function_results[index_instance][index_q]);
            } else {
                for (std::size_t index_energy = 0; index_energy < dielectric_function_results[index_instance][index_q].size();
                     ++index_energy) {
                    total_dielectric_function[index_q][index_energy] += dielectric_function_results[index_instance][index_q][index_energy];
                }
            }
        }
    }
    // Renormalization
    double       renormalization = 1.0 / static_cast<double>(total_number_kpoints);
    std::cout << "Renormalization: " << renormalization << std::endl;
    for (std::size_t index_q = 0; index_q < total_dielectric_function.size(); ++index_q) {
        Vector3D<double> q_vect    = RootDielectricFunction.m_qpoints[index_q];
        double           q_squared = pow(q_vect.Length(), 2);
        // TO IMPROVE: This is a hack to avoid division by zero.
        if (q_squared == 0.0) {
            q_squared = 1e-14;
        }
        for (std::size_t index_energy = 0; index_energy < total_dielectric_function[index_q].size(); ++index_energy) {
            total_dielectric_function[index_q][index_energy] *= renormalization;
            total_dielectric_function[index_q][index_energy] =
                1.0 + (2.0 * M_PI / q_squared) * total_dielectric_function[index_q][index_energy];
        }
    }

    DielectricFunction dielectric_function         = RootDielectricFunction;
    dielectric_function.m_dielectric_function_real = total_dielectric_function;
    return dielectric_function;
}

Eigen::MatrixXd create_kramers_matrix(std::size_t N) {
    Eigen::MatrixXd kramers_matrix(N, N);
    for (std::size_t idx_line = 0; idx_line < N; ++idx_line) {
        for (std::size_t idx_col = 0; idx_col < N; ++idx_col) {
            if (idx_line == idx_col) {
                kramers_matrix(idx_line, idx_col) = 0.0;
            } else {
                double a = double(idx_line) / double(static_cast<long>(idx_col * idx_col) - static_cast<long>(idx_line * idx_line));
                kramers_matrix(idx_line, idx_col) = a;
            }
        }
    }
    kramers_matrix *= -2.0 / M_PI;
    return kramers_matrix;
}

void DielectricFunction::apply_kramers_kronig() {
    std::cout << "Applying Kramers-Kronig" << std::endl;
    std::fstream kramers_matrix_file;
    kramers_matrix_file.open("Bkramers_matrix.dat", std::ios::out);
    Eigen::MatrixXd kramers_matrix2 = create_kramers_matrix(6);
    kramers_matrix_file << kramers_matrix2 << std::endl;
    kramers_matrix_file.close();
    std::cout << "Kramers matrix created" << std::endl;
    m_dielectric_function_imag.clear();
    m_dielectric_function_imag.resize(m_dielectric_function_real.size());
    Eigen::MatrixXd kramers_matrix = create_kramers_matrix(m_energies.size());
    for (std::size_t idx_q = 0; idx_q < m_qpoints.size(); ++idx_q) {
        Eigen::VectorXd epsilon(m_energies.size());
        for (std::size_t idx_energy = 0; idx_energy < m_energies.size(); ++idx_energy) {
            epsilon(idx_energy) = m_dielectric_function_real[idx_q][idx_energy] - 1.0;
        }
        Eigen::VectorXd epsilon_imag = kramers_matrix * epsilon;
        m_dielectric_function_imag[idx_q].resize(m_energies.size());
        for (std::size_t idx_energy = 0; idx_energy < m_energies.size(); ++idx_energy) {
            m_dielectric_function_imag[idx_q][idx_energy] = epsilon_imag(idx_energy);
        }
    }
}

void DielectricFunction::export_dielectric_function_at_q(const std::string& filename, std::size_t idx_q, bool name_auto) const {
    std::string outname;
    if (name_auto) {
        // outname = m_export_prefix + '_' + std::to_string(idx_q) + '_' + std::to_string(m_qpoints[idx_q].X) + "_" +
        // std::to_string(m_qpoints[idx_q].Y) + "_" +
        //           std::to_string(m_qpoints[idx_q].Z) + ".csv";
        outname = fmt::format("{}_{:05}_{:.6f}_{:.6f}_{:.6f}.csv",
                              m_export_prefix,
                              idx_q,
                              m_qpoints[idx_q].X,
                              m_qpoints[idx_q].Y,
                              m_qpoints[idx_q].Z);
    } else {
        outname = filename;
    }
    std::ofstream outfile(outname);
    std::cout << m_energies[0] << " " << m_dielectric_function_real[idx_q][0] << std::endl;
    outfile << "Energy (eV),EpsilonReal,EpsilonImaginary" << std::endl;
    for (std::size_t index_energy = 0; index_energy < m_energies.size(); ++index_energy) {
        outfile << m_energies[index_energy] << "," << m_dielectric_function_real[idx_q][index_energy] << ","
                << m_dielectric_function_imag[idx_q][index_energy] << std::endl;
    }
    outfile.close();
}

void DielectricFunction::export_dielectric_function(const std::string& filename, bool name_auto) const {
    for (std::size_t index_q = 0; index_q < m_qpoints.size(); ++index_q) {
        export_dielectric_function_at_q(filename, index_q, name_auto);
    }
}

void DielectricFunction::export_kpoints(const std::string& filename) const {
    std::ofstream file(filename);
    file << "X,Y,Z" << std::endl;
    for (const auto& k : m_kpoints) {
        file << k.X << "," << k.Y << "," << k.Z << std::endl;
    }
}

}  // namespace uepm::pseudopotential