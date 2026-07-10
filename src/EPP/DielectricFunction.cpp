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
#include <stdexcept>
#include <vector>

#include "Hamiltonian.h"
#include "epm_material.hpp"
#include "physical_constants.hpp"

namespace uepm::pseudopotential {

bool is_in_irreducible_wedge(const Vector3D<double>& k) {
    return (k.Z >= 0.0) && (k.Y >= k.Z) && (k.X >= k.Y) && (k.X <= 1.0) && (k.X + k.Y + k.Z <= 3.0 / 2.0);
}

bool is_in_first_BZ(const Vector3D<double>& k, bool positive_octant = false) {
    bool cond_1      = fabs(k.X) <= 1.0 && fabs(k.Y) <= 1.0 && fabs(k.Z) <= 1.0;
    bool cond_2      = fabs(k.X) + fabs(k.Y) + fabs(k.Z) <= 3.0 / 2.0;
    bool cond_octant = (k.X >= 0.0 && k.Y >= 0.0 && k.Z >= 0.0);
    return cond_1 && cond_2 && (positive_octant ? cond_octant : true);
}

DielectricFunction::DielectricFunction(const epm_material&               material,
                                       const std::vector<Vector3D<int>>& basisVectors,
                                       const int                         nb_bands)
    : m_basisVectors(basisVectors),
      m_material(material),
      m_nb_bands(nb_bands) {}

void DielectricFunction::set_optical_direction(const Vector3D<double>& direction) {
    const double norm = direction.Length();
    if (!(norm > 0.0) || !std::isfinite(norm)) {
        throw std::invalid_argument("DielectricFunction::set_optical_direction requires a finite non-zero direction.");
    }
    m_optical_direction = direction / norm;
}

void DielectricFunction::set_optical_derivative_step(double step) {
    if (!(step > 0.0) || !std::isfinite(step)) {
        throw std::invalid_argument("DielectricFunction::set_optical_derivative_step requires a positive finite step.");
    }
    m_optical_derivative_step = step;
}

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

void DielectricFunction::generate_k_points_grid(std::size_t              Nx,
                                                std::size_t              Ny,
                                                std::size_t              Nz,
                                                double                   shift,
                                                DielectricKPointSampling sampling) {
    if (Nx == 0 || Ny == 0 || Nz == 0) {
        throw std::invalid_argument("DielectricFunction::generate_k_points_grid requires non-zero grid dimensions.");
    }
    m_kpoints.clear();
    const bool       positive_octant   = sampling == DielectricKPointSampling::q100_octant;
    const bool       irreducible_wedge = sampling == DielectricKPointSampling::fcc_irreducible_wedge;
    constexpr double min               = -1.0;
    constexpr double max               = 1.0;
    for (std::size_t i = 0; i < Nx; ++i) {
        for (std::size_t j = 0; j < Ny; ++j) {
            for (std::size_t k = 0; k < Nz; ++k) {
                Vector3D<double> k_vect(
                    min + (max - min) * (static_cast<double>(i) + 0.5) / static_cast<double>(Nx) + shift,
                    min + (max - min) * (static_cast<double>(j) + 0.5) / static_cast<double>(Ny) + shift,
                    min + (max - min) * (static_cast<double>(k) + 0.5) / static_cast<double>(Nz) + shift);
                if (is_in_first_BZ(k_vect, positive_octant) &&
                    (!irreducible_wedge || is_in_irreducible_wedge(k_vect))) {
                    m_kpoints.push_back(k_vect);
                }
            }
        }
    }
}

void DielectricFunction::generate_k_points_grid(std::size_t Nx,
                                                std::size_t Ny,
                                                std::size_t Nz,
                                                double      shift,
                                                bool        irreducible_wedge) {
    generate_k_points_grid(
        Nx,
        Ny,
        Nz,
        shift,
        irreducible_wedge ? DielectricKPointSampling::fcc_irreducible_wedge : DielectricKPointSampling::full_bz);
}

std::complex<double> local_optical_dH_element(const epm_material&               material,
                                              const std::vector<Vector3D<int>>& basis_vectors,
                                              const Vector3D<double>&           k,
                                              const Vector3D<double>&           direction,
                                              const Eigen::MatrixXcd&           eigenvectors,
                                              int                               idx_conduction_band,
                                              int                               idx_valence_band) {
    const Eigen::Index basis_size    = static_cast<Eigen::Index>(basis_vectors.size());
    const double       two_pi_over_a = 2.0 * M_PI / material.get_lattice_constant_meter();
    constexpr double   kinetic_prefactor =
        (uepm::constants::h_bar * uepm::constants::h_bar) / (2.0 * uepm::constants::m_e * uepm::constants::q_e);
    const double d_kinetic_prefactor = 2.0 * kinetic_prefactor * two_pi_over_a * two_pi_over_a;

    const Eigen::VectorXcd conduction = eigenvectors.col(idx_conduction_band);
    const Eigen::VectorXcd valence    = eigenvectors.col(idx_valence_band);
    std::complex<double>   element{0.0, 0.0};

    for (Eigen::Index i = 0; i < basis_size; ++i) {
        const Vector3D<double> k_plus_g = k + basis_vectors[static_cast<std::size_t>(i)];
        const double           dH_diag  = d_kinetic_prefactor * (k_plus_g * direction);
        element += std::conj(conduction(i)) * dH_diag * valence(i);
    }

    if (conduction.size() == 2 * basis_size) {
        for (Eigen::Index i = 0; i < basis_size; ++i) {
            const Vector3D<double> k_plus_g = k + basis_vectors[static_cast<std::size_t>(i)];
            const double           dH_diag  = d_kinetic_prefactor * (k_plus_g * direction);
            element += std::conj(conduction(i + basis_size)) * dH_diag * valence(i + basis_size);
        }
    }

    return element;
}

Eigen::MatrixXcd finite_difference_dH(const epm_material&               material,
                                      const std::vector<Vector3D<int>>& basis_vectors,
                                      const Vector3D<double>&           k,
                                      const Vector3D<double>&           direction,
                                      double                            step,
                                      bool                              nonlocal_epm) {
    Hamiltonian hamiltonian_plus(material, basis_vectors);
    Hamiltonian hamiltonian_minus(material, basis_vectors);
    hamiltonian_plus.SetMatrix(k + direction * step, nonlocal_epm);
    hamiltonian_minus.SetMatrix(k - direction * step, nonlocal_epm);
    return (hamiltonian_plus.get_matrix() - hamiltonian_minus.get_matrix()) / (2.0 * step);
}

double optical_transition_strength(const epm_material&               material,
                                   const std::vector<Vector3D<int>>& basis_vectors,
                                   const Vector3D<double>&           k,
                                   const Vector3D<double>&           direction,
                                   const Eigen::MatrixXcd&           eigenvectors,
                                   const Eigen::MatrixXcd*           finite_difference_derivative,
                                   int                               idx_conduction_band,
                                   int                               idx_valence_band,
                                   double                            delta_energy) {
    if (!(delta_energy > 0.0) || !std::isfinite(delta_energy)) {
        return 0.0;
    }

    std::complex<double> dH_element{0.0, 0.0};
    if (finite_difference_derivative != nullptr) {
        dH_element = eigenvectors.col(idx_conduction_band)
                         .dot((*finite_difference_derivative) * eigenvectors.col(idx_valence_band));
    } else {
        dH_element = local_optical_dH_element(material,
                                              basis_vectors,
                                              k,
                                              direction,
                                              eigenvectors,
                                              idx_conduction_band,
                                              idx_valence_band);
    }
    return std::norm(dH_element) / (delta_energy * delta_energy);
}

/**
 * @brief Compute the energy and wave vector dependent dielectric function.
 * The formula used is the one from the paper "
 *
 * @param eta_smearing
 */
void DielectricFunction::compute_dielectric_function(double eta_smearing, int mpi_rank) {
    const int index_first_conduction_band = 4;
    if (m_energies.empty()) {
        throw std::invalid_argument("DielectricFunction::compute_dielectric_function requires an energy grid.");
    }
    if (!(eta_smearing > 0.0) || !std::isfinite(eta_smearing)) {
        throw std::invalid_argument("DielectricFunction::compute_dielectric_function requires a positive finite eta.");
    }
    if (m_qpoints.empty()) {
        throw std::invalid_argument("DielectricFunction::compute_dielectric_function requires at least one q-point.");
    }
    if (m_nb_bands <= index_first_conduction_band) {
        throw std::invalid_argument("DielectricFunction::compute_dielectric_function needs conduction bands.");
    }
    const bool keep_eigenvectors = true;

    std::size_t nb_kpoints = m_kpoints.size();
    if (m_offset_k_index + m_nb_kpoints > nb_kpoints) {
        throw std::out_of_range("DielectricFunction::compute_dielectric_function k-point range is invalid.");
    }
    m_dielectric_function_real.clear();
    m_dielectric_function_imag.clear();
    m_eigenvalues_k.resize(nb_kpoints);
    m_eigenvectors_k.resize(nb_kpoints);
    auto        start         = std::chrono::high_resolution_clock::now();
    const bool  optical_limit = m_response_mode == DielectricResponseMode::optical_limit;
    Hamiltonian hamiltonian_k(m_material, m_basisVectors);
    Hamiltonian hamiltonian_k_plus_q(m_material, m_basisVectors);
    for (std::size_t index_q = 0; index_q < m_qpoints.size(); ++index_q) {
        Vector3D<double> q_vect            = m_qpoints[index_q];
        Vector3D<double> optical_direction = m_optical_direction;
        if (optical_limit && q_vect.Length() > 0.0) {
            optical_direction = q_vect / q_vect.Length();
        }
        if (!optical_limit && (!(q_vect.Length() > 0.0) || !std::isfinite(q_vect.Length()))) {
            throw std::invalid_argument("DielectricFunction::compute_dielectric_function requires finite non-zero q.");
        }
        std::vector<Vector3D<double>> k_plus_q_vects(m_kpoints.size());
        if (!optical_limit) {
            std::transform(m_kpoints.begin(),
                           m_kpoints.end(),
                           k_plus_q_vects.begin(),
                           [&q_vect](const Vector3D<double>& k) { return k + q_vect; });
        }
        std::vector<double> list_total_sum_real(m_energies.size());
        std::vector<double> list_total_sum_imag(m_energies.size());
        for (std::size_t index_k = m_offset_k_index; index_k < m_offset_k_index + m_nb_kpoints; ++index_k) {
            const auto k_vect = m_kpoints[index_k];
            if (index_q == 0) {
                hamiltonian_k.SetMatrix(k_vect, m_nonlocal_epm);
                hamiltonian_k.Diagonalize(keep_eigenvectors);
                m_eigenvalues_k[index_k]  = hamiltonian_k.eigenvalues();
                m_eigenvectors_k[index_k] = hamiltonian_k.get_eigenvectors();
                // Keep only firsts columns
                auto nb_rows = m_eigenvectors_k[index_k].rows();
                m_eigenvectors_k[index_k].conservativeResize(nb_rows, m_nb_bands);
            }
            Eigen::VectorXd  eigenvalues_k_plus_q;
            Eigen::MatrixXcd eigenvectors_k_plus_q;
            if (!optical_limit) {
                auto k_plus_q_vect = k_plus_q_vects[index_k];
                hamiltonian_k_plus_q.SetMatrix(k_plus_q_vect, m_nonlocal_epm);
                hamiltonian_k_plus_q.Diagonalize(keep_eigenvectors);
                eigenvalues_k_plus_q  = hamiltonian_k_plus_q.eigenvalues();
                eigenvectors_k_plus_q = hamiltonian_k_plus_q.get_eigenvectors();
            }
            Eigen::MatrixXcd finite_difference_derivative;
            const bool       use_finite_difference_derivative = optical_limit && m_nonlocal_epm;
            if (use_finite_difference_derivative) {
                finite_difference_derivative = finite_difference_dH(m_material,
                                                                    m_basisVectors,
                                                                    k_vect,
                                                                    optical_direction,
                                                                    m_optical_derivative_step,
                                                                    m_nonlocal_epm);
            }
            std::vector<double> list_k_sum_real(m_energies.size());
            std::vector<double> list_k_sum_imag(m_energies.size());
            for (int idx_conduction_band = index_first_conduction_band; idx_conduction_band < m_nb_bands;
                 ++idx_conduction_band) {
                for (int idx_valence_band = 0; idx_valence_band < index_first_conduction_band; ++idx_valence_band) {
                    double transition_strength = 0.0;
                    double delta_energy        = 0.0;
                    if (optical_limit) {
                        delta_energy =
                            m_eigenvalues_k[index_k][idx_conduction_band] - m_eigenvalues_k[index_k][idx_valence_band];
                        transition_strength = optical_transition_strength(
                            m_material,
                            m_basisVectors,
                            k_vect,
                            optical_direction,
                            m_eigenvectors_k[index_k],
                            use_finite_difference_derivative ? &finite_difference_derivative : nullptr,
                            idx_conduction_band,
                            idx_valence_band,
                            delta_energy);
                    } else {
                        transition_strength = std::norm(eigenvectors_k_plus_q.col(idx_conduction_band)
                                                            .dot(m_eigenvectors_k[index_k].col(idx_valence_band)));
                        delta_energy =
                            eigenvalues_k_plus_q[idx_conduction_band] - m_eigenvalues_k[index_k][idx_valence_band];
                    }
                    for (std::size_t index_energy = 0; index_energy < m_energies.size(); ++index_energy) {
                        double energy        = m_energies[index_energy];
                        double delta_minus_e = delta_energy - energy;
                        double delta_plus_e  = delta_energy + energy;
                        double denom_minus   = delta_minus_e * delta_minus_e + eta_smearing * eta_smearing;
                        double denom_plus    = delta_plus_e * delta_plus_e + eta_smearing * eta_smearing;
                        double real_factor   = delta_minus_e / denom_minus + delta_plus_e / denom_plus;
                        double imag_factor   = eta_smearing / denom_minus - eta_smearing / denom_plus;
                        list_k_sum_real[index_energy] += transition_strength * real_factor;
                        list_k_sum_imag[index_energy] += transition_strength * imag_factor;
                    }
                }
            }
            if (m_qpoints.size() <= 1) {
                // if there is only one q point in the list, we don't keep the eigenvectors, to save memory.
                m_eigenvectors_k[index_k].resize(1, 1);
                m_eigenvalues_k[index_k].resize(1);
            }
            for (std::size_t index_energy = 0; index_energy < m_energies.size(); ++index_energy) {
                list_total_sum_real[index_energy] += list_k_sum_real[index_energy];
                list_total_sum_imag[index_energy] += list_k_sum_imag[index_energy];
            }
        }
        std::vector<double> list_epsilon_real(m_energies.size());
        std::vector<double> list_epsilon_imag(m_energies.size());
        for (std::size_t index_energy = 0; index_energy < m_energies.size(); ++index_energy) {
            list_epsilon_real[index_energy] = list_total_sum_real[index_energy];
            list_epsilon_imag[index_energy] = list_total_sum_imag[index_energy];
        }
        auto end     = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
        if (mpi_rank == 0) {
            if (optical_limit) {
                std::cout << "Computed optical-limit dielectric function for direction = " << optical_direction
                          << " -> " << index_q + 1 << "/" << m_qpoints.size() << " in " << elapsed << " s" << std::endl;
            } else {
                std::cout << "Computed dielectric function for q = " << m_qpoints[index_q] << " -> " << index_q + 1
                          << "/" << m_qpoints.size() << " in " << elapsed << " s" << std::endl;
            }
        }
        start = std::chrono::high_resolution_clock::now();
        m_dielectric_function_real.push_back(list_epsilon_real);
        m_dielectric_function_imag.push_back(list_epsilon_imag);
    }
}

DielectricFunction DielectricFunction::merge_results(
    DielectricFunction                                   RootDielectricFunction,
    const std::vector<std::vector<std::vector<double>>>& dielectric_function_real_results,
    const std::vector<std::vector<std::vector<double>>>& dielectric_function_imag_results,
    std::vector<int>                                     nb_kpoints_per_instance) {
    std::vector<std::vector<double>> total_dielectric_function_real;
    std::vector<std::vector<double>> total_dielectric_function_imag;
    if (dielectric_function_real_results.empty() || dielectric_function_imag_results.empty()) {
        throw std::runtime_error("No results to merge");
    }
    if (dielectric_function_real_results.size() != nb_kpoints_per_instance.size() ||
        dielectric_function_imag_results.size() != nb_kpoints_per_instance.size()) {
        throw std::runtime_error("Number of results and number of k-points per instance do not match");
    }
    std::size_t total_number_kpoints = 0;
    for (std::size_t index_instance = 0; index_instance < nb_kpoints_per_instance.size(); ++index_instance) {
        if (nb_kpoints_per_instance[index_instance] < 0) {
            throw std::runtime_error("Negative k-point count while merging dielectric-function results.");
        }
        total_number_kpoints += nb_kpoints_per_instance[index_instance];
    }
    if (total_number_kpoints == 0) {
        throw std::runtime_error("Cannot merge dielectric-function results without k-points.");
    }
    std::cout << "Number total kpoint to the merge : " << total_number_kpoints << std::endl;
    const auto merge_component = [](const std::vector<std::vector<std::vector<double>>>& component_results) {
        std::vector<std::vector<double>> total_component;
        for (std::size_t index_instance = 0; index_instance < component_results.size(); ++index_instance) {
            for (std::size_t index_q = 0; index_q < component_results[index_instance].size(); ++index_q) {
                if (index_instance == 0) {
                    total_component.push_back(component_results[index_instance][index_q]);
                } else {
                    for (std::size_t index_energy = 0; index_energy < component_results[index_instance][index_q].size();
                         ++index_energy) {
                        total_component[index_q][index_energy] +=
                            component_results[index_instance][index_q][index_energy];
                    }
                }
            }
        }
        return total_component;
    };

    total_dielectric_function_real = merge_component(dielectric_function_real_results);
    total_dielectric_function_imag = merge_component(dielectric_function_imag_results);

    // Renormalization
    double renormalization = 1.0 / static_cast<double>(total_number_kpoints);
    std::cout << "Renormalization: " << renormalization << std::endl;
    for (std::size_t index_q = 0; index_q < total_dielectric_function_real.size(); ++index_q) {
        Vector3D<double> q_vect    = RootDielectricFunction.m_qpoints[index_q];
        double           q_squared = pow(q_vect.Length(), 2);
        if (RootDielectricFunction.m_response_mode == DielectricResponseMode::optical_limit) {
            q_squared = 1.0;
        }
        if (!(q_squared > 0.0) || !std::isfinite(q_squared)) {
            throw std::runtime_error("Cannot merge dielectric-function results at singular q = 0.");
        }
        const double dielectric_prefactor = 2.0 * M_PI / q_squared;
        for (std::size_t index_energy = 0; index_energy < total_dielectric_function_real[index_q].size();
             ++index_energy) {
            total_dielectric_function_real[index_q][index_energy] *= renormalization;
            total_dielectric_function_imag[index_q][index_energy] *= renormalization;
            total_dielectric_function_real[index_q][index_energy] =
                1.0 + dielectric_prefactor * total_dielectric_function_real[index_q][index_energy];
            total_dielectric_function_imag[index_q][index_energy] =
                dielectric_prefactor * total_dielectric_function_imag[index_q][index_energy];
        }
    }

    DielectricFunction dielectric_function         = RootDielectricFunction;
    dielectric_function.m_dielectric_function_real = total_dielectric_function_real;
    dielectric_function.m_dielectric_function_imag = total_dielectric_function_imag;
    return dielectric_function;
}

Eigen::MatrixXd create_kramers_matrix(const std::vector<double>& energies) {
    const std::size_t N = energies.size();
    if (N < 2) {
        throw std::invalid_argument("Kramers-Kronig transform requires at least two energy samples.");
    }
    const double d_energy = energies[1] - energies[0];
    if (!(d_energy > 0.0)) {
        throw std::invalid_argument("Kramers-Kronig energy grid must be strictly increasing.");
    }
    for (std::size_t idx = 2; idx < N; ++idx) {
        const double step = energies[idx] - energies[idx - 1];
        if (std::abs(step - d_energy) > 1e-9 * std::max(1.0, std::abs(d_energy))) {
            throw std::invalid_argument("Kramers-Kronig transform currently requires a uniform energy grid.");
        }
    }

    Eigen::MatrixXd kramers_matrix(N, N);
    for (std::size_t idx_line = 0; idx_line < N; ++idx_line) {
        for (std::size_t idx_col = 0; idx_col < N; ++idx_col) {
            if (idx_line == idx_col) {
                kramers_matrix(idx_line, idx_col) = 0.0;
            } else {
                const double omega_j              = energies[idx_line];
                const double omega_k              = energies[idx_col];
                const double denominator          = omega_j * omega_j - omega_k * omega_k;
                kramers_matrix(idx_line, idx_col) = (2.0 * omega_j / M_PI) * d_energy / denominator;
            }
        }
    }
    return kramers_matrix;
}

void DielectricFunction::apply_kramers_kronig() {
    std::cout << "Applying Kramers-Kronig" << std::endl;
    m_dielectric_function_imag.clear();
    m_dielectric_function_imag.resize(m_dielectric_function_real.size());
    Eigen::MatrixXd kramers_matrix = create_kramers_matrix(m_energies);
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

void DielectricFunction::export_dielectric_function_at_q(const std::string& filename,
                                                         std::size_t        idx_q,
                                                         bool               name_auto) const {
    std::string outname;
    if (name_auto) {
        // outname = m_export_prefix + '_' + std::to_string(idx_q) + '_' + std::to_string(m_qpoints[idx_q].X) + "_" +
        // std::to_string(m_qpoints[idx_q].Y) + "_" +
        //           std::to_string(m_qpoints[idx_q].Z) + ".csv";
        outname = fmt::format("{}_{:05}_{:.9g}_{:.9g}_{:.9g}.csv",
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
