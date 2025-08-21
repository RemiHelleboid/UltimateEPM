/**
 * @file elelectron_phonon.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-02-09
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "electron_phonon.hpp"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <tuple>

#include "Constants.hpp"
#include "Vector3D.h"
#include "bz_states.hpp"
#include "gmsh.h"
#include "omp.h"

namespace bz_mesh {

// double ElectronPhonon::bose_einstein_distribution(double energy, double temperature) {
//     if (energy < 1e-9) {
//         return 0.0;
//     }
//     double f = 1.0 / (std::expm1(energy / (EmpiricalPseudopotential::Constants::k_b_eV * temperature)));
//     if (f < 0 || std::isnan(f) || std::isinf(f)) {
//         std::cout << "Energy: " << energy << " Temperature: " << temperature << " f: " << f << std::endl;
//         throw std::runtime_error("Bose-Einstein distribution is negative or NaN");
//     }
//     return f;
// }

double ElectronPhonon::bose_einstein_distribution(double energy, double temperature) {
    double x = energy / (EmpiricalPseudopotential::Constants::k_b_eV * temperature);
    return 1.0 / std::expm1(x);
}

double ElectronPhonon::electron_overlap_integral(const Vector3D<double>& k1, const Vector3D<double>& k2) {
    constexpr double R_Wigner_Seitz = 2.122e-10;
    const double     qRws           = (k1 - k2).Length() * R_Wigner_Seitz;
    double           integral       = 3.0 * (std::sin(qRws) - qRws * std::cos(qRws)) / (qRws * qRws * qRws);
    return integral;
}

double ElectronPhonon::hole_overlap_integral(int n1, const Vector3D<double>& k1, int n2, const Vector3D<double>& k2) {
    const double cos_angle_k1_k2_2 = compte_cos_angle(k1, k2) * compte_cos_angle(k1, k2);
    auto         A_B_params        = m_hole_overlap_int_params.get_params(n1, n2);
    double       integral          = (1.0 / 2.0) * std::sqrt(A_B_params[0] + A_B_params[1] * cos_angle_k1_k2_2);
    return integral;
}

/**
 * @brief Compute the electron-phonon rate for a given electron band and k-point.
 *
 * \begin{aligned}
P_{\eta, \nu}(k)=\sum_{\nu^{\prime}, m} \frac{\pi}{\rho \omega_\eta\left(q_m\right)} \Delta_\eta^2 & \times\left|\mathcal{I}\left(k,
k_m^{\prime}\right)\right|^2 \times\left[n_\eta\left(q_m\right)+1 / 2 \pm 1 / 2\right]  \times \mathcal{D}_m\left(E_{\nu^{\prime}}^{\prime}\right)
\end{aligned}
 *
 * @param idx_n1
 * @param idx_k1
 * @return RateValues
 */
// RateValues ElectronPhonon::compute_electron_phonon_rate(int idx_n1, std::size_t idx_k1) {
//     RateValues                rates_k1_n1;
//     const std::vector<Tetra>& list_tetrahedra          = m_list_tetrahedra;
//     auto                      indices_conduction_bands = m_indices_conduction_bands;
//     double                    energy_n1_k1             = m_list_vertices[idx_k1].get_energy_at_band(idx_n1);
    
//     auto             k_1 = m_list_vertices[idx_k1].get_position();
//     for (auto&& idx_n2 : indices_conduction_bands) {
//         for (auto&& tetra : list_tetrahedra) {
//             auto             k_2 = tetra.compute_barycenter();
//             Vector3D<double> k1{k_1.x(), k_1.y(), k_1.z()};
//             Vector3D<double> k2{k_2.x(), k_2.y(), k_2.z()};

//             double           overlap_integral = electron_overlap_integral(k1, k2);
//             auto             q                = k2 - k1;
//             auto             initial_q        = q;
//             Vector3D<double> q_ph{q.X, q.Y, q.Z};
//             bool             is_in_bz = is_inside_mesh_geometry(q_ph);
//             // std::cout << "Qxyz : " << q_ph  << " is in BZ: " << is_in_bz << std::endl;
//             // No Umklapp process for now.
//             if (!is_in_bz) {
//                 // continue;
//                 auto new_q = retrieve_k_inside_mesh_geometry(q_ph);
//                 q_ph       = Vector3D<double>{new_q.x(), new_q.y(), new_q.z()};
//                 // std::cout << "q: " << q_ph / m_material.get_fourier_factor() << " new_q: " << new_q / m_material.get_fourier_factor() <<
//                 // std::endl;
//             }
//             is_in_bz = is_inside_mesh_geometry(q_ph);
//             if (!is_in_bz) {
//                 throw std::runtime_error("Q is not inside the BZ");
//             }

//             for (auto&& mode_phonon : m_phonon_dispersion) {
//                 PhononModeDirection mode_direction = mode_phonon.first;
//                 double e_ph = mode_phonon.second.get_phonon_dispersion(q_ph.Length()) * EmpiricalPseudopotential::Constants::h_bar_eV;
//                 if (e_ph > 100e-3) {
//                     std::cout << "Energy phonon: " << e_ph << std::endl << std::endl;
//                     throw std::runtime_error("Energy phonon too high");
//                 }
//                 if (e_ph <= 0.0) {
//                     std::cout << "Energy phonon: " << e_ph << std::endl << std::endl;
//                     std::cout << "Q: " << q_ph << std::endl << std::endl;
//                     std::cout << "Q: " << q_ph.Length() / m_material.get_fourier_factor() << std::endl << std::endl;
//                     std::cout << "Energy (n1, k1): " << energy_n1_k1 << std::endl << std::endl;
//                     bool is_in_bz = is_inside_mesh_geometry(q_ph);
//                     std::cout << "Is in BZ: " << is_in_bz << std::endl << std::endl;
//                     throw std::runtime_error("Energy phonon negative");
//                     continue;
//                 }
//                 // std::cout << "Energy phonon: " << e_ph << std::endl;
//                 for (auto sign_phonon : {-1.0, 1.0}) {
//                     const double N0           = bose_einstein_distribution(e_ph, m_temperature);
//                     const double bose_part    = (sign_phonon == 1.0) ? N0           // absorption
//                                                                      : (N0 + 1.0);  // emission
//                     double energy_final = energy_n1_k1 + e_ph * sign_phonon;
//                     // std::cout << "Energy (n1, k1): " << energy_n1_k1 << " Energy phonon: " << e_ph << " Energy final: " << energy_final
//                     //           << std::endl;
//                     if (!tetra.is_energy_inside_band(energy_final, idx_n2)) {
//                         continue;
//                     }
//                     // std::cout << "Volume: " << this->get_volume() << std::endl;
//                     double dos_tetra = tetra.compute_tetra_dos_energy_band(energy_final, idx_n2);
//                     double dos_joules = dos_tetra / EmpiricalPseudopotential::Constants::q_e;
//                     // std::cout << "Energy final: " << energy_final << " dos tetra: " << dos_tetra
//                     //           << std::endl;
//                     DeformationPotential deformation_potential       = (mode_direction.first == PhononMode::acoustic)
//                                                                            ? m_acoustic_deformation_potential_e
//                                                                            : m_optical_deformation_potential_e;
//                     double               deformation_potential_value = deformation_potential.get_deformation_potential(q_ph, energy_final) * EmpiricalPseudopotential::Constants::q_e;
//                     PhononMode           phonon_mode                 = mode_direction.first;
//                     PhononDirection      phonon_direction            = mode_direction.second;
//                     PhononEvent          phonon_event = (sign_phonon == 1.0) ? PhononEvent::absorption : PhononEvent::emission;

//                     double raw_rate_value = (dos_joules / (e_ph * EmpiricalPseudopotential::Constants::q_e)) *
//                                             (bose_part)*overlap_integral * overlap_integral * deformation_potential_value *
//                                             deformation_potential_value;
//                     // std::cout << "Raw rate value: " << raw_rate_value << std::endl;
//                     double rate_value = (EmpiricalPseudopotential::Constants::pi / (m_rho * e_ph)) * deformation_potential_value *
//                                         deformation_potential_value * overlap_integral * overlap_integral * bose_part * dos_joules;
//                     // std::cout << this->get_volume() << std::endl;
//                     // rate_value /= this->get_volume();
//                     // rate_value *= EmpiricalPseudopotential::Constants::q_e;

//                     if (0 || rate_value < 0.0 || rate_value > 1e50 || std::isnan(rate_value) || std::isinf(rate_value)) {
//                         std::cout << "k1 n1: " << idx_k1 << " " << idx_n1 << std::endl;
//                         std::cout << "Rate value: " << rate_value << std::endl;
//                         std::cout << "Overlap integral: " << overlap_integral << std::endl;
//                         std::cout << "Deformation potential: " << deformation_potential_value << std::endl;
//                         std::cout << "Bose part: " << bose_part << std::endl;
//                         std::cout << "DOS tetra: " << dos_tetra << std::endl;
//                         std::cout << "Energy final: " << energy_final << std::endl;
//                         std::cout << "q: " << q_ph << std::endl;
//                         std::cout << "initial q: " << initial_q / m_material.get_fourier_factor() << std::endl;
//                         std::cout << "Energy phonon: " << e_ph << std::endl;
//                         std::cout << "Energy (n1, k1): " << energy_n1_k1 << std::endl;
//                         std::cout << "--------------------------------------------------------------------------------" << std::endl;
//                     }
//                     // std::cout << "Rate value: " << rate_value << std::endl;
//                     RateValue rate(phonon_mode, phonon_direction, phonon_event, rate_value);
//                     rates_k1_n1.add_rate(rate);
//                 }
//             }
//         }
//     }
//     rates_k1_n1.print_rates();
//     std::cout << "*************************************************************************************************************" <<
//     std::endl;
//     return rates_k1_n1;
// }

RateValues ElectronPhonon::compute_electron_phonon_rate(int idx_n1, std::size_t idx_k1) {
    RateValues                rates_k1_n1;
    const std::vector<Tetra>& list_tetrahedra          = m_list_tetrahedra;
    const auto&               indices_conduction_bands = m_indices_conduction_bands;

    // Initial state energy (eV) and k (in SI units 1/m)
    const double           energy_n1_k1 = m_list_vertices[idx_k1].get_energy_at_band(idx_n1);
    const auto             k1_v3        = m_list_vertices[idx_k1].get_position();
    const Vector3D<double> k1{k1_v3.x(), k1_v3.y(), k1_v3.z()};
    std::fstream file_error_k_points("error_k_points.txt", std::ios::out | std::ios::app);

    for (int idx_n2 : indices_conduction_bands) {
        for (const auto& tetra : list_tetrahedra) {
            // Final k chosen as tetra barycenter (SI)
            const auto             k2_v3 = tetra.compute_barycenter();
            const Vector3D<double> k2{k2_v3.x(), k2_v3.y(), k2_v3.z()};

            // Overlap (dimensionless with your definition)
            const double overlap2 = std::pow(electron_overlap_integral(k1, k2), 2);

            // Phonon wavevector q = k2 - k1 (SI)
            const Vector3D<double> q = k2 - k1;
            Vector3D<double>       q_ph{q.X, q.Y, q.Z};

            // Fold q back to first BZ if needed (no Umklapp yet)
            if (!is_inside_mesh_geometry(q_ph)) {
                try {
                    // Fold q to the first BZ
                    const auto q_folded = retrieve_k_inside_mesh_geometry({q_ph.X, q_ph.Y, q_ph.Z});
                    q_ph                = Vector3D<double>(q_folded.x(), q_folded.y(), q_folded.z());
                } catch (const std::runtime_error& e) {
                    std::cerr << "Error folding q: " << e.what() << "\n";
                    file_error_k_points << q_ph.X << " " << q_ph.Y << " " << q_ph.Z << "\n";
                    continue;  // Skip this tetrahedron if folding fails
                }
            }
            if (!is_inside_mesh_geometry(q_ph)) throw std::runtime_error("Q is not inside the BZ");

            // Loop phonon branches
            for (const auto& kv : m_phonon_dispersion) {
                const PhononModeDirection mode_direction = kv.first;
                const auto&               disp           = kv.second;

                // ---- Phonon quantities
                // omega [1/s] from dispersion(|q|)
                const double omega = disp.get_phonon_dispersion(q_ph.Length());
                if (omega <= 0.0) continue;  // skip unphysical/zero frequency points

                // E_ph in eV for Bose factor
                const double Eph_eV = EmpiricalPseudopotential::Constants::h_bar_eV * omega;

                // Basic sanity
                if (Eph_eV > 0.1)  // 100 meV cap as you had
                    throw std::runtime_error("Energy phonon too high");

                // ---- Loop emission/absorption
                for (double sign_phonon : {-1.0, 1.0}) {
                    // Bose factor in eV units
                    const double N0        = bose_einstein_distribution(Eph_eV, m_temperature);  // dimensionless
                    const double bose_part = (sign_phonon > 0.0) ? (N0 + 1.0) : N0;              // +1: emission, plain N0: absorption

                    // Final electronic energy (eV)
                    const double E_final_eV = energy_n1_k1 + sign_phonon * Eph_eV;
                    if (!tetra.is_energy_inside_band(E_final_eV, idx_n2)) continue;

                    // Electronic DOS from tetra method:
                    //  - your function returns per eV (assumed) — convert to per J
                    //  - if your DOS is *total* (not per m^3), keep it; if it is per volume, that’s fine too:
                    //    the formula here expects the same DOS convention you used in the derivation.
                    const double dos_eV = tetra.compute_tetra_dos_energy_band(E_final_eV, idx_n2);
                    const double dos_per_J =
                        dos_eV / EmpiricalPseudopotential::Constants::q_e;  // [1/J] (or [1/(J·m^3)] if it already had /m^3)

                    // Deformation potential in SI Joules
                    const DeformationPotential defpot = (mode_direction.first == PhononMode::acoustic) ? m_acoustic_deformation_potential_e
                                                                                                       : m_optical_deformation_potential_e;
                    const double               Delta_J =
                        defpot.get_deformation_potential(q_ph, E_final_eV) * EmpiricalPseudopotential::Constants::q_e;  // eV -> J

                    // ---- Final rate (your formula):
                    // P = (pi / (rho * omega)) * Delta^2 * |I|^2 * (bose) * DOS(E)
                    // Units will be s^-1 if DOS matches the convention used in the derivation (see our notes).
                    const double rate_value = (EmpiricalPseudopotential::Constants::pi / (m_rho * omega)) * (Delta_J * Delta_J) * overlap2 *
                                              bose_part * dos_per_J;

                    // Store
                    const PhononMode      phonon_mode      = mode_direction.first;
                    const PhononDirection phonon_direction = mode_direction.second;
                    const PhononEvent     phonon_event     = (sign_phonon > 0.0) ? PhononEvent::emission : PhononEvent::absorption;

                    // (Optional) sanity guard
                    if (!(rate_value >= 0.0) || !std::isfinite(rate_value)) {
                        std::cerr << "Bad rate: " << rate_value << "  (omega=" << omega << " s^-1, Eph=" << Eph_eV
                                  << " eV, DOS_J=" << dos_per_J << ", Delta_J=" << Delta_J << ", |I|^2=" << overlap2
                                  << ", bose=" << bose_part << ")\n";
                        continue;
                    }

                    rates_k1_n1.add_rate(RateValue(phonon_mode, phonon_direction, phonon_event, rate_value));
                }
            }
        }
    }

    // rates_k1_n1.print_rates();
    // std::cout << "*************************************************************************************************************\n";
    return rates_k1_n1;
}

RateValues ElectronPhonon::compute_hole_phonon_rate(int idx_n1, std::size_t idx_k1) {
    RateValues                rates_k1_n1;
    const std::vector<Tetra>& list_tetrahedra       = m_list_tetrahedra;
    auto                      indices_valence_bands = m_indices_valence_bands;
    double                    energy_n1_k1          = m_list_vertices[idx_k1].get_energy_at_band(idx_n1);

    for (auto&& idx_n2 : indices_valence_bands) {
        for (auto&& tetra : list_tetrahedra) {
            auto             k_1 = m_list_vertices[idx_k1].get_position();
            auto             k_2 = tetra.compute_barycenter();
            Vector3D<double> k1{k_1.x(), k_1.y(), k_1.z()};
            Vector3D<double> k2{k_2.x(), k_2.y(), k_2.z()};

            double           overlap_integral = hole_overlap_integral(idx_n1, k1, idx_n2, k2);
            auto             q                = k2 - k1;
            auto             initial_q        = q;
            Vector3D<double> q_ph{q.X, q.Y, q.Z};
            bool             is_in_bz = is_inside_mesh_geometry(q_ph / m_material.get_fourier_factor());

            // No Umklapp process for now.
            if (!is_in_bz) {
                auto new_q = retrieve_k_inside_mesh_geometry(q_ph / m_material.get_fourier_factor()) * m_material.get_fourier_factor();
                q_ph       = Vector3D<double>{new_q.x(), new_q.y(), new_q.z()};
                // std::cout << "q: " << q_ph / m_material.get_fourier_factor() << " new_q: " << new_q / m_material.get_fourier_factor() <<
                // std::endl;
            }
            is_in_bz = is_inside_mesh_geometry(q_ph / m_material.get_fourier_factor());
            if (!is_in_bz) {
                throw std::runtime_error("Q is not inside the BZ");
            }

            for (auto&& mode_phonon : m_phonon_dispersion) {
                PhononModeDirection mode_direction = mode_phonon.first;
                double e_ph = mode_phonon.second.get_phonon_dispersion(q_ph.Length()) * EmpiricalPseudopotential::Constants::h_bar_eV;
                if (e_ph > 100e-3) {
                    std::cout << "Energy phonon: " << e_ph << std::endl << std::endl;
                    throw std::runtime_error("Energy phonon too high");
                }
                if (e_ph <= 0.0) {
                    // std::cout << "Energy phonon: " << e_ph << std::endl << std::endl;
                    // std::cout << "Q: " << q_ph << std::endl << std::endl;
                    // std::cout << "Q: " << q_ph.Length() / m_material.get_fourier_factor() << std::endl << std::endl;
                    // std::cout << "Energy (n1, k1): " << energy_n1_k1 << std::endl << std::endl;
                    // bool is_in_bz = is_inside_mesh_geometry(q_ph / m_material.get_fourier_factor());
                    // std::cout << "Is in BZ: " << is_in_bz << std::endl << std::endl;
                    // throw std::runtime_error("Energy phonon negative");
                    continue;
                }
                // std::cout << "Energy phonon: " << e_ph << std::endl;
                for (auto sign_phonon : {-1.0, 1.0}) {
                    double bose_part    = (bose_einstein_distribution(e_ph, m_temperature) + 0.5 + 0.5 * std::pow(-1.0, sign_phonon));
                    double energy_final = energy_n1_k1 + e_ph * sign_phonon;
                    // std::cout << "Energy (n1, k1): " << energy_n1_k1 << " Energy phonon: " << e_ph << " Energy final: " << energy_final
                    //           << std::endl;
                    if (!tetra.is_energy_inside_band(energy_final, idx_n2)) {
                        continue;
                    }
                    // std::cout << "Volume: " << this->get_volume() << std::endl;
                    double dos_tetra = tetra.compute_tetra_dos_energy_band(energy_final, idx_n2);
                    // double dos_tetra = 1.0;
                    // std::cout << "DOS: " << dos_tetra << std::endl;
                    DeformationPotential deformation_potential       = (mode_direction.first == PhononMode::acoustic)
                                                                           ? m_acoustic_deformation_potential_h
                                                                           : m_optical_deformation_potential_h;
                    double               deformation_potential_value = deformation_potential.get_deformation_potential(q_ph, energy_final);
                    PhononMode           phonon_mode                 = mode_direction.first;
                    PhononDirection      phonon_direction            = mode_direction.second;
                    PhononEvent          phonon_event = (sign_phonon == 1.0) ? PhononEvent::absorption : PhononEvent::emission;
                    double rate_value = (EmpiricalPseudopotential::Constants::pi / (m_rho * e_ph)) * deformation_potential_value *
                                        deformation_potential_value * overlap_integral * overlap_integral * bose_part * dos_tetra;
                    // std::cout << this->get_volume() << std::endl;
                    // rate_value /= this->get_volume();
                    // rate_value *= EmpiricalPseudopotential::Constants::q_e;

                    if (rate_value < 0.0 || rate_value > 1e50 || std::isnan(rate_value) || std::isinf(rate_value)) {
                        std::cout << "Rate value: " << rate_value << std::endl;
                        std::cout << "Overlap integral: " << overlap_integral << std::endl;
                        std::cout << "Deformation potential: " << deformation_potential_value << std::endl;
                        std::cout << "Bose part: " << bose_part << std::endl;
                        std::cout << "DOS tetra: " << dos_tetra << std::endl;
                        std::cout << "Energy final: " << energy_final << std::endl;
                        std::cout << "q: " << q_ph << std::endl;
                        std::cout << "initial q: " << initial_q / m_material.get_fourier_factor() << std::endl;
                        std::cout << "Energy phonon: " << e_ph << std::endl;
                        std::cout << "Energy (n1, k1): " << energy_n1_k1 << std::endl;
                        std::cout << "--------------------------------------------------------------------------------" << std::endl;
                    }
                    // std::cout << "Rate value: " << (deformation_potential_value) << std::endl;
                    RateValue rate(phonon_mode, phonon_direction, phonon_event, rate_value);
                    rates_k1_n1.add_rate(rate);
                }
            }
        }
    }
    // rates_k1_n1.print_rates();
    // std::cout << "*************************************************************************************************************" <<
    // std::endl;
    return rates_k1_n1;
}

void ElectronPhonon::compute_electron_phonon_rates_over_mesh() {
    auto indices_conduction_bands = m_indices_conduction_bands;
    auto min_idx_conduction_band  = *std::min_element(indices_conduction_bands.begin(), indices_conduction_bands.end());
    auto max_idx_conduction_band  = *std::max_element(indices_conduction_bands.begin(), indices_conduction_bands.end());
    std::cout << "Min index conduction band: " << min_idx_conduction_band << std::endl;
    std::ofstream file("rates.txt");

    std::random_device                     rd;
    std::mt19937                           gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    double                                 p_compute_rate = 1.0;

    std::cout << "Computing electron-phonon rates over mesh for " << m_list_vertices.size() * p_compute_rate << " k-points." << std::endl;

#pragma omp parallel for schedule(dynamic)
    for (std::size_t idx_k1 = 0; idx_k1 < m_list_vertices.size(); ++idx_k1) {
        double r = dis(gen);
        for (std::size_t idx_n1 = 0; idx_n1 < min_idx_conduction_band; ++idx_n1) {
            // DEBUG PROVISOIRE
            if (1e10 > p_compute_rate) {
                std::array<double, 8> array = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                m_list_vertices[idx_k1].add_electron_phonon_rates(array);
                continue;
            }
            auto hole_rate = compute_hole_phonon_rate(idx_n1, idx_k1);
            auto array     = hole_rate.to_array();
            m_list_vertices[idx_k1].add_electron_phonon_rates(array);
        }

        for (std::size_t idx_n1 = min_idx_conduction_band; idx_n1 <= max_idx_conduction_band; ++idx_n1) {
            if (r > p_compute_rate) {
                std::array<double, 8> array = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                m_list_vertices[idx_k1].add_electron_phonon_rates(array);
                continue;
            } else {
                auto rate  = compute_electron_phonon_rate(idx_n1, idx_k1);
                auto array = rate.to_array();
                m_list_vertices[idx_k1].add_electron_phonon_rates(array);
            }
        }
        if (omp_get_thread_num() == 0) {
            std::cout << "\rComputing rates for k-point " << idx_k1 << " / " << m_list_vertices.size() << std::flush;
        }
    }
    std::cout << std::endl;
    file.close();
    // Add rate to the mesh
}

void ElectronPhonon::compute_plot_electron_phonon_rates_vs_energy_over_mesh(int                nb_bands,
                                                                            double             max_energy,
                                                                            double             energy_step,
                                                                            const std::string& filename) {
    std::ofstream file(filename);
    if (!file) throw std::runtime_error("Cannot open " + filename);

    // Choose band range: here I use the actual conduction band indices.
    // If you want valence, swap to m_indices_valence_bands.
    if (m_indices_conduction_bands.empty()) throw std::runtime_error("No conduction bands indexed.");
    const int band_min = *std::min_element(m_indices_conduction_bands.begin(), m_indices_conduction_bands.end());
    const int band_max = *std::max_element(m_indices_conduction_bands.begin(), m_indices_conduction_bands.end());
    // Optionally clamp to nb_bands if caller insists:
    const int last_band = std::min(band_max, nb_bands - 1);

    // Optional: CSV-like header
    file << "# E[eV]  DOS(E) "
         << " hole_LO_em  hole_LO_ab  hole_TR_em  hole_TR_ab"
         << " elec_LO_em  elec_LO_ab  elec_TR_em  elec_TR_ab\n";

    for (double E = 0.0; E < max_energy; E += energy_step) {
        std::cout << "\rEnergy: " << E << " / " << max_energy << std::flush;

        double                dos_E = 0.0;
        std::array<double, 8> num_accum{};  // numerator accumulators
        num_accum.fill(0.0);

        // Sum over tetrahedra and bands
        for (const auto& tetra : m_list_tetrahedra) {
            for (int b = band_min; b <= last_band; ++b) {
                const double dos_tetra = tetra.compute_tetra_dos_energy_band(E, b);
                if (!std::isfinite(dos_tetra)) {
                    std::ostringstream oss;
                    oss << "DOS tetra is NaN/Inf at E=" << E << " band=" << b;
                    throw std::runtime_error(oss.str());
                }
                if (dos_tetra <= 0.0) continue;

                dos_E += dos_tetra;

                const std::array<double, 8> r = tetra.get_mean_electron_phonon_rates(b);
                for (std::size_t i = 0; i < r.size(); ++i)
                    num_accum[i] += r[i] * dos_tetra;  // DOS-weighted sum
            }
        }

        // Per-energy normalization (NOT by total across energies)
        std::array<double, 8> mean_rates{};
        if (dos_E > 0.0) {
            for (std::size_t i = 0; i < mean_rates.size(); ++i)
                mean_rates[i] = num_accum[i] / dos_E;
        } else {
            mean_rates.fill(0.0);  // or continue; to skip the line
        }

        // write line
        file << std::setprecision(10) << E << " " << dos_E << " ";
        for (double v : mean_rates)
            file << v << " ";
        file << "\n";
    }

    std::cout << std::endl;
    file.close();
}

void ElectronPhonon::plot_phonon_dispersion(const std::string& filename) const {
    std::ofstream file(filename);
    for (auto&& vtx : m_list_vertices) {
        auto k = vtx.get_position();
        file << k.x() << " " << k.y() << " " << k.z() << " ";
        for (auto&& mode_direction : m_phonon_dispersion) {
            auto q = k;
            q /= m_material.get_fourier_factor();
            auto e_ph = mode_direction.second.get_phonon_dispersion(q.norm()) * EmpiricalPseudopotential::Constants::h_bar_eV;
            file << e_ph << " ";
        }
        file << std::endl;
    }
}

void ElectronPhonon::add_electron_phonon_rates_to_mesh(const std::string& initial_filename, const std::string& final_filename) {
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 0);
    gmsh::model::add("bz_mesh");
    gmsh::open(initial_filename);

    std::string model_file_name;
    gmsh::model::getCurrent(model_file_name);

    std::vector<std::size_t> node_tags;
    std::vector<double>      nodeCoords;
    std::vector<double>      nodeParams;
    gmsh::model::mesh::reclassifyNodes();
    gmsh::model::mesh::getNodes(node_tags, nodeCoords, nodeParams, -1, -1, false, false);

    auto indices_conduction_bands = m_indices_conduction_bands;
    auto min_idx_conduction_band  = *std::min_element(indices_conduction_bands.begin(), indices_conduction_bands.end());
    auto max_idx_conduction_band  = *std::max_element(indices_conduction_bands.begin(), indices_conduction_bands.end());
    std::cout << "Min index conduction band: " << min_idx_conduction_band << std::endl;
    int nb_bands = m_indices_conduction_bands.size() + m_indices_valence_bands.size();

    for (int idx_val_band = 0; idx_val_band < max_idx_conduction_band; ++idx_val_band) {
        std::vector<double> rates_hole_lo_em(m_list_vertices.size());
        std::vector<double> rates_hole_lo_ab(m_list_vertices.size());
        std::vector<double> rates_hole_tr_em(m_list_vertices.size());
        std::vector<double> rates_hole_tr_ab(m_list_vertices.size());
        std::vector<double> rates_elec_lo_em(m_list_vertices.size());
        std::vector<double> rates_elec_lo_ab(m_list_vertices.size());
        std::vector<double> rates_elec_tr_em(m_list_vertices.size());
        std::vector<double> rates_elec_tr_ab(m_list_vertices.size());

        for (std::size_t idx_k1 = 0; idx_k1 < m_list_vertices.size(); ++idx_k1) {
            auto rates               = m_list_vertices[idx_k1].get_electron_phonon_rates(idx_val_band);
            rates_hole_lo_em[idx_k1] = rates[0];
            rates_hole_lo_ab[idx_k1] = rates[1];
            rates_hole_tr_em[idx_k1] = rates[2];
            rates_hole_tr_ab[idx_k1] = rates[3];
            rates_elec_lo_em[idx_k1] = rates[4];
            rates_elec_lo_ab[idx_k1] = rates[5];
            rates_elec_tr_em[idx_k1] = rates[6];
            rates_elec_tr_ab[idx_k1] = rates[7];
        }
        std::string name_rate_hole_lo_em = "hole_lo_em_" + std::to_string(idx_val_band);
        std::string name_rate_hole_lo_ab = "hole_lo_ab_" + std::to_string(idx_val_band);
        std::string name_rate_hole_tr_em = "hole_tr_em_" + std::to_string(idx_val_band);
        std::string name_rate_hole_tr_ab = "hole_tr_ab_" + std::to_string(idx_val_band);
        std::string name_rate_elec_lo_em = "elec_lo_em_" + std::to_string(idx_val_band);
        std::string name_rate_elec_lo_ab = "elec_lo_ab_" + std::to_string(idx_val_band);
        std::string name_rate_elec_tr_em = "elec_tr_em_" + std::to_string(idx_val_band);
        std::string name_rate_elec_tr_ab = "elec_tr_ab_" + std::to_string(idx_val_band);

        int data_tag_hole_lo_em = gmsh::view::add(name_rate_hole_lo_em);
        int data_tag_hole_lo_ab = gmsh::view::add(name_rate_hole_lo_ab);
        int data_tag_hole_tr_em = gmsh::view::add(name_rate_hole_tr_em);
        int data_tag_hole_tr_ab = gmsh::view::add(name_rate_hole_tr_ab);
        int data_tag_elec_lo_em = gmsh::view::add(name_rate_elec_lo_em);
        int data_tag_elec_lo_ab = gmsh::view::add(name_rate_elec_lo_ab);
        int data_tag_elec_tr_em = gmsh::view::add(name_rate_elec_tr_em);
        int data_tag_elec_tr_ab = gmsh::view::add(name_rate_elec_tr_ab);

        gmsh::view::addHomogeneousModelData(data_tag_hole_lo_em, 0, model_file_name, "NodeData", node_tags, rates_hole_lo_em);
        gmsh::view::addHomogeneousModelData(data_tag_hole_lo_ab, 0, model_file_name, "NodeData", node_tags, rates_hole_lo_ab);
        gmsh::view::addHomogeneousModelData(data_tag_hole_tr_em, 0, model_file_name, "NodeData", node_tags, rates_hole_tr_em);
        gmsh::view::addHomogeneousModelData(data_tag_hole_tr_ab, 0, model_file_name, "NodeData", node_tags, rates_hole_tr_ab);
        gmsh::view::addHomogeneousModelData(data_tag_elec_lo_em, 0, model_file_name, "NodeData", node_tags, rates_elec_lo_em);
        gmsh::view::addHomogeneousModelData(data_tag_elec_lo_ab, 0, model_file_name, "NodeData", node_tags, rates_elec_lo_ab);
        gmsh::view::addHomogeneousModelData(data_tag_elec_tr_em, 0, model_file_name, "NodeData", node_tags, rates_elec_tr_em);
        gmsh::view::addHomogeneousModelData(data_tag_elec_tr_ab, 0, model_file_name, "NodeData", node_tags, rates_elec_tr_ab);

        gmsh::view::write(data_tag_hole_lo_em, final_filename, true);
        gmsh::view::write(data_tag_hole_lo_ab, final_filename, true);
        gmsh::view::write(data_tag_hole_tr_em, final_filename, true);
        gmsh::view::write(data_tag_hole_tr_ab, final_filename, true);
        gmsh::view::write(data_tag_elec_lo_em, final_filename, true);
        gmsh::view::write(data_tag_elec_lo_ab, final_filename, true);
        gmsh::view::write(data_tag_elec_tr_em, final_filename, true);
        gmsh::view::write(data_tag_elec_tr_ab, final_filename, true);
    }

    // for (int index_band = 0; index_band < number_bands; ++index_band) {
    //     std::string         band_name = "band_" + std::to_string(index_band);
    //     std::vector<double> current_band_values(m_node_tags.size());
    //     for (std::size_t index_node = 0; index_node < m_node_tags.size(); ++index_node) {
    //         current_band_values[index_node] = band_values[index_node * number_bands + index_band];
    //         // std::cout << "band_values[" << index_node << "]: " << band_values[index_node * number_bands + index_band] << std::endl;
    //     }
    //     int data_tag = gmsh::view::add(band_name);
    //     if (m_node_tags.size() != current_band_values.size()) {
    //         std::cout << "number of nodes: " << m_node_tags.size() << std::endl;
    //         std::cout << "number of values: " << current_band_values.size() << std::endl;
    //         throw std::runtime_error("Number of nodes and number of values are not the same. Abort.");
    //     }
    //     gmsh::view::addHomogeneousModelData(data_tag, 0, model_file_name, "NodeData", m_node_tags, current_band_values);
    //     const int   index_view             = gmsh::view::getIndex(data_tag);
    //     std::string name_object_visibility = "View[" + std::to_string(index_view) + "].Visible";
    //     gmsh::option::setNumber(name_object_visibility, 0);
    //     gmsh::view::write(data_tag, out_filename, true);
    // }
}

void ElectronPhonon::load_phonon_parameters(const std::string& filename) {
    YAML::Node config = YAML::LoadFile(filename);

    // Check if the file is empty
    if (config.IsNull()) {
        throw std::runtime_error("File " + filename + " is empty");
    }
    // Print the file
    std::cout << "File " << filename << " contains: " << std::endl;
    std::cout << config << std::endl;

    auto list_materials = config["materials"];

    const std::string& my_material = m_material.get_name();

    auto same_material = [my_material](const YAML::Node& node) { return node["name"].as<std::string>() == my_material; };
    auto it_material   = std::find_if(list_materials.begin(), list_materials.end(), same_material);

    if (it_material == list_materials.end()) {
        throw std::runtime_error("Material " + my_material + " not found in file " + filename);
    }
    auto material = *it_material;

    double radiusWS      = material["Radius-WS"].as<double>();
    m_radii_wigner_seitz = radiusWS;

    // Example for reading nested dispersion properties
    auto dispersion = material["dispersion"];
    for (const auto& type : {"longitudinal", "transverse"}) {
        auto dispType = dispersion[type];
        for (const auto& wave : {"acoustic", "optic"}) {
            auto   waveType = dispType[wave];
            double w0       = waveType["w0"].as<double>();
            double vs       = waveType["vs"].as<double>();
            double c        = waveType["c"].as<double>();
            std::cout << "w0: " << w0 << " vs: " << vs << " c: " << c << std::endl;

            PhononDirection  direction = (type == "longitudinal") ? PhononDirection::longitudinal : PhononDirection::transverse;
            PhononMode       mode      = (wave == "acoustic") ? PhononMode::acoustic : PhononMode::optical;
            PhononDispersion phononDispersion(mode, direction, w0, vs, c);
            m_phonon_dispersion[std::make_pair(mode, direction)] = phononDispersion;
        }
    }

    // Example for reading deformation potential
    auto node_deformationPotential = material["deformation-potential"];
    for (const auto& carrierType : {"electron", "hole"}) {
        auto carrier = node_deformationPotential[carrierType];
        for (const auto& wave : {"acoustic", "optic"}) {
            auto   waveType = carrier[wave];
            double A        = waveType["A"].as<double>();
            double B        = waveType["B"].as<double>();
            std::cout << "A: " << A << " B: " << B << std::endl;

            PhononMode           mode = (wave == "acoustic") ? PhononMode::acoustic : PhononMode::optical;
            DeformationPotential deformationPotential(mode, A, B);
            if (carrierType == "electron") {
                if (mode == PhononMode::acoustic) {
                    m_acoustic_deformation_potential_e = deformationPotential;
                } else {
                    m_optical_deformation_potential_e = deformationPotential;
                }
            } else if (carrierType == "hole") {
                if (mode == PhononMode::acoustic) {
                    m_acoustic_deformation_potential_h = deformationPotential;
                } else {
                    m_optical_deformation_potential_h = deformationPotential;
                }
            } else {
                throw std::runtime_error("Unknown carrier type.");
            }
        }
    }
}

void ElectronPhonon::export_rate_values(const std::string& filename) const {
    std::ofstream file(filename);
    for (auto&& vertex : m_list_vertices) {
        std::vector<std::array<double, 8>> all_rates = vertex.get_electron_phonon_rates_all_bands();
        for (std::size_t idx_band = 0; idx_band < all_rates.size(); ++idx_band) {
            double energy = vertex.get_energy_at_band(idx_band);
            std::cout << "Energy at band " << idx_band << ": " << energy << std::endl;
            file << energy << " ";
            for (auto&& rate : all_rates[idx_band]) {
                file << rate << " ";
            }
            file << std::endl;
        }
    }
    file.close();
}

}  // namespace bz_mesh