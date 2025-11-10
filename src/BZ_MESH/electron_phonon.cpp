/**
 * @file elelectron_phonon.cpp
 * @brief Electron–phonon implementation (matches refactored header: m_phonon_dispersion[4], RateValues::add)
 * @date 2024-02-09
 */

#include "electron_phonon.hpp"

#include <csv.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <gmsh.h>
#include <omp.h>

#include <Eigen/Dense>
#include <array>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <tuple>
#include <vector>

#include "bz_states.hpp"
#include "physical_constants.hpp"
#include "physical_functions.hpp"
#include "vector.hpp"
#include "yaml-cpp/yaml.h"

namespace uepm::mesh_bz {
/**
 * @brief Scale the q-norm into the first BZ using a Wigner-Seitz scheme.
 * (Unused currently.)
 *
 * @param q_norm
 * @return double
 */
double ElectronPhonon::scale_q_norm(double q_norm) const {
    const double qmax = 1.0 / si_to_reduced_scale();
    double       qm   = q_norm / qmax;
    if (qm > 1.0) {
        if (qm <= 2.0) {
            qm = 2.0 - qm;
        } else {
            qm -= 2.0;
        }
    }
    qm *= qmax;
    return qm;
}

double PGamma::interpolate_P_Gamma(double energy_eV) const {
    if (m_energies_eV.empty() || m_P_Gamma_values.empty()) {
        throw std::runtime_error("interpolate_P_Gamma: empty data.");
    }
    if (energy_eV <= m_energies_eV.front()) {
        return m_P_Gamma_values.front();
    }
    if (energy_eV >= m_energies_eV.back()) {
        return m_P_Gamma_values.back();
    }
    // Linear interpolation
    auto        it  = std::upper_bound(m_energies_eV.begin(), m_energies_eV.end(), energy_eV);
    std::size_t idx = std::distance(m_energies_eV.begin(), it);
    double      t   = (energy_eV - m_energies_eV[idx - 1]) / (m_energies_eV[idx] - m_energies_eV[idx - 1]);
    return (1.0 - t) * m_P_Gamma_values[idx - 1] + t * m_P_Gamma_values[idx];
}

void ElectronPhonon::clean_all_elph_data() {
    m_phonon_rates_transport.clear();
    m_list_phonon_scattering_rates.clear();
    m_count_weight_tetra_per_vertex.clear();
}

double ElectronPhonon::get_max_phonon_energy() const {
    double max_energy = std::numeric_limits<double>::lowest();
    for (const auto& disp : m_phonon_dispersion) {
        double max_w = disp.max_omega();  // ω_max [1/s]
        if (max_w > max_energy) {
            max_energy = max_w;
        }
    }
    return max_energy * uepm::constants::h_bar_eV;
}

void ElectronPhonon::export_phonon_dispersion(const std::string& filename) const {
    for (const auto& disp : m_phonon_dispersion) {
        std::string   mode_str = (disp.mode == PhononMode::acoustic) ? "acoustic" : "optical";
        std::string   dir_str  = (disp.direction == PhononDirection::longitudinal) ? "longitudinal" : "transverse";
        std::string   file_out = fmt::format("{}_{}_dispersion.csv", mode_str, dir_str);
        std::ofstream disp_file(file_out);
        if (!disp_file.is_open()) {
            throw std::runtime_error(fmt::format("export_phonon_dispersion: cannot open file {}", file_out));
        }
        disp_file << "q (1/m),omega (1/s),omega (eV)\n";
        const double      hbar_eV = uepm::constants::h_bar_eV;
        const std::size_t N       = disp.N;
        for (std::size_t i = 0; i < N; ++i) {
            const double q        = disp.q0 + i * (1.0 / disp.inv_dq);
            const double omega    = disp.omega_samples[i];
            const double omega_eV = hbar_eV * omega;
            disp_file << fmt::format("{:.6e},{:.6e},{:.6e}\n", q, omega, omega_eV);
        }
        disp_file.close();
        fmt::print("Exported phonon dispersion to {}\n", file_out);
    }
}

/**
 * @brief Compute the electron-phonon transition rates for a pair of states.
 *
 * @param idx_n1 Index of the initial band.
 * @param idx_k1 Index of the initial k-point.
 * @param idx_n2 Index of the final band.
 * @param idx_tetra_final Index of the final tetrahedron.
 * @param push_nk_npkp Whether to push the (n,k) → (n',k') transition rates.
 * @return Rate8 The computed transition rates.
 * @param push_nk_npkp
 * @return Rate8
 */
Rate8 ElectronPhonon::compute_electron_phonon_transition_rates_pair(std::size_t idx_n1,
                                                                    std::size_t idx_k1,
                                                                    std::size_t idx_n2,
                                                                    std::size_t idx_tetra_final) {
    Rate8 rates_n1k1_n2kT{};

    const auto& vtx1  = m_list_vertices[idx_k1];
    const auto& tetra = m_list_tetrahedra[idx_tetra_final];

    const double   Ei_eV = vtx1.get_energy_at_band(idx_n1);
    const vector3& k1    = vtx1.get_position();
    const vector3& k2    = tetra.compute_barycenter();

    std::size_t local_n1 = get_local_band_index(idx_n1);
    std::size_t local_n2 = get_local_band_index(idx_n2);

    // Overlap integral |I|^2
    const double I  = electron_overlap_integral(k1, k2, m_radius_wigner_seitz_m);
    const double I2 = I * I;

    // q = k2 - k1
    vector3 q = k2 - k1;

    if (!is_inside_mesh_geometry(q)) {
        q = retrieve_k_inside_mesh_geometry(q);
    }
    if (!is_inside_mesh_geometry(q)) {
        throw std::runtime_error("q not in BZ");
    }

    const double qn = q.norm();

    constexpr double SMALL_OMEGA_CUTOFF = 1.0;  // [1/s]
    const double     pi                 = uepm::constants::pi;
    const double     qe                 = uepm::constants::q_e;
    const double     hbar_eV            = uepm::constants::h_bar_eV;

    double        inv_mrta_rate          = 0.0;
    const vector3 vnk                    = vtx1.get_energy_gradient_at_band(idx_n1) * (1.0 / hbar_eV);   // m/s
    const vector3 v_npkp                 = tetra.get_gradient_energy_at_band(idx_n2) * (1.0 / hbar_eV);  // m/s
    const double  transport_weight_value = transport_weight_RTA(vnk, v_npkp);

    // Loop 4 branches: md=0..3 → (ac/op)×(L/T)
    for (int md = 0; md < 4; ++md) {
        const auto&           disp = m_phonon_dispersion[md];
        const PhononMode      mode = ((md >> 1) == 0) ? PhononMode::acoustic : PhononMode::optical;
        const PhononDirection dir  = ((md & 1) == 0) ? PhononDirection::longitudinal : PhononDirection::transverse;

        // ω(|q|) [1/s] — use your lookup or analytic
        const double omega = disp.omega_analytic(qn);
        if (omega <= SMALL_OMEGA_CUTOFF) {
            continue;
        }

        const double Eph_eV = hbar_eV * omega;
        const double N0     = bose_einstein_distribution(Eph_eV, m_temperature_K);

        // Deformation potential (J)
        const DeformationPotential& defpot  = (mode == PhononMode::acoustic) ? m_ac_defpot_e : m_op_defpot_e;
        double                      Delta_J = defpot.get_fischetti_deformation_potential(qn, local_n1) * qe;

        const double pref = (pi / (m_rho_kg_m3 * omega)) * (Delta_J * Delta_J) * I2;

        // --- Emission (Ef = Ei - ħω), bose = N0 + 1 ---
        {
            const double Ef_eV = Ei_eV - Eph_eV;

            // // draw k' on iso-triangle at Ef_sel_eV
            // vector3 k_final = Tsel.draw_random_uniform_point_at_energy(Ef_sel_eV, static_cast<std::size_t>(n2_sel), rng);

            // const double interpolated_Ef = Tsel.interpolate_energy_at_band(k_final, static_cast<std::size_t>(n2_sel));
            // if (std::abs(interpolated_Ef - Ef_sel_eV) > 1e-9) {
            //     throw std::runtime_error("select_final_state: sampled k_final energy mismatch.");
            // }
            // // Return final state data

            if (tetra.is_energy_inside_band(Ef_eV, idx_n2)) {
                const double dos_eV = tetra.compute_tetra_dos_energy_band(Ef_eV, idx_n2);
                if (dos_eV > 0.0) {
                    double    val      = pref * (N0 + 1.0) * (dos_eV / qe);
                    const int mode_idx = rate_index(mode, dir, PhononEvent::emission);
                    rates_n1k1_n2kT[mode_idx] += val;
                    inv_mrta_rate += val * transport_weight_value;  // for 1/τ_tr
                }
            }
        }
        // --- Absorption (Ef = Ei + ħω), bose = N0 ---
        {
            const double Ef_eV = Ei_eV + Eph_eV;
            if (tetra.is_energy_inside_band(Ef_eV, idx_n2)) {
                const double dos_eV = tetra.compute_tetra_dos_energy_band(Ef_eV, idx_n2);
                if (dos_eV > 0.0) {
                    double    val      = pref * (N0) * (dos_eV / qe);
                    const int mode_idx = rate_index(mode, dir, PhononEvent::absorption);
                    rates_n1k1_n2kT[mode_idx] += val;
                    inv_mrta_rate += val * transport_weight_value;  // for 1/τ_tr
                }
            }
        }
    }
    m_phonon_rates_transport[local_n1][idx_k1] += inv_mrta_rate;  // 1/tau_transport
    return rates_n1k1_n2kT;
}

/**
 * @brief Compute the out-scattering electron-phonon rate for a given state (n1,k1).
 *
 * @param idx_n1 Index of the initial band.
 * @param idx_k1 Index of the initial k-point.
 * @return RateValues The computed out-scattering rates.
 */
RateValues ElectronPhonon::compute_electron_phonon_rate(std::size_t idx_n1, std::size_t idx_k1) {
    RateValues       total_out_rates_n1k1{};
    std::size_t      count_nnz  = 0;
    const double     Ei_eV      = m_list_vertices[idx_k1].get_energy_at_band(idx_n1);
    const double     Eph_max    = get_max_phonon_energy();
    constexpr double eps_energy = 1e-2;
    const double     Ef_min     = Ei_eV - Eph_max - eps_energy;
    const double     Ef_max     = Ei_eV + Eph_max + eps_energy;

    auto list_bands_n2 = get_band_indices(MeshParticleType::conduction);

    std::size_t count_bands_n2_spaned = 0;
    std::size_t total_spaned          = 0;

    constexpr double threshold = 1e1;

    for (auto idx_n2 : list_bands_n2) {
        // Quick reject band window
        if (Ef_min > m_max_band[idx_n2] || Ef_max < m_min_band[idx_n2]) {
            continue;
        }
        count_bands_n2_spaned++;
        const auto& ordered_tetra_indices = m_tetra_ordered_energy_min[idx_n2];
        total_spaned += ordered_tetra_indices.candidate_indices(Ef_min, Ef_max).size();

        for (auto idx_tetra : ordered_tetra_indices.candidate_indices(Ef_min, Ef_max)) {
            const auto& tetra = m_list_tetrahedra[idx_tetra];

            if (!tetra.does_intersect_band_energy_range(Ef_min, Ef_max, idx_n2)) {
                continue;
            }

            const Rate8 rates_n1k1_n2kT = compute_electron_phonon_transition_rates_pair(idx_n1, idx_k1, idx_n2, idx_tetra);
            for (int ph_branch = 0; ph_branch < 8; ++ph_branch) {
                total_out_rates_n1k1.m_rate_values[ph_branch] += rates_n1k1_n2kT[ph_branch];
                if (rates_n1k1_n2kT[ph_branch] > threshold) {
                    count_nnz++;
                }
            }
        }
    }
    // DEBUG
    double ratio_nnz = static_cast<double>(count_nnz) / 8.0 / total_spaned * 100.0;
    // fmt::print("compute_electron_phonon_rate: (n1={}, k1={}) nonzero transitions: {} ({:.2}%)\n", idx_n1, idx_k1, count_nnz, ratio_nnz);
    return total_out_rates_n1k1;
}

/**
 * @brief Compute the hole-phonon transition rates for a given state (n1,k1).
 *
 * @param idx_n1 Index of the initial band.
 * @param idx_k1 Index of the initial k-point.
 * @return RateValues The computed hole-phonon rates.
 */
RateValues ElectronPhonon::compute_hole_phonon_rate(std::size_t idx_n1, std::size_t idx_k1) {
    RateValues  rates_k1_n1;
    const auto& list_tetrahedra = m_list_tetrahedra;

    const double  Ei_eV = m_list_vertices[idx_k1].get_energy_at_band(idx_n1);
    const vector3 k1    = m_list_vertices[idx_k1].get_position();

    constexpr double SMALL_OMEGA_CUTOFF = 1.0;

    for (auto idx_n2 : get_band_indices(MeshParticleType::valence)) {
        for (const auto& tetra : list_tetrahedra) {
            const vector3 k2 = tetra.compute_barycenter();

            const double overlap  = hole_overlap_integral(idx_n1, k1, idx_n2, k2, m_hole_overlap_int_params);
            const double overlap2 = overlap * overlap;

            vector3 q = k2 - k1;
            if (!is_inside_mesh_geometry(q)) {
                q = retrieve_k_inside_mesh_geometry(q);
            }
            if (!is_inside_mesh_geometry(q)) {
                throw std::runtime_error("q not in BZ");
            }
            double q_norm = q.norm();
            for (int md = 0; md < 4; ++md) {
                const auto&           disp = m_phonon_dispersion[md];
                const PhononMode      mode = (md < 2) ? PhononMode::acoustic : PhononMode::optical;
                const PhononDirection dir  = (md & 1) ? PhononDirection::transverse : PhononDirection::longitudinal;

                const double omega = disp.omega_analytic(q_norm);
                if (omega <= SMALL_OMEGA_CUTOFF) {
                    continue;
                }

                const double Eph_eV = uepm::constants::h_bar_eV * omega;
                const double N0     = bose_einstein_distribution(Eph_eV, m_temperature_K);

                const DeformationPotential& defpot  = (mode == PhononMode::acoustic) ? m_ac_defpot_h : m_op_defpot_h;
                const double                Delta_J = defpot.get_fischetti_deformation_potential(q_norm, idx_n1) * uepm::constants::q_e;

                // Emission
                {
                    const double Ef_eV  = Ei_eV - Eph_eV;
                    const double dos_eV = tetra.interpolate_dos_at_energy_per_band(Ef_eV, static_cast<std::size_t>(idx_n2));
                    if (dos_eV > 0.0) {
                        const double dos_per_J = dos_eV / uepm::constants::q_e;
                        double       rate_value =
                            (uepm::constants::pi / (m_rho_kg_m3 * omega)) * (Delta_J * Delta_J) * overlap2 * (N0 + 1.0) * dos_per_J;
                        rates_k1_n1.add(mode, dir, PhononEvent::emission, rate_value);
                    }
                }
                // Absorption
                {
                    const double Ef_eV  = Ei_eV + Eph_eV;
                    const double dos_eV = tetra.interpolate_dos_at_energy_per_band(Ef_eV, static_cast<std::size_t>(idx_n2));
                    if (dos_eV > 0.0) {
                        const double dos_per_J = dos_eV / uepm::constants::q_e;
                        double rate_value = (uepm::constants::pi / (m_rho_kg_m3 * omega)) * (Delta_J * Delta_J) * overlap2 * (N0)*dos_per_J;
                        rates_k1_n1.add(mode, dir, PhononEvent::absorption, rate_value);
                    }
                }
            }  // md
        }  // tetra
    }  // bands

    return rates_k1_n1;
}

/**
 * @brief Compute the out-scattering electron-phonon rates over the entire k-point mesh.
 *
 * @param energy_max Maximum energy (eV) to compute rates for.
 * @param irreducible_wedge_only Whether to use only the irreducible wedge of the BZ.
 */
void ElectronPhonon::compute_electron_phonon_rates_over_mesh(double energy_max, bool irreducible_wedge_only) {
    fmt::print("Computing electron-phonon rates over mesh...\n");
    fmt::print("  Energy max: {:.2f} eV\n", energy_max);
    fmt::print("  Irreducible wedge only: {}\n", irreducible_wedge_only);
    auto start_time = std::chrono::high_resolution_clock::now();

    m_phonon_rates_transport.clear();
    m_phonon_rates_transport.resize(get_number_bands(MeshParticleType::conduction));
    fmt::print("Allocating transport rates for {} conduction bands and {} k-points.\n",
               m_phonon_rates_transport.size(),
               m_list_vertices.size());
    for (auto& vec : m_phonon_rates_transport) {
        vec.resize(m_list_vertices.size(), 0.0);
    }

    std::size_t nb_vertices_to_compute = irreducible_wedge_only ? m_list_vtx_in_iwedge.size() : m_list_vertices.size();
    fmt::print("Computing electron-phonon rates for {} k-points{}\n",
               nb_vertices_to_compute,
               irreducible_wedge_only ? " (irreducible wedge only)" : "");
    fmt::print("Using {} threads.\n", m_nb_threads_mesh_ops);

    auto list_bands = get_band_indices(MeshParticleType::conduction);
    fmt::print("Computing electron-phonon rates for conduction bands: ");
    for (auto b : list_bands) {
        fmt::print("{} ", b);
    }
    fmt::print("\n");
    std::size_t              skipped_bc_energy = 0;
    std::atomic<std::size_t> done{0};
    std::size_t              total = nb_vertices_to_compute;
    std::size_t              step  = std::max<std::size_t>(1, total / 1000);  // ~1% steps
    fmt::print("Starting computation...\n");

    const std::size_t chunk_size = 32;
#pragma omp parallel for schedule(dynamic) num_threads(m_nb_threads_mesh_ops) reduction(+ : skipped_bc_energy)
    for (std::size_t idx = 0; idx < nb_vertices_to_compute; ++idx) {
        const std::size_t idx_k1 = irreducible_wedge_only ? m_list_vtx_in_iwedge[idx] : idx;
        Vertex&           k1     = m_list_vertices[idx_k1];
        for (auto idx_n1 : list_bands) {
            if (k1.get_energy_at_band(idx_n1) > energy_max) {
                ++skipped_bc_energy;
                k1.add_electron_phonon_rates(std::array<double, 8>{0, 0, 0, 0, 0, 0, 0, 0});
                continue;
            } else {
                auto rate = compute_electron_phonon_rate(idx_n1, idx_k1);
                k1.add_electron_phonon_rates(rate.as_array());
            }
        }
        std::size_t d = ++done;
        if (d % step == 0 || d == total) {
#pragma omp critical
            {
                fmt::print("\rDone {}/{} ({:.1f}%)", d, total, 100.0 * double(d) / double(total));
                std::cout.flush();
            }
        }
    }
    fmt::print("\n");

    const std::size_t total_states = nb_vertices_to_compute * list_bands.size();
    fmt::print("Skipped {} k-points above {} eV : ({:.1f}% of total states)\n",
               skipped_bc_energy,
               energy_max,
               100.0 * skipped_bc_energy / total_states);

    if (irreducible_wedge_only) {
        done  = 0;
        total = m_list_vtx_in_iwedge.size();
        fmt::print("Set electron-phonon rates for all mesh vertices.\n");
#pragma omp parallel for schedule(static, chunk_size) num_threads(m_nb_threads_mesh_ops)
        for (std::size_t idx_iw : m_list_vtx_in_iwedge) {
            const auto& vtx           = m_list_vertices[idx_iw];
            const auto& equiv_indices = get_all_equivalent_indices_in_bz(vtx.get_position());
            for (auto idx_k1 : equiv_indices) {
                if (idx_k1 == idx_iw) {
                    continue;
                }
                for (auto idx_n1 : get_band_indices(MeshParticleType::conduction)) {
                    std::size_t local_idx_n1 = get_local_band_index(idx_n1);
                    auto        rates_symm   = vtx.get_electron_phonon_rates(local_idx_n1);
                    m_list_vertices[idx_k1].add_electron_phonon_rates(rates_symm);
                    m_phonon_rates_transport[local_idx_n1][idx_k1] = m_phonon_rates_transport[local_idx_n1][idx_iw];
                }
            }
            std::size_t d = ++done;
            if (d % step == 0 || d == total) {
#pragma omp critical
                {
                    fmt::print("\rDone {}/{} ({:.1f}%)", d, total, 100.0 * double(d) / double(total));
                    std::cout.flush();
                }
            }
        }
        fmt::print("\nSet rates for all k-points.\n\n");
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = end_time - start_time;
    fmt::print("Completed electron-phonon rates computation in {:.2f} seconds.\n\n", std::chrono::duration<double>(duration).count());
}

SelectedFinalState ElectronPhonon::select_electron_phonon_final_state(std::size_t     idx_band_initial,
                                                                      const vector3&  k_initial,
                                                                      PhononMode      mode,
                                                                      PhononDirection direction,
                                                                      PhononEvent     event,
                                                                      std::mt19937&   rng) const {
    const int md = md_index(mode, direction);
    if (md < 0) {
        throw std::runtime_error("select_final_state: invalid mode/direction.");
    }
    const auto& disp = m_phonon_dispersion[md];

    // Initial tetra & energy
    Tetra* init_tetra = find_tetra_at_location(k_initial);
    if (!init_tetra) {
        throw std::runtime_error("select_final_state: initial k not inside any tetrahedron.");
    }
    const double Ei_eV = init_tetra->interpolate_energy_at_band(k_initial, idx_band_initial);

    const double      sign_ph  = (event == PhononEvent::emission) ? -1.0 : +1.0;
    const std::size_t nb_bands = m_list_vertices.empty() ? 0 : m_list_vertices.front().get_number_bands();
    if (nb_bands == 0) {
        return SelectedFinalState{idx_band_initial, init_tetra->get_index(), init_tetra, k_initial, Ei_eV};
    }

    using Cand = std::pair<std::size_t, std::size_t>;  // (band n2, global tetra index)
    std::vector<Cand>   candidates;
    std::vector<double> probs;

    const double pi      = uepm::constants::pi;
    const double qe      = uepm::constants::q_e;
    const double hbar_eV = uepm::constants::h_bar_eV;

    constexpr double eps_energy = 1e-5;
    const double     Eph_max_eV = get_max_phonon_energy() + eps_energy;

    const double Ef_min_win = Ei_eV - Eph_max_eV;
    const double Ef_max_win = Ei_eV + Eph_max_eV;

    std::size_t local_idx_band_initial = get_local_band_index(idx_band_initial);

    for (std::size_t n2 = 0; n2 < nb_bands; ++n2) {
        // quick band-window reject
        if (Ef_min_win > m_max_band[n2] || Ef_max_win < m_min_band[n2]) {
            continue;
        }

        const auto& list_idx_relevent_tetra = get_ordered_tetra_indices_at_band(n2);

        // for (std::size_t idx_tetra = 0; idx_tetra < list_idx_relevent_tetra.size(); ++idx_tetra) {
        //     const std::size_t idx_tetra   = list_idx_relevent_tetra[idx_tetra];  // global tetra index
        for (auto idx_tetra : list_idx_relevent_tetra) {
            const auto& final_tetra = m_list_tetrahedra[idx_tetra];

            if (!final_tetra.does_intersect_band_energy_range(Ef_min_win, Ef_max_win, n2)) {
                continue;
            }
            // The tetra list is ordered by increasing min energy, so we can break early.
            if (final_tetra.get_min_energy_at_band(n2) > Ef_max_win) {
                break;
            }

            // representative q* → ω*, Ef*
            const vector3 k2_bary = final_tetra.compute_barycenter();
            vector3       qstar   = k2_bary - k_initial;

            if (!is_inside_mesh_geometry(qstar)) {
                qstar = retrieve_k_inside_mesh_geometry(qstar);
            }
            if (!is_inside_mesh_geometry(qstar)) {
                continue;
            }

            const double omega_star = disp.omega_analytic(qstar.norm());
            if (!(omega_star > 0.0)) {
                continue;
            }

            const double Eph_star_eV = hbar_eV * omega_star;
            const double Ef_eV       = Ei_eV + sign_ph * Eph_star_eV;

            if (!final_tetra.is_energy_inside_band(Ef_eV, n2)) {
                continue;
            }

            // DOS(Ef)
            const double dos_eV = final_tetra.compute_tetra_dos_energy_band(Ef_eV, n2);
            if (!(dos_eV > 0.0)) {
                continue;
            }

            // centroid of iso-triangle for k-dependent parts
            const auto iso_pts = final_tetra.compute_band_iso_energy_surface(Ef_eV, n2);
            if (iso_pts.size() < 3) {
                continue;
            }
            const vector3 k2_centroid = (iso_pts[0] + iso_pts[1] + iso_pts[2]) / 3.0;

            vector3 qloc = k2_centroid - k_initial;
            if (!is_inside_mesh_geometry(qloc)) {
                qloc = retrieve_k_inside_mesh_geometry(qloc);
                // std::cout << "Qn : " << qloc.norm() << "\n";
            }
            if (!is_inside_mesh_geometry(qloc)) {
                continue;
            }
            double qn = qloc.norm();

            const double omega_loc = disp.omega_analytic(qloc.norm());
            if (!(omega_loc > 0.0)) {
                continue;
            }

            const double Eph_loc_eV = hbar_eV * omega_loc;
            const double N0         = bose_einstein_distribution(Eph_loc_eV, m_temperature_K);
            const double bose       = (event == PhononEvent::emission) ? (N0 + 1.0) : N0;

            const double dos_per_J = dos_eV / qe;

            const double I  = electron_overlap_integral(k_initial, k2_centroid, m_radius_wigner_seitz_m);
            const double I2 = I * I;

            const DeformationPotential& defpot  = (mode == PhononMode::acoustic) ? m_ac_defpot_e : m_op_defpot_e;
            double                      Delta_J = defpot.get_fischetti_deformation_potential(qn, local_idx_band_initial) * qe;

            // fit factor
            double factor_rates = (mode == PhononMode::acoustic) ? m_fit_acoustic : m_fit_optical;

            double P = (pi / (m_rho_kg_m3 * omega_loc)) * (Delta_J * Delta_J) * I2 * bose * dos_per_J;
            P *= factor_rates;

            if (P > 0.0 && std::isfinite(P)) {
                candidates.emplace_back(n2, idx_tetra);  // keep the **global** tetra index
                probs.push_back(P);
            }
        }
    }

    // --- sample from probs ---
    double total = std::accumulate(probs.begin(), probs.end(), 0.0);

    if (!(total > 0.0) || !std::isfinite(total)) {
        // no admissible final state: stay put
        return SelectedFinalState{idx_band_initial, init_tetra->get_index(), init_tetra, k_initial, Ei_eV};
    }

    std::uniform_real_distribution<double> U(0.0, 1.0);
    const double                           threshold = U(rng) * total;

    double      acc  = 0.0;
    std::size_t pick = 0;
    for (; pick < probs.size(); ++pick) {
        acc += probs[pick];
        if (acc >= threshold) {
            break;
        }
    }
    if (pick >= probs.size()) {
        pick = probs.size() - 1;
    }

    const std::size_t n2_sel  = candidates[pick].first;
    const std::size_t tet_sel = candidates[pick].second;  // **global** tetra index
    const auto&       Tsel    = m_list_tetrahedra[tet_sel];

    // rebuild Ef using same recipe as scoring
    const vector3 k2_bary_sel = Tsel.compute_barycenter();

    vector3 qstar_sel = k2_bary_sel - k_initial;

    if (!is_inside_mesh_geometry(qstar_sel)) {
        qstar_sel = retrieve_k_inside_mesh_geometry(qstar_sel);
    }
    const double omega_star_sel  = disp.omega_analytic(qstar_sel.norm());
    const double Eph_star_sel_eV = hbar_eV * omega_star_sel;
    const double Ef_sel_eV       = Ei_eV + sign_ph * Eph_star_sel_eV;

    // draw k' on iso-triangle at Ef_sel_eV
    vector3 k_final = Tsel.draw_random_uniform_point_at_energy(Ef_sel_eV, static_cast<std::size_t>(n2_sel), rng);

    const double interpolated_Ef = Tsel.interpolate_energy_at_band(k_final, static_cast<std::size_t>(n2_sel));
    if (std::abs(interpolated_Ef - Ef_sel_eV) > 1e-9) {
        throw std::runtime_error("select_final_state: sampled k_final energy mismatch.");
    }
    if (!Tsel.is_location_inside(k_final)) {
        throw std::runtime_error("select_final_state: sampled k_final not inside tetrahedron.");
    }

    return SelectedFinalState{n2_sel, tet_sel, const_cast<Tetra*>(&Tsel), k_final, Ef_sel_eV};
}

SelectedFinalState ElectronPhonon::select_electron_phonon_final_state(std::size_t    idx_band_initial,
                                                                      const vector3& k_initial,
                                                                      int            idx_phonon_branch,
                                                                      std::mt19937&  rng) const {
    PhononScatteringEvent PhBranch = inverse_rate_index(idx_phonon_branch);
    return select_electron_phonon_final_state(idx_band_initial, k_initial, PhBranch.mode, PhBranch.direction, PhBranch.event, rng);
}

/**
 * @brief Export the computed electron-phonon rates to a CSV file.
 *
 * @param filename The output CSV filename.
 */
void ElectronPhonon::export_rate_values(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("export_rate_values: cannot open '" + filename + "' for writing");
    }
    file.setf(std::ios::scientific);
    file.precision(6);
    // Header
    file << "vertex_index,local_band_index,energy_eV";
    file << ",rate_ac_L_ab,rate_ac_T_ab,rate_op_L_ab,rate_op_T_ab";
    file << ",rate_ac_L_em,rate_ac_T_em,rate_op_L_em,rate_op_T_em\n";
    const std::size_t nb_band_elph = get_nb_bands_elph();
    for (const auto& vertex : m_list_vertices) {
        const auto all_rates = vertex.get_electron_phonon_rates_all_bands();
        if (all_rates.size() != nb_band_elph) {
            throw std::runtime_error("export_rate_values: size mismatch in rates");
        }

        for (std::size_t local = 0; local < nb_band_elph; ++local) {
            const std::size_t idx_vtx         = vertex.get_index();
            std::size_t       global_band_idx = get_global_band_index(local, m_elph_particle_type);
            const double      E               = vertex.get_energy_at_band(global_band_idx);
            const auto&       r               = all_rates[local];
            // Only export if at least one rate is above threshold
            constexpr double rates_threshold = 1e-12;  // s^-1
            if (std::find_if(r.begin(), r.end(), [rates_threshold](double v) { return v > rates_threshold; }) != r.end()) {
                file << idx_vtx << ',' << local << ',' << E;
                for (double v : r) {
                    file << ',' << v;
                }
                file << '\n';
            }
        }
    }
}

/**
 * @brief Compute and plot the electron-phonon rates as a function of energy over the mesh.
 *
 * @param nb_bands The number of bands to consider.
 * @param max_energy The maximum energy to consider.
 * @param energy_step The energy step size.
 * @param filename The output filename for the plot.
 */
void ElectronPhonon::compute_plot_electron_phonon_rates_vs_energy_over_mesh(double             max_energy,
                                                                            double             energy_step,
                                                                            const std::string& filename) {
    if (energy_step <= 0.0) {
        throw std::invalid_argument("energy_step must be > 0");
    }
    if (max_energy < 0.0) {
        throw std::invalid_argument("max_energy must be >= 0");
    }

    if (m_list_vertices.empty()) {
        throw std::runtime_error("No vertices in mesh.");
    }

    static constexpr std::array<const char*, 8>
        mode_phonons{"ac_L_ab", "ac_T_ab", "op_L_ab", "op_T_ab", "ac_L_em", "ac_T_em", "op_L_em", "op_T_em"};

    // out << "energy(eV),dos(ev^-1m^-3)";

    const std::size_t nb_bands_elph = get_nb_bands_elph();
    // Get the energy minimum of the el-ph bands
    double energy_min = std::numeric_limits<double>::max();
    for (std::size_t idx_band = 0; idx_band < nb_bands_elph; ++idx_band) {
        std::size_t gband = get_global_band_index(idx_band, m_elph_particle_type);
        if (m_min_band[gband] < energy_min) {
            energy_min = m_min_band[gband];
        }
    }
    max_energy += energy_min;  // shift to min energy of el-ph bands

    const std::size_t n_steps = static_cast<std::size_t>(std::floor((max_energy - energy_min) / energy_step)) + 1;
    fmt::print("Computing electron-phonon rates vs energy over mesh...\n");
    std::vector<double>                energies(n_steps);
    std::vector<double>                dos_values(n_steps, 0.0);
    std::vector<std::array<double, 8>> rate_values(n_steps, std::array<double, 8>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

    std::atomic<std::size_t> done{0};
    std::size_t              total = n_steps;
    std::size_t              step  = std::max<std::size_t>(1, total / 100);  // ~1% steps
#pragma omp parallel for schedule(static) num_threads(m_nb_threads_mesh_ops)
    for (std::size_t istep = 0; istep < n_steps; ++istep) {
        const double E                = energy_step * static_cast<double>(istep) + energy_min;
        energies[istep]               = energy_step * static_cast<double>(istep);
        double                dos_sum = 0.0;
        std::array<double, 8> num{};
        num.fill(0.0);
        for (const auto& tetra : m_list_tetrahedra) {
            for (std::size_t idx_band = 0; idx_band < nb_bands_elph; ++idx_band) {
                std::size_t  idx_global_band = get_global_band_index(idx_band, m_elph_particle_type);
                const double dos_t           = tetra.compute_tetra_dos_energy_band(E, idx_global_band);
                if (!std::isfinite(dos_t)) {
                    throw std::runtime_error("DOS is NaN/Inf at E=" + std::to_string(E) + " band=" + std::to_string(idx_global_band));
                }
                if (dos_t <= 0.0) {
                    continue;
                }

                dos_sum += dos_t;

                const std::array<double, 8> rates = tetra.get_tetra_electron_phonon_rates(idx_band);
                for (int i = 0; i < 8; ++i) {
                    num[i] += rates[i] * dos_t;
                }
            }
        }
        std::size_t d = ++done;
        if (d % step == 0 || d == total) {
#pragma omp critical
            {
                fmt::print("\rEnergy: {:.3f} / {:.3f} eV -> DOS: {:.2e} ev^-1m^-3", E, max_energy, dos_sum);
                std::cout.flush();
            }
        }
        dos_values[istep] = dos_sum;

        std::array<double, 8> mean{};
        if (dos_sum > 0.0) {
            const double inv_dos = 1.0 / dos_sum;
            for (int i = 0; i < 8; ++i) {
                mean[i] = num[i] * inv_dos;
            }
        } else {
            mean.fill(0.0);
        }
        rate_values[istep] = mean;
    }

    fmt::print("\nWriting data to {}...\n", filename);

    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Cannot open " + filename + " for writing.");
    }

    fmt::print(out, "energy(eV),dos(ev^-1m^-3)");
    for (auto* mode : mode_phonons) {
        fmt::print(out, ",rate_{}", mode);
    }
    fmt::print(out, "\n");
    for (std::size_t istep = 0; istep < n_steps; ++istep) {
        fmt::print(out, "{:.6e},{:.6e}", energies[istep], dos_values[istep]);
        for (int i = 0; i < 8; ++i) {
            fmt::print(out, ",{:.6e}", rate_values[istep][i]);
        }
        fmt::print(out, "\n");
    }

    m_P_Gamma_data.m_energies_eV = energies;
    m_P_Gamma_data.m_P_Gamma_values.resize(energies.size());
    for (std::size_t i = 0; i < energies.size(); ++i) {
        m_P_Gamma_data.m_P_Gamma_values[i] = std::accumulate(rate_values[i].begin(), rate_values[i].end(), 0.0);
    }

    std::cout << std::endl;
    out.close();
}

/**
 * @brief Interpolate the phonon scattering rate at a given location for a given band.
 *
 * @param location The location in k-space.
 * @param idx_band The band index.
 * @return Rate8
 */
Rate8 ElectronPhonon::interpolate_phonon_scattering_rate_at_location(const vector3& location, const std::size_t& idx_band) const {
    // Find the tetrahedron containing the location
    const Tetra* tetra = find_tetra_at_location(location);
    if (!tetra) {
        throw std::runtime_error("Location is not inside any tetrahedron");
    }

    // Get the vertex indices of the tetrahedron
    const auto& vertex_indices = tetra->get_list_indices_vertices();

    // Interpolate the scattering rates at the vertices
    Rate8 rates;
    for (int idx_mode = 0; idx_mode < 8; ++idx_mode) {
        std::vector<double> vertex_rates(4);
        for (auto idx_vtx = 0; idx_vtx < 4; ++idx_vtx) {
            const auto& vertex    = m_list_vertices[vertex_indices[idx_vtx]];
            vertex_rates[idx_vtx] = vertex.get_electron_phonon_rates(idx_band)[idx_mode];
        }
        rates[idx_mode] = tetra->interpolate_scalar_at_position(location, vertex_rates);
    }
    return rates;
}

/**
 * @brief Plot the phonon dispersion to a file.
 *
 * @param filename The output filename.
 */
void ElectronPhonon::plot_phonon_dispersion(const std::string& filename) const {
    std::ofstream file(filename);
    for (auto&& vtx : m_list_vertices) {
        auto k = vtx.get_position();
        file << k.x() << " " << k.y() << " " << k.z() << " ";
        for (int md = 0; md < 4; ++md) {
            auto q = k;
            q /= m_material.get_fourier_factor();
            const auto&  disp = m_phonon_dispersion[md];
            const double e_ph = disp.omega_analytic(q.norm()) * uepm::constants::h_bar_eV;
            file << e_ph << " ";
        }
        file << '\n';
    }
}

/**
 * @brief Add the computed electron-phonon rates to a Gmsh mesh file.
 *
 * @param initial_filename The input Gmsh mesh filename.
 * @param final_filename The output Gmsh mesh filename with added rates.
 */
void ElectronPhonon::add_electron_phonon_rates_to_mesh(const std::string& initial_filename, const std::string& final_filename) {
    // If the file exists, remove it to avoid appending to an old file
    if (std::ifstream(final_filename)) {
        std::remove(final_filename.c_str());
    }

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

    for (int idx_band = 0; idx_band < m_nb_bands_elph; ++idx_band) {
        std::vector<double> rates_ac_lo_em(m_list_vertices.size());
        std::vector<double> rates_ac_lo_ab(m_list_vertices.size());
        std::vector<double> rates_ac_tr_em(m_list_vertices.size());
        std::vector<double> rates_ac_tr_ab(m_list_vertices.size());
        std::vector<double> rates_opt_lo_em(m_list_vertices.size());
        std::vector<double> rates_opt_lo_ab(m_list_vertices.size());
        std::vector<double> rates_opt_tr_em(m_list_vertices.size());
        std::vector<double> rates_opt_tr_ab(m_list_vertices.size());

        for (std::size_t idx_k1 = 0; idx_k1 < m_list_vertices.size(); ++idx_k1) {
            auto rates              = m_list_vertices[idx_k1].get_electron_phonon_rates(static_cast<std::size_t>(idx_band));
            rates_ac_lo_ab[idx_k1]  = rates[0];
            rates_ac_tr_ab[idx_k1]  = rates[1];
            rates_opt_lo_ab[idx_k1] = rates[2];
            rates_opt_tr_ab[idx_k1] = rates[3];
            rates_ac_lo_em[idx_k1]  = rates[4];
            rates_ac_tr_em[idx_k1]  = rates[5];
            rates_opt_lo_em[idx_k1] = rates[6];
            rates_opt_tr_em[idx_k1] = rates[7];
        }

        std::string name_rate_ac_lo_em  = "ac_lo_em_" + std::to_string(idx_band);
        std::string name_rate_ac_lo_ab  = "ac_lo_ab_" + std::to_string(idx_band);
        std::string name_rate_ac_tr_em  = "ac_tr_em_" + std::to_string(idx_band);
        std::string name_rate_ac_tr_ab  = "ac_tr_ab_" + std::to_string(idx_band);
        std::string name_rate_opt_lo_em = "opt_lo_em_" + std::to_string(idx_band);
        std::string name_rate_opt_lo_ab = "opt_lo_ab_" + std::to_string(idx_band);
        std::string name_rate_opt_tr_em = "opt_tr_em_" + std::to_string(idx_band);
        std::string name_rate_opt_tr_ab = "opt_tr_ab_" + std::to_string(idx_band);

        int data_tag_ac_lo_em  = gmsh::view::add(name_rate_ac_lo_em);
        int data_tag_ac_lo_ab  = gmsh::view::add(name_rate_ac_lo_ab);
        int data_tag_ac_tr_em  = gmsh::view::add(name_rate_ac_tr_em);
        int data_tag_ac_tr_ab  = gmsh::view::add(name_rate_ac_tr_ab);
        int data_tag_opt_lo_em = gmsh::view::add(name_rate_opt_lo_em);
        int data_tag_opt_lo_ab = gmsh::view::add(name_rate_opt_lo_ab);
        int data_tag_opt_tr_em = gmsh::view::add(name_rate_opt_tr_em);
        int data_tag_opt_tr_ab = gmsh::view::add(name_rate_opt_tr_ab);

        gmsh::view::addHomogeneousModelData(data_tag_ac_lo_em, 0, model_file_name, "NodeData", node_tags, rates_ac_lo_em);
        gmsh::view::addHomogeneousModelData(data_tag_ac_lo_ab, 0, model_file_name, "NodeData", node_tags, rates_ac_lo_ab);
        gmsh::view::addHomogeneousModelData(data_tag_ac_tr_em, 0, model_file_name, "NodeData", node_tags, rates_ac_tr_em);
        gmsh::view::addHomogeneousModelData(data_tag_ac_tr_ab, 0, model_file_name, "NodeData", node_tags, rates_ac_tr_ab);
        gmsh::view::addHomogeneousModelData(data_tag_opt_lo_em, 0, model_file_name, "NodeData", node_tags, rates_opt_lo_em);
        gmsh::view::addHomogeneousModelData(data_tag_opt_lo_ab, 0, model_file_name, "NodeData", node_tags, rates_opt_lo_ab);
        gmsh::view::addHomogeneousModelData(data_tag_opt_tr_em, 0, model_file_name, "NodeData", node_tags, rates_opt_tr_em);
        gmsh::view::addHomogeneousModelData(data_tag_opt_tr_ab, 0, model_file_name, "NodeData", node_tags, rates_opt_tr_ab);

        gmsh::option::setNumber("PostProcessing.SaveMesh", 1);  // Save mesh only once
        gmsh::view::write(data_tag_ac_lo_em, final_filename, true);
        gmsh::option::setNumber("PostProcessing.SaveMesh", 0);
        gmsh::view::write(data_tag_ac_lo_ab, final_filename, true);
        gmsh::view::write(data_tag_ac_tr_em, final_filename, true);
        gmsh::view::write(data_tag_ac_tr_ab, final_filename, true);
        gmsh::view::write(data_tag_opt_lo_em, final_filename, true);
        gmsh::view::write(data_tag_opt_lo_ab, final_filename, true);
        gmsh::view::write(data_tag_opt_tr_em, final_filename, true);
        gmsh::view::write(data_tag_opt_tr_ab, final_filename, true);
    }
    gmsh::finalize();
}

/**
 * @brief Load phonon parameters from a YAML configuration file.
 *
 * @param filename The path to the YAML configuration file.
 */
void ElectronPhonon::load_phonon_parameters(const std::string& filename) {
    fmt::print("Loading phonon parameters from file {} ...\n", filename);

    YAML::Node config = YAML::LoadFile(filename);
    if (config.IsNull()) {
        throw std::runtime_error("File " + filename + " is empty");
    }

    auto               list_materials = config["materials"];
    const std::string& my_material    = m_material.get_name();

    auto same_material = [&](const YAML::Node& node) { return node["name"].as<std::string>() == my_material; };
    auto it_material   = std::find_if(list_materials.begin(), list_materials.end(), same_material);
    if (it_material == list_materials.end()) {
        throw std::runtime_error("Material " + my_material + " not found in file " + filename);
    }
    auto material = *it_material;

    double radiusWS         = material["Radius-WS"].as<double>();
    m_radius_wigner_seitz_m = radiusWS;

    // Dispersion → m_phonon_dispersion[md]
    auto dispersion = material["dispersion"];
    for (const auto& type : {"longitudinal", "transverse"}) {
        auto dispType = dispersion[type];
        for (const auto& wave : {"acoustic", "optic"}) {
            auto   waveType = dispType[wave];
            double w0       = waveType["w0"].as<double>();
            double vs       = waveType["vs"].as<double>();
            double c        = waveType["c"].as<double>();
            // std::cout << "w0: " << w0 << " vs: " << vs << " c: " << c << std::endl;

            PhononDirection direction = (std::string(type) == "longitudinal") ? PhononDirection::longitudinal : PhononDirection::transverse;
            PhononMode      mode      = (std::string(wave) == "acoustic") ? PhononMode::acoustic : PhononMode::optical;

            PhononDispersion phononDispersion(mode, direction, w0, vs, c);
            double           q_max_norm = 1.0 / m_si2red;
            std::size_t      points     = 200;
            // std::cout << "Max q norm in reduced units: " << q_max_norm << std::endl;
            phononDispersion.build_lookup(q_max_norm, points);

            const int md = md_index(mode, direction);
            if (md < 0) {
                throw std::runtime_error("load_phonon_parameters: bad mode/direction index.");
            }
            m_phonon_dispersion[md] = std::move(phononDispersion);
        }
    }

    // Deformation potentials
    auto node_deformationPotential = material["deformation-potential"];
    for (const auto& carrierType : {"electron", "hole"}) {
        auto   carrier          = node_deformationPotential[carrierType];
        double energy_threshold = carrier["energy-threshold"].as<double>();
        for (const auto& wave : {"acoustic", "optic"}) {
            auto   waveType = carrier[wave];
            double A        = waveType["A"].as<double>();
            double B        = waveType["B"].as<double>();

            PhononMode           mode = (std::string(wave) == "acoustic") ? PhononMode::acoustic : PhononMode::optical;
            DeformationPotential deformationPotential(mode, A, B, energy_threshold);
            if (std::string(carrierType) == "electron") {
                if (mode == PhononMode::acoustic) {
                    m_ac_defpot_e = deformationPotential;
                } else {
                    m_op_defpot_e = deformationPotential;
                }
            } else {
                if (mode == PhononMode::acoustic) {
                    m_ac_defpot_h = deformationPotential;
                } else {
                    m_op_defpot_h = deformationPotential;
                }
            }
        }
    }
    fmt::print("Finished loading phonon parameters for material {}.\n", my_material);
}

/**
 * @brief Read phonon scattering rates from a CSV file and populate the internal data structure.
 *
 * The CSV file is expected to have the following format:
 * - Columns: vertex_index, local_band_index, energy_eV, rate_ac_L_ab, rate_ac_T_ab, rate_op_L_ab, rate_op_T_ab,
 *           rate_ac_L_em, rate_ac_T_em, rate_op_L_em, rate_op_T_em
 *
 * @param path The path to the CSV file.
 */
void ElectronPhonon::read_phonon_scattering_rates_from_file(const std::filesystem::path& path) {
    fmt::print("Reading phonon scattering rates from file {} ...\n", path.string());

    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open file " + path.string());
    }

    m_list_phonon_scattering_rates.clear();
    m_list_phonon_scattering_rates.resize(m_list_vertices.size());
    for (auto& perVertex : m_list_phonon_scattering_rates) {
        perVertex.resize(m_nb_bands_elph);
    }

    constexpr std::size_t  nb_cols = 11;  // vtx_index, band_index, energy, 8 rates
    io::CSVReader<nb_cols> csv_reader(path.string());

    csv_reader.read_header(io::ignore_extra_column,
                           "vertex_index",
                           "local_band_index",
                           "energy_eV",
                           "rate_ac_L_ab",
                           "rate_ac_T_ab",
                           "rate_op_L_ab",
                           "rate_op_T_ab",
                           "rate_ac_L_em",
                           "rate_ac_T_em",
                           "rate_op_L_em",
                           "rate_op_T_em");
    std::size_t vertex_index;
    std::size_t band_index;
    double      energy_eV;
    double      rate_ac_L_ab;
    double      rate_ac_T_ab;
    double      rate_op_L_ab;
    double      rate_op_T_ab;
    double      rate_ac_L_em;
    double      rate_ac_T_em;
    double      rate_op_L_em;
    double      rate_op_T_em;

    // Initialize all vertices to have the correct number of bands. The rates are zero by default.
    for (auto&& vtx : m_list_vertices) {
        vtx.set_nb_electron_phonon_rates(m_nb_bands_elph);
    }
    m_list_phonon_scattering_rates.clear();
    m_list_phonon_scattering_rates.resize(m_list_vertices.size());
    for (auto& perVertex : m_list_phonon_scattering_rates) {
        perVertex.resize(m_nb_bands_elph);
    }
    std::vector<std::vector<double>> energies(m_list_vertices.size(), std::vector<double>(m_nb_bands_elph, 0.0));

    while (csv_reader.read_row(vertex_index,
                               band_index,
                               energy_eV,
                               rate_ac_L_ab,
                               rate_ac_T_ab,
                               rate_op_L_ab,
                               rate_op_T_ab,
                               rate_ac_L_em,
                               rate_ac_T_em,
                               rate_op_L_em,
                               rate_op_T_em)) {
        if (vertex_index >= m_list_vertices.size()) {
            throw std::runtime_error("Vertex index out of range in phonon scattering rates file.");
        }
        if (band_index >= m_nb_bands_elph) {
            throw std::runtime_error("Band index out of range in phonon scattering rates file.");
        }
        Rate8 rates = {rate_ac_L_ab, rate_ac_T_ab, rate_op_L_ab, rate_op_T_ab, rate_ac_L_em, rate_ac_T_em, rate_op_L_em, rate_op_T_em};

        // RAW FIT
        rates[2] *= m_fit_optical;
        rates[3] *= m_fit_optical;
        rates[6] *= m_fit_optical;
        rates[7] *= m_fit_optical;

        rates[0] *= m_fit_acoustic;
        rates[1] *= m_fit_acoustic;
        rates[4] *= m_fit_acoustic;
        rates[5] *= m_fit_acoustic;

        m_list_vertices[vertex_index].set_electron_phonon_rates(band_index, rates);
        m_list_phonon_scattering_rates[vertex_index][band_index] = rates;
        energies[vertex_index][band_index]                       = energy_eV;
    }
    fmt::print("Finished reading phonon scattering rates from file {}.\n", path.string());
    // CHECK: verify energies match
    for (std::size_t idx_vtx = 0; idx_vtx < m_list_vertices.size(); ++idx_vtx) {
        for (std::size_t idx_band = 0; idx_band < m_nb_bands_elph; ++idx_band) {
            double E_file = energies[idx_vtx][idx_band];
            double E_mesh = m_list_vertices[idx_vtx].get_energy_at_band(get_global_band_index(idx_band, m_elph_particle_type));
            if (std::abs(E_file) > 0) {
                if (std::abs(E_file - E_mesh) > 1e-6) {
                    throw std::runtime_error("Energy mismatch at vertex " + std::to_string(idx_vtx) + " band " + std::to_string(idx_band) +
                                             ": file=" + std::to_string(E_file) + " mesh=" + std::to_string(E_mesh));
                }
            }
        }
    }
}

/**
 * @brief Compute the maximum total phonon scattering rate (P_Gamma) over all vertices and bands.
 *
 * @return The maximum P_Gamma value.
 */
double ElectronPhonon::compute_P_Gamma() const {
    // Compute the maximum total phonon scattering rate over all vertices and bands
    double P_Gamma_max = 0.0;
    for (const auto& vertex_rates : m_list_phonon_scattering_rates) {
        for (const auto& rates_per_band : vertex_rates) {
            double total_rate = 0.0;
            for (double rate : rates_per_band) {
                total_rate += rate;
            }
            if (total_rate > P_Gamma_max) {
                P_Gamma_max = total_rate;
            }
        }
    }
    return P_Gamma_max;
}

/**
 * @brief Compute the electron mobility tensor using the MRTA approach.
 *
 * @param fermi_level_eV
 * @param temperature_K
 * @param conduction_only
 * @return Eigen::Matrix3d
 */
Eigen::Matrix3d ElectronPhonon::compute_electron_MRTA_mobility_tensor(double fermi_level_eV, double temperature_K, bool conduction_only) {
    using namespace uepm::physics;

    // Sanity checks
    if (m_list_vertices.empty() || m_list_tetrahedra.empty()) {
        throw std::runtime_error("MRTA: mesh is empty.");
    }
    if (m_phonon_rates_transport.empty()) {
        throw std::runtime_error(
            "MRTA: m_phonon_rates_transport is empty. "
            "Run compute_electron_phonon_rates_over_mesh() first.");
    }
    if (m_phonon_rates_transport.size() != m_nb_bands_elph) {
        throw std::runtime_error("MRTA: m_phonon_rates_transport size mismatch vs number of bands.");
    }
    for (const auto& v : m_phonon_rates_transport) {
        if (v.size() != m_list_vertices.size()) {
            throw std::runtime_error("MRTA: m_phonon_rates_transport[*] size mismatch vs number of k-points.");
        }
    }

    // 1) Per-vertex k-space weights (share of k-volume per vertex), scaled so that
    //    Σ_k w_k (per band) ≈ g_s / V_cell (as in your DOS normalization).
    const double inv_2pi3 = 1.0 / std::pow(2.0 * uepm::constants::pi, 3);  // avoid M_PI
    if (m_count_weight_tetra_per_vertex.size() != m_list_vertices.size()) {
        m_count_weight_tetra_per_vertex.assign(m_list_vertices.size(), 0.0);
        for (const auto& T : m_list_tetrahedra) {
            const double Vt  = std::fabs(T.get_signed_volume());
            const auto   ids = T.get_list_indices_vertices();  // 4 vertex indices
            // Equal-share to vertices (robust on fine meshes).
            for (int i = 0; i < 4; ++i) {
                m_count_weight_tetra_per_vertex[ids[i]] += Vt * 0.25 * inv_2pi3;  // [m^-3]
            }
        }
        // Scale by reduce-BZ factor and spin degeneracy
        const double scale = m_spin_degeneracy;
        for (double& w : m_count_weight_tetra_per_vertex) {
            w *= scale;
        }
    }

    // 2) σ = q^2 ∑_{n,k} w_k τ_tr(n,k) (-df0/dE) v v^T  with  v=(1/ħ)∇_k E
    //    n = ∑_{n∈cond,k} w_k f0(E)
    const double q  = uepm::constants::q_e;       // Coulomb
    const double hE = uepm::constants::h_bar_eV;  // eV·s

    Eigen::Matrix3d sigma = Eigen::Matrix3d::Zero();  // S/m
    double          n_e   = 0.0;                      // m^-3

    auto bands = get_band_indices(MeshParticleType::conduction);
    if (!conduction_only) {
        // Extend here if you later want to add holes (σ_h) too.
        // For now, keep electrons only.
    }

    for (auto b : bands) {
        const std::size_t idx_band_local = get_local_band_index(b);
        const auto&       inv_tau_at_k   = m_phonon_rates_transport[idx_band_local];  // [1/s] per vertex

        for (std::size_t k = 0; k < m_list_vertices.size(); ++k) {
            const double wk = m_count_weight_tetra_per_vertex[k];
            if (wk <= 0.0) {
                continue;
            }

            const double E    = m_list_vertices[k].get_energy_at_band(b);  // eV
            const double f0   = fermi_dirac_distribution(E, fermi_level_eV, temperature_K);
            const double dfdE = -d_de_fermi_dirac_dE(E, fermi_level_eV, temperature_K) / uepm::constants::q_e;  // [1/J]

            const double inv_tau = inv_tau_at_k[k];
            if (!(inv_tau > 0.0) || !std::isfinite(inv_tau)) {
                continue;
            }
            const double tau = 1.0 / inv_tau;  // s

            // v = (1/ħ) ∇_k E
            const vector3   gE = m_list_vertices[k].get_energy_gradient_at_band(b);  // eV/m
            Eigen::Vector3d v;
            v << (gE.x() / hE), (gE.y() / hE), (gE.z() / hE);  // m/s

            // **Include q^2** to get σ in S/m
            sigma.noalias() += (q * q) * (wk * tau * dfdE) * (v * v.transpose());

            // electron density
            n_e += wk * f0;
        }
    }

    if (!(n_e > 0.0) || !std::isfinite(n_e)) {
        throw std::runtime_error(
            "MRTA: computed carrier density n_e is zero/invalid. "
            "Check EF, T, and that rates were computed for conduction bands.");
    }

    // 3) μ = σ / (n q)  (returns m^2/(V·s))
    const double    denom = n_e * q;
    Eigen::Matrix3d mu    = sigma / denom;
    return mu;
}

double ElectronPhonon::compute_electron_MRTA_mobility_isotropic(double fermi_level_eV, double temperature_K, bool conduction_only) {
    const Eigen::Matrix3d mu = compute_electron_MRTA_mobility_tensor(fermi_level_eV, temperature_K, conduction_only);
    return mu.trace() / 3.0;  // isotropic average, m^2/(V·s)
}

double ElectronPhonon::mean_electron_energy_equilibrium(double fermi_level_eV, double temperature_K, bool excess_above_cbm) const {
    using namespace uepm::physics;

    if (m_list_vertices.empty() || m_list_tetrahedra.empty()) {
        throw std::runtime_error("Mean energy (e): mesh is empty.");
    }

    // Build per-vertex k-space weights (same recipe as MRTA, but local copy)
    const double        inv_2pi3 = 1.0 / std::pow(2.0 * uepm::constants::pi, 3);
    std::vector<double> w_k(m_list_vertices.size(), 0.0);
    for (const auto& T : m_list_tetrahedra) {
        const double Vt  = std::fabs(T.get_signed_volume());
        const auto   ids = T.get_list_indices_vertices();
        for (int i = 0; i < 4; ++i) {
            w_k[ids[i]] += Vt * 0.25 * inv_2pi3;
        }
    }
    const double scale = get_reduce_bz_factor() * m_spin_degeneracy;
    for (double& w : w_k) {
        w *= scale;
    }

    // Reference band edge (CBM) if requested
    double Ec_min = std::numeric_limits<double>::infinity();
    for (auto b : get_band_indices(MeshParticleType::conduction)) {
        Ec_min = std::min(Ec_min, m_min_band[b]);
    }

    double num = 0.0;  // eV * m^-3
    double den = 0.0;  // m^-3

    auto bands = get_band_indices(MeshParticleType::conduction);
    for (auto b : bands) {
        for (std::size_t k = 0; k < m_list_vertices.size(); ++k) {
            const double wk = w_k[k];
            if (wk <= 0.0) {
                continue;
            }
            const double E  = m_list_vertices[k].get_energy_at_band(b);  // eV
            const double f0 = fermi_dirac_distribution(E, fermi_level_eV, temperature_K);
            if (f0 <= 0.0) {
                continue;
            }

            const double Euse = excess_above_cbm ? (E - Ec_min) : E;  // eV
            num += wk * Euse * f0;
            den += wk * f0;
        }
    }

    if (!(den > 0.0) || !std::isfinite(den)) {
        throw std::runtime_error("Mean energy (e): carrier density is zero/invalid for given EF,T.");
    }
    return num / den;  // eV
}

void ElectronPhonon::test_elph() const {
    std::cout << "Running electron-phonon test: sampling rates over energy..." << std::endl;
    std::string   filename = "eelph_test_output.txt";
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + filename + " for writing.");
    }
    std::size_t Nsample = 100;
    // Uniform distribution in energy range
    std::size_t                            glob_band  = 0;
    double                                 min_energy = m_min_band[glob_band] + 0.01;
    double                                 max_energy = m_max_band[glob_band] - 0.01;
    std::mt19937                           rng(42);  // fixed seed for reproducibility
    std::uniform_real_distribution<double> energy_dist(min_energy, max_energy);
    std::vector<double>                    energies(Nsample);
    std::vector<double>                    sum_rates(Nsample, 0.0);
#pragma omp parallel for
    for (std::size_t i = 0; i < Nsample; ++i) {
        double  energy  = energy_dist(rng);
        vector3 k_point = draw_random_k_point_at_energy(energy, glob_band, rng);
        Tetra*  tetra   = find_tetra_at_location(k_point);
        if (tetra) {
            double e_interp = tetra->interpolate_energy_at_band(k_point, glob_band);
            auto   rates    = tetra->interpolate_phonon_scattering_rate_at_location(k_point, 0);
            double sum_rate = 0.0;
            sum_rate        = std::accumulate(rates.begin(), rates.end(), 0.0);
            energies[i]     = energy;
            sum_rates[i]    = sum_rate;
        }
    }
    // Write to file
    file << "energy(eV),sum_rate(s^-1)\n";
    for (std::size_t i = 0; i < Nsample; ++i) {
        file << energies[i] << "," << sum_rates[i] << "\n";
    }
    file.close();
}

}  // namespace uepm::mesh_bz
