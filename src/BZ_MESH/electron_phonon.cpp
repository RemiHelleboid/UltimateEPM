/**
 * @file elelectron_phonon.cpp
 * @brief Electron–phonon implementation (matches refactored header: m_phonon_dispersion[4], RateValues::add)
 * @date 2024-02-09
 */

#include "electron_phonon.hpp"

#include <Eigen/Dense>
#include <array>
#include <atomic>
#include <cassert>
#include <cmath>
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
#include "gmsh.h"
#include "omp.h"
#include "physical_constants.hpp"
#include "physical_functions.hpp"
#include "vector.hpp"
#include "yaml-cpp/yaml.h"

namespace uepm::mesh_bz {

void ElectronPhonon::clean_all_elph_data() {
    m_phonon_rates_transport.clear();
    m_rates_nk_npkp.clear();
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
    return max_energy * uepm::Constants::h_bar_eV;  // ħω → eV
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
                                                                    std::size_t idx_tetra_final,
                                                                    bool        push_nk_npkp) {
    Rate8 out{};

    const auto& vtx1  = m_list_vertices[idx_k1];
    const auto& tetra = m_list_tetrahedra[idx_tetra_final];

    const double   Ei_eV = vtx1.get_energy_at_band(idx_n1);
    const vector3& k1    = vtx1.get_position();
    const vector3& k2    = tetra.compute_barycenter();

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
    const double     pi                 = uepm::Constants::pi;
    const double     qe                 = uepm::Constants::q_e;
    const double     hbar_eV            = uepm::Constants::h_bar_eV;

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
        const double omega = disp.omega_lookup(qn);
        if (omega <= SMALL_OMEGA_CUTOFF) {
            continue;
        }

        const double Eph_eV = hbar_eV * omega;
        const double N0     = bose_einstein_distribution(Eph_eV, m_temperature_K);

        // Deformation potential (J)
        const DeformationPotential& defpot  = (mode == PhononMode::acoustic) ? m_ac_defpot_e : m_op_defpot_e;
        const double                Delta_J = defpot.get_fischetti_deformation_potential(q, idx_n1) * qe;

        const double pref = m_spin_degeneracy * (pi / (m_rho_kg_m3 * omega)) * (Delta_J * Delta_J) * I2 / m_reduce_bz_factor;

        // --- Emission (Ef = Ei - ħω), bose = N0 + 1 ---
        {
            const double Ef_eV = Ei_eV - Eph_eV;
            if (tetra.is_energy_inside_band(Ef_eV, idx_n2)) {
                // const double dos_eV = tetra.interpolate_dos_at_energy_per_band(Ef_eV, idx_n2);
                const double dos_eV = tetra.compute_tetra_dos_energy_band(Ef_eV, idx_n2);
                if (dos_eV > 0.0) {
                    const double val      = pref * (N0 + 1.0) * (dos_eV / qe);
                    const int    mode_idx = rate_index(mode, dir, PhononEvent::emission);
                    out[mode_idx] += val;
                    inv_mrta_rate += val * transport_weight_value;  // for 1/τ_tr
                }
            }
        }
        // --- Absorption (Ef = Ei + ħω), bose = N0 ---
        {
            const double Ef_eV = Ei_eV + Eph_eV;
            if (tetra.is_energy_inside_band(Ef_eV, idx_n2)) {
                // const double dos_eV = tetra.interpolate_dos_at_energy_per_band(Ef_eV, idx_n2);
                const double dos_eV = tetra.compute_tetra_dos_energy_band(Ef_eV, idx_n2);
                if (dos_eV > 0.0) {
                    const double val      = pref * (N0) * (dos_eV / qe);
                    const int    mode_idx = rate_index(mode, dir, PhononEvent::absorption);
                    out[mode_idx] += val;
                    inv_mrta_rate += val * transport_weight_value;  // for 1/τ_tr
                }
            }
        }
    }
    std::size_t local_n1 = get_local_band_index(idx_n1);
    m_phonon_rates_transport[local_n1][idx_k1] += inv_mrta_rate;  // accumulate 1/τ_tr(E) on uniform grid
    return out;
}

/**
 * @brief Compute the out-scattering electron-phonon rate for a given state (n1,k1).
 *
 * @param idx_n1 Index of the initial band.
 * @param idx_k1 Index of the initial k-point.
 * @param populate_nk_npkp Whether to populate the (n,k) → (n',k') transition rates.
 * @return RateValues The computed out-scattering rates.
 */
RateValues ElectronPhonon::compute_electron_phonon_rate(std::size_t idx_n1, std::size_t idx_k1, bool populate_nk_npkp) {
    RateValues acc{};

    const double     Ei_eV      = m_list_vertices[idx_k1].get_energy_at_band(idx_n1);
    const double     Eph_max    = get_max_phonon_energy();
    constexpr double eps_energy = 1e-5;
    const double     Ef_min     = Ei_eV - Eph_max - eps_energy;
    const double     Ef_max     = Ei_eV + Eph_max + eps_energy;

    for (auto idx_n2 : get_band_indices(MeshParticleType::conduction)) {
        // Quick reject band window
        if (Ef_min > m_max_band[idx_n2] || Ef_max < m_min_band[idx_n2]) {
            continue;
        }

        for (std::size_t t = 0; t < m_list_tetrahedra.size(); ++t) {
            const auto& tetra = m_list_tetrahedra[t];
            if (!tetra.does_intersect_band_energy_range(Ef_min, Ef_max, idx_n2)) {
                continue;
            }

            const Rate8 r = compute_electron_phonon_transition_rates_pair(idx_n1, idx_k1, idx_n2, t, populate_nk_npkp);
            for (int i = 0; i < 8; ++i) {
                acc.v[i] += r[i];
                if (r[i] > 1e12) {
                    std::cout << std::scientific << std::setprecision(3) << r[i] << std::defaultfloat << std::endl;
                }
            }
        }
    }

    return acc;
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

            const double q_norm = q.norm();

            for (int md = 0; md < 4; ++md) {
                const auto&           disp = m_phonon_dispersion[md];
                const PhononMode      mode = (md < 2) ? PhononMode::acoustic : PhononMode::optical;
                const PhononDirection dir  = (md & 1) ? PhononDirection::transverse : PhononDirection::longitudinal;

                const double omega = disp.omega_lookup(q_norm);
                if (omega <= SMALL_OMEGA_CUTOFF) {
                    continue;
                }

                const double Eph_eV = uepm::Constants::h_bar_eV * omega;
                const double N0     = bose_einstein_distribution(Eph_eV, m_temperature_K);

                const DeformationPotential& defpot  = (mode == PhononMode::acoustic) ? m_ac_defpot_h : m_op_defpot_h;
                const double                Delta_J = defpot.get_fischetti_deformation_potential(q, idx_n1) * uepm::Constants::q_e;

                // Emission
                {
                    const double Ef_eV  = Ei_eV - Eph_eV;
                    const double dos_eV = tetra.interpolate_dos_at_energy_per_band(Ef_eV, static_cast<std::size_t>(idx_n2));
                    if (dos_eV > 0.0) {
                        const double dos_per_J = dos_eV / uepm::Constants::q_e;
                        double       rate_value =
                            (uepm::Constants::pi / (m_rho_kg_m3 * omega)) * (Delta_J * Delta_J) * overlap2 * (N0 + 1.0) * dos_per_J;
                        rate_value /= m_reduce_bz_factor;
                        rate_value *= m_spin_degeneracy;

                        rates_k1_n1.add(mode, dir, PhononEvent::emission, rate_value);
                    }
                }
                // Absorption
                {
                    const double Ef_eV  = Ei_eV + Eph_eV;
                    const double dos_eV = tetra.interpolate_dos_at_energy_per_band(Ef_eV, static_cast<std::size_t>(idx_n2));
                    if (dos_eV > 0.0) {
                        const double dos_per_J = dos_eV / uepm::Constants::q_e;
                        double rate_value = (uepm::Constants::pi / (m_rho_kg_m3 * omega)) * (Delta_J * Delta_J) * overlap2 * (N0)*dos_per_J;
                        rate_value /= m_reduce_bz_factor;
                        rate_value *= m_spin_degeneracy;

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
 * @param populate_nk_npkp Whether to populate the (n,k) → (n',k') transition rates.
 */
void ElectronPhonon::compute_electron_phonon_rates_over_mesh(double energy_max, bool irreducible_wedge_only, bool populate_nk_npkp) {
    std::cout << "Progress: 0%";
    std::atomic<std::size_t> counter{0};
    constexpr int            chunk_size = 32;

    m_phonon_rates_transport.clear();
    m_phonon_rates_transport.resize(get_number_bands(MeshParticleType::conduction));
    std::cout << "Allocating transport rates for " << m_phonon_rates_transport.size() << " conduction bands and " << m_list_vertices.size()
              << " k-points.\n";
    for (auto& vec : m_phonon_rates_transport) {
        vec.resize(m_list_vertices.size(), 0.0);
    }
    // Create a shuffled list of indices to balance load when using irreducible wedge only
    std::vector<std::size_t> random_indices(m_list_vertices.size());
    for (std::size_t i = 0; i < m_list_vertices.size(); ++i) {
        random_indices[i] = i;
    }

    std::size_t nb_vertices_to_compute = irreducible_wedge_only ? m_list_vtx_in_iwedge.size() : m_list_vertices.size();
    std::cout << "Computing electron-phonon rates for " << nb_vertices_to_compute << " k-points"
              << (irreducible_wedge_only ? " (irreducible wedge only)" : "") << ".\n";
    std::cout << "Using " << m_nb_threads << " threads.\n";

    auto list_bands = get_band_indices(MeshParticleType::conduction);
    std::cout << "Computing electron-phonon rates for conduction bands: ";
    for (auto b : list_bands) {
        std::cout << b << " ";
    }
    std::cout << "\n";
    std::size_t              skipped_bc_energy = 0;
    std::atomic<std::size_t> done{0};
    std::size_t              total = nb_vertices_to_compute;
    std::size_t              step  = std::max<std::size_t>(1, total / 100);  // ~1% steps

#pragma omp parallel for schedule(dynamic, chunk_size) num_threads(m_nb_threads) reduction(+ : skipped_bc_energy)
    for (std::size_t idx = 0; idx < nb_vertices_to_compute; ++idx) {
        const std::size_t idx_k1 = irreducible_wedge_only ? m_list_vtx_in_iwedge[idx] : idx;
        Vertex&           k1     = m_list_vertices[idx_k1];
        for (auto idx_n1 : get_band_indices(MeshParticleType::conduction)) {
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
#pragma omp critical(cout)
            {
                std::cout << "\rDone " << d << "/" << total << " (" << std::fixed << std::setprecision(1)
                          << (100.0 * double(d) / double(total)) << "%)" << std::flush;
            }
        }
    }
    std::cout << std::defaultfloat;

    const std::size_t total_states = nb_vertices_to_compute * list_bands.size();
    std::cout << "Skipped " << skipped_bc_energy << " k-points above " << energy_max << " eV : "
              << "(" << std::fixed << std::setprecision(1) << (100.0 * skipped_bc_energy / total_states) << "% of total states)\n";

    if (irreducible_wedge_only) {
        done  = 0;
        total = m_list_vtx_in_iwedge.size();
        step  = std::max<std::size_t>(1, total / 100);
        std::cout << "Set electron-phonon rates for all mesh vertices.\n";
#pragma omp parallel for schedule(dynamic) num_threads(m_nb_threads)
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
#pragma omp critical(cout)
                {
                    std::cout << "\rDone " << d << "/" << total << " (" << std::fixed << std::setprecision(1)
                              << (100.0 * double(d) / double(total)) << "%)" << std::flush;
                }
            }
        }
        std::cout << "\nSet rates for all k-points." << std::endl;
    }
}

/**
 * @brief Select the final state for an electron-phonon interaction.
 *
 * @param idx_band_initial The initial band index.
 * @param k_initial The initial k-point.
 * @param mode The phonon mode.
 * @param direction The phonon direction.
 * @param event The phonon event.
 * @param rng The random number generator.
 * @return std::pair<int, std::size_t> The selected final band index and k-point index.
 */
std::pair<int, std::size_t> ElectronPhonon::select_electron_phonon_final_state(std::size_t     idx_band_initial,
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

    // Initial state
    Tetra* init_tetra = find_tetra_at_location(k_initial);
    if (!init_tetra) {
        throw std::runtime_error("select_final_state: initial k not inside any tetrahedron.");
    }
    const double Ei_eV = init_tetra->interpolate_energy_at_band(k_initial, idx_band_initial);

    const double sign_ph = (event == PhononEvent::emission) ? -1.0 : +1.0;

    const size_t        nb_bands = m_list_vertices.front().get_number_bands();
    const size_t        nb_tetra = m_list_tetrahedra.size();
    std::vector<double> probs_flat(nb_bands * nb_tetra, 0.0);
    auto                P_ref = [&](int n2, size_t t) -> double& { return probs_flat[static_cast<size_t>(n2) * nb_tetra + t]; };

    const double pi      = uepm::Constants::pi;
    const double qe      = uepm::Constants::q_e;
    const double hbar_eV = uepm::Constants::h_bar_eV;

    const double Eph_max_eV = get_max_phonon_energy();

    for (std::size_t n2 = 0; n2 < nb_bands; ++n2) {
        // Band window
        const double Ef_min = Ei_eV - Eph_max_eV;
        const double Ef_max = Ei_eV + Eph_max_eV;
        if (Ef_min > m_max_band[n2] || Ef_max < m_min_band[n2]) {
            continue;
        }

        for (size_t t = 0; t < nb_tetra; ++t) {
            const auto& tetra = m_list_tetrahedra[t];
            if (!tetra.does_intersect_band_energy_range(Ef_min, Ef_max, n2)) {
                continue;
            }

            const vector3 k2_bary = tetra.compute_barycenter();

            vector3 q = k2_bary - k_initial;
            if (!is_inside_mesh_geometry(q)) {
                q = retrieve_k_inside_mesh_geometry(q);
            }
            if (!is_inside_mesh_geometry(q)) {
                continue;
            }

            const double qn    = q.norm();
            const double omega = disp.omega_lookup(qn);
            if (!(omega > 0.0) || omega < 1e-12) {
                continue;
            }

            const double Eph_eV = hbar_eV * omega;
            const double N0     = bose_einstein_distribution(Eph_eV, m_temperature_K);
            const double bose   = (sign_ph < 0.0) ? (N0 + 1.0) : N0;

            const double Ef_eV  = Ei_eV + sign_ph * Eph_eV;
            const double dos_eV = tetra.interpolate_dos_at_energy_per_band(Ef_eV, static_cast<std::size_t>(n2));
            if (dos_eV <= 0.0) {
                continue;
            }

            const double dos_per_J = dos_eV / qe;

            const double I  = electron_overlap_integral(k_initial, k2_bary, m_radius_wigner_seitz_m);
            const double I2 = I * I;

            const DeformationPotential& defpot  = (mode == PhononMode::acoustic) ? m_ac_defpot_e : m_op_defpot_e;
            const double                Delta_J = defpot.get_fischetti_deformation_potential(q, static_cast<int>(idx_band_initial)) * qe;

            double P = (pi / (m_rho_kg_m3 * omega)) * (Delta_J * Delta_J) * I2 * bose * dos_per_J;
            P /= m_reduce_bz_factor;
            P *= m_spin_degeneracy;

            if (P > 0.0 && std::isfinite(P)) {
                P_ref(n2, t) = P;
            }
        }
    }

    double total = 0.0;
    for (double p : probs_flat) {
        total += p;
    }
    if (!(total > 0.0) || !std::isfinite(total)) {
        throw std::runtime_error("select_final_state: no admissible final states (total probability = 0).");
    }

    std::uniform_real_distribution<double> U(0.0, 1.0);
    const double                           threshold = U(rng) * total;

    double acc = 0.0;
    for (size_t flat = 0; flat < probs_flat.size(); ++flat) {
        acc += probs_flat[flat];
        if (acc >= threshold) {
            const int    n2 = static_cast<int>(flat / nb_tetra);
            const size_t t  = static_cast<size_t>(flat % nb_tetra);
            return {n2, t};  // band, tetra index (k' = barycenter(t))
        }
    }

    throw std::runtime_error("select_final_state: internal sampling error.");
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
    file.precision(8);

    const std::size_t nb_elph = get_nb_bands_elph();
    for (const auto& vertex : m_list_vertices) {
        const auto all_rates = vertex.get_electron_phonon_rates_all_bands();
        if (all_rates.size() != nb_elph) {
            throw std::runtime_error("export_rate_values: size mismatch in rates");
        }

        for (std::size_t local = 0; local < nb_elph; ++local) {
            const std::size_t g = get_global_band_index(local, m_elph_particle_type);  // map to global
            const double      E = vertex.get_energy_at_band(g);
            file << g << ',' << E;
            const auto& r = all_rates[local];
            for (double v : r) {
                file << ',' << v;
            }
            file << '\n';
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
void ElectronPhonon::compute_plot_electron_phonon_rates_vs_energy_over_mesh(int                nb_bands,
                                                                            double             max_energy,
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

    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Cannot open " + filename + " for writing.");
    }

    static constexpr std::array<const char*, 8>
        mode_phonons{"ac_L_ab", "ac_T_ab", "op_L_ab", "op_T_ab", "ac_L_em", "ac_T_em", "op_L_em", "op_T_em"};

    out << "energy(eV),dos(ev^-1m^-3)";
    for (auto* mode : mode_phonons) {
        out << ',' << mode;
    }
    out << '\n';

    const std::size_t n_steps = static_cast<std::size_t>(std::floor(max_energy / energy_step)) + 1;
    std::cout << std::setprecision(8) << "Computing electron-phonon rates vs energy over mesh..." << std::endl;
    for (std::size_t istep = 0; istep < n_steps; ++istep) {
        const double E = std::min(max_energy, istep * energy_step);
        std::cout << "\rEnergy: " << E << " / " << max_energy << std::flush;

        double                dos_sum = 0.0;
        std::array<double, 8> num{};  // accumulators
        num.fill(0.0);
        const std::size_t nb_bands_elph = get_nb_bands_elph();
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

        std::array<double, 8> mean{};
        if (dos_sum > 0.0) {
            const double inv_dos = 1.0 / dos_sum;
            for (int i = 0; i < 8; ++i) {
                mean[i] = num[i] * inv_dos;
            }
        } else {
            mean.fill(0.0);
        }

        out << std::scientific << std::setprecision(10) << E << ',' << dos_sum;
        for (double v : mean) {
            out << ',' << v;
        }
        out << '\n';
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
    const auto& vertex_indices = tetra->get_index_vertices_with_sorted_energy_at_band(idx_band);

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
            const double e_ph = disp.omega_lookup(q.norm()) * uepm::Constants::h_bar_eV;
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

    for (int idx_band = 0; idx_band < get_number_bands_total(); ++idx_band) {
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
            rates_ac_lo_em[idx_k1]  = rates[0];
            rates_ac_lo_ab[idx_k1]  = rates[1];
            rates_ac_tr_em[idx_k1]  = rates[2];
            rates_ac_tr_ab[idx_k1]  = rates[3];
            rates_opt_lo_em[idx_k1] = rates[4];
            rates_opt_lo_ab[idx_k1] = rates[5];
            rates_opt_tr_em[idx_k1] = rates[6];
            rates_opt_tr_ab[idx_k1] = rates[7];
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
    YAML::Node config = YAML::LoadFile(filename);
    if (config.IsNull()) {
        throw std::runtime_error("File " + filename + " is empty");
    }

    // std::cout << "File " << filename << " contains:\n" << config << std::endl;

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
            double           q_max_norm = 1.5 / m_si2red;
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
}

/**
 * @brief Read phonon scattering rates from a CSV file and populate the internal data structure.
 *
 * The CSV file is expected to have the following format:
 * band_index, energy (eV), ac_lo_em, ac_lo_ab, ac_tr_em, ac_tr_ab, op_lo_em, op_lo_ab, op_tr_em, op_tr_ab
 *
 * @param path The path to the CSV file.
 */
void ElectronPhonon::read_phonon_scattering_rates_from_file(const std::filesystem::path& path) {
    std::cout << "Reading phonon scattering rates (CSV) from file " << path.string() << " ...\n";

    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open file " + path.string());
    }

    m_list_phonon_scattering_rates.clear();
    m_list_phonon_scattering_rates.resize(m_list_vertices.size());

    std::string line;
    std::size_t line_no = 0;

    for (std::size_t idx_vtx = 0; idx_vtx < m_list_vertices.size(); ++idx_vtx) {
        const auto&       vertex   = m_list_vertices[idx_vtx];
        const std::size_t nbands   = vertex.get_number_bands();
        auto&             per_band = m_list_phonon_scattering_rates[idx_vtx];
        per_band.resize(nbands);

        for (std::size_t idx_band = 0; idx_band < nbands; ++idx_band) {
            do {
                if (!std::getline(in, line)) {
                    std::cout << "Line no: " << line_no << " " << line << std::endl;
                    throw std::runtime_error("Unexpected EOF at vertex " + std::to_string(idx_vtx) + ", band " + std::to_string(idx_band));
                }
                ++line_no;
            } while (line.empty() || line[0] == '#' || line[0] == ';');

            for (char& c : line) {
                if (c == ',') {
                    c = ' ';
                }
            }
            std::istringstream iss(line);

            std::size_t band_idx_file{};
            double      energy_file{};
            Rate8       rates{};

            // Keep your CSV layout if you still emit ALO..ETA in this order elsewhere
            if (!(iss >> band_idx_file >> energy_file >> rates[0] >> rates[1] >> rates[2] >> rates[3] >> rates[4] >> rates[5] >> rates[6] >>
                  rates[7])) {
                throw std::runtime_error("Malformed CSV line " + std::to_string(line_no));
            }
            if (band_idx_file == nbands) {
                // More bands in the file than in the mesh: ignore extra bands
                break;
            }
            per_band[idx_band] = rates;
        }
    }
    in.close();

    std::cout << "Finished reading phonon scattering rates for " << m_list_vertices.size() << " vertices.\n";
}

/**
 * @brief Compute the maximum total phonon scattering rate (P_Gamma) over all vertices and bands.
 *
 * @return The maximum P_Gamma value.
 */
double ElectronPhonon::compute_P_Gamma() const {
    double pgamma_max = 0.0;
    for (const auto& perVertex : m_list_phonon_scattering_rates) {
        for (const auto& rate8 : perVertex) {
            const double tot = std::accumulate(rate8.begin(), rate8.end(), 0.0);
            if (tot > pgamma_max) {
                pgamma_max = tot;
            }
        }
    }
    return pgamma_max;
}
// -------------------- TRANSPORT ---------------------
// elelectron_phonon.cpp  (near the bottom, before namespace close)
#include "physical_functions.hpp"  // fermi_dirac_distribution, d_de_fermi_dirac_dE

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
    const double inv_2pi3 = 1.0 / std::pow(2.0 * uepm::Constants::pi, 3);  // avoid M_PI
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
        const double scale = get_reduce_bz_factor() * m_spin_degeneracy;
        for (double& w : m_count_weight_tetra_per_vertex) {
            w *= scale;
        }
    }

    // 2) σ = q^2 ∑_{n,k} w_k τ_tr(n,k) (-df0/dE) v v^T  with  v=(1/ħ)∇_k E
    //    n = ∑_{n∈cond,k} w_k f0(E)
    const double q  = uepm::Constants::q_e;       // Coulomb
    const double hE = uepm::Constants::h_bar_eV;  // eV·s

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
            const double dfdE = -d_de_fermi_dirac_dE(E, fermi_level_eV, temperature_K) / uepm::Constants::q_e;  // [1/J]

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

}  // namespace uepm::mesh_bz
