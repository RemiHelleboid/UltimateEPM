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
#include <random>
#include <sstream>
#include <tuple>
#include <vector>

#include "Constants.hpp"
#include "Vector3D.h"
#include "bz_states.hpp"
#include "gmsh.h"
#include "omp.h"
#include "yaml-cpp/yaml.h"

namespace bz_mesh {

// ----- OpenMP custom reduction for std::array<double,8> -----
static inline void arr8_add(std::array<double, 8>& __restrict out, const std::array<double, 8>& __restrict in) noexcept {
    for (int i = 0; i < 8; ++i)
        out[i] += in[i];
}
#pragma omp declare reduction(merge : std::array<double, 8> : arr8_add(omp_out, omp_in)) initializer(omp_priv = std::array<double, 8>{})

// -------------------- Small inline helpers --------------------

inline double ElectronPhonon::bose_einstein_distribution(double energy_eV, double temperature_K) const {
    // N0 = 1 / (exp(E / kT) - 1)
    const double x = energy_eV / (EmpiricalPseudopotential::Constants::k_b_eV * temperature_K);
    return 1.0 / std::expm1(x);  // stable for small x
}

inline double ElectronPhonon::electron_overlap_integral(const vector3& k1, const vector3& k2) const {
    // Spherical j1-like: 3 (sin x - x cos x)/x^3 with x = |k1-k2| * R_WS
    constexpr double R_Wigner_Seitz = 2.122e-10;  // m
    const double     x              = (k1 - k2).norm() * R_Wigner_Seitz;

    if (std::abs(x) < 1e-6) {
        const double x2 = x * x;
        return 1.0 - 0.1 * x2;  // 1 - x^2/10 + O(x^4)
    }
    return 3.0 * (std::sin(x) - x * std::cos(x)) / (x * x * x);
}

double ElectronPhonon::hole_overlap_integral(int n1, const vector3& k1, int n2, const vector3& k2) const {
    const double cos_angle_k1_k2   = compte_cos_angle(k1, k2);
    const double cos_angle_k1_k2_2 = cos_angle_k1_k2 * cos_angle_k1_k2;
    auto         A_B_params        = m_hole_overlap_int_params.get_params(n1, n2);
    double       integral          = 0.5 * std::sqrt(A_B_params[0] + A_B_params[1] * cos_angle_k1_k2_2);
    return integral;
}

double ElectronPhonon::get_max_phonon_energy() const {
    double max_energy = std::numeric_limits<double>::lowest();
    for (const auto& disp : m_phonon_dispersion) {
        double max_w = disp.max_omega();  // ω_max [1/s]
        if (max_w > max_energy) max_energy = max_w;
    }
    return max_energy * EmpiricalPseudopotential::Constants::h_bar_eV;  // ħω → eV
}

// --- Pairwise kernel: one (n1,k1) → (n2, bary(t)) transition, returns 8 channels ---
// Uses tetra barycenter for k2 and your DOS(Ef) per band at that tetra.
Rate8 ElectronPhonon::compute_transition_rates_pair(int         idx_n1,
                                                    std::size_t idx_k1,
                                                    int         idx_n2,
                                                    std::size_t idx_t /* tetra index */,
                                                    bool        push_nk_npkp /* optional sparse fill */) {
    Rate8 out{};

    const auto& vtx1  = m_list_vertices[idx_k1];
    const auto& tetra = m_list_tetrahedra[idx_t];

    const double  Ei_eV = vtx1.get_energy_at_band(idx_n1);
    const vector3 k1    = vtx1.get_position();
    const vector3 k2    = tetra.compute_barycenter();

    // Overlap once
    const double I  = electron_overlap_integral(k1, k2);
    const double I2 = I * I;

    // q = k2 - k1, fold to 1st BZ if needed (normal processes only)
    vector3 q = k2 - k1;
    if (!is_inside_mesh_geometry(q)) q = retrieve_k_inside_mesh_geometry(q);
    if (!is_inside_mesh_geometry(q)) return out;  // nothing to do

    const double qn = q.norm();

    // Row/col for optional sparse insert (compact conduction indexing assumed)
    Eigen::Index Nk = static_cast<Eigen::Index>(m_list_vertices.size());
    Eigen::Index Nt = static_cast<Eigen::Index>(m_list_tetrahedra.size());

    // If you compact bands for matrices, map full band → compact band index
    auto band_to_comp = [&](int n) -> int {
        auto it = std::find(m_indices_conduction_bands.begin(), m_indices_conduction_bands.end(), n);
        return (it == m_indices_conduction_bands.end()) ? -1 : int(it - m_indices_conduction_bands.begin());
    };
    const int          n1c = band_to_comp(idx_n1);
    const int          n2c = band_to_comp(idx_n2);
    const Eigen::Index row = (n1c >= 0) ? (Eigen::Index)n1c * Nk + (Eigen::Index)idx_k1 : -1;
    const Eigen::Index col = (n2c >= 0) ? (Eigen::Index)n2c * Nt + (Eigen::Index)idx_t : -1;

    constexpr double SMALL_OMEGA_CUTOFF = 1.0;  // [1/s]
    const double     pi                 = EmpiricalPseudopotential::Constants::pi;
    const double     qe                 = EmpiricalPseudopotential::Constants::q_e;       // J/eV
    const double     hbar_eV            = EmpiricalPseudopotential::Constants::h_bar_eV;  // eV·s

    // Loop 4 branches: md=0..3 → (ac/op)×(L/T)
    for (int md = 0; md < 4; ++md) {
        const auto&           disp = m_phonon_dispersion[md];
        const PhononMode      mode = ((md >> 1) == 0) ? PhononMode::acoustic : PhononMode::optical;
        const PhononDirection dir  = ((md & 1) == 0) ? PhononDirection::longitudinal : PhononDirection::transverse;

        // ω(|q|) [1/s] — use your lookup or analytic
        const double omega = disp.omega_lookup(qn);
        if (omega <= SMALL_OMEGA_CUTOFF) continue;

        const double Eph_eV = hbar_eV * omega;
        const double N0     = bose_einstein_distribution(Eph_eV, m_temperature);

        // Deformation potential (J)
        const DeformationPotential& defpot  = (mode == PhononMode::acoustic) ? m_ac_defpot_e : m_op_defpot_e;
        const double                Delta_J = defpot.get_fischetti_deformation_potential(q, idx_n1) * qe;

        // Common prefactor
        const double pref = (pi / (m_rho * omega)) * (Delta_J * Delta_J) * I2 / m_reduce_bz_factor * m_spin_degeneracy;

        // --- Emission (Ef = Ei - ħω), bose = N0 + 1 ---
        {
            const double Ef_eV = Ei_eV - Eph_eV;
            if (tetra.is_energy_inside_band(Ef_eV, idx_n2)) {
                const double dos_eV = tetra.interpolate_dos_at_energy_per_band(Ef_eV, static_cast<std::size_t>(idx_n2));
                if (dos_eV > 0.0) {
                    const double val = pref * (N0 + 1.0) * (dos_eV / qe);
                    const int    b   = rate_index(mode, dir, PhononEvent::emission);
                    out[static_cast<std::size_t>(b)] += val;

                    // optional sparse push
                    // if (push_nk_npkp && row >= 0 && col >= 0) {
                    //     push_nk_npkp(b, row, col, val);
                    // }
                }
            }
        }
        // --- Absorption (Ef = Ei + ħω), bose = N0 ---
        {
            const double Ef_eV = Ei_eV + Eph_eV;
            if (tetra.is_energy_inside_band(Ef_eV, idx_n2)) {
                const double dos_eV = tetra.interpolate_dos_at_energy_per_band(Ef_eV, static_cast<std::size_t>(idx_n2));
                if (dos_eV > 0.0) {
                    const double val = pref * (N0) * (dos_eV / qe);
                    const int    b   = rate_index(mode, dir, PhononEvent::absorption);
                    out[static_cast<std::size_t>(b)] += val;

                    // if (push_nk_npkp && row >= 0 && col >= 0) {
                    //     push_nk_npkp(b, row, col, val);
                    // }
                }
            }
        }
    }

    return out;
}
RateValues ElectronPhonon::compute_electron_phonon_rate(int idx_n1, std::size_t idx_k1, bool populate_nk_npkp) {
    RateValues acc;

    const double Ei_eV   = m_list_vertices[idx_k1].get_energy_at_band(idx_n1);
    const double Eph_max = get_max_phonon_energy();
    const double Ef_min  = Ei_eV - Eph_max;
    const double Ef_max  = Ei_eV + Eph_max;

    std::size_t nnz = 0;

    for (int idx_n2 : m_indices_conduction_bands) {
        // Quick reject band window
        if (Ef_min > m_max_band[idx_n2] || Ef_max < m_min_band[idx_n2]) continue;

        for (std::size_t t = 0; t < m_list_tetrahedra.size(); ++t) {
            const auto& tetra = m_list_tetrahedra[t];
            if (!tetra.does_intersect_band_energy_range(Ef_min, Ef_max, idx_n2)) continue;

            const Rate8 r = compute_transition_rates_pair(idx_n1, idx_k1, idx_n2, t, /*push=*/populate_nk_npkp);
            nnz += (r != Rate8{});  // count non-zero contributions
            for (int i = 0; i < 8; ++i)
                acc.v[i] += r[i];
        }
    }
    std::cout << "Computed rates for (n,k)=(" << idx_n1 << "," << idx_k1 << ") Ei=" << std::setprecision(6) << Ei_eV
              << " eV, non-zero contributions: " << nnz << " / " << m_indices_conduction_bands.size() * m_list_tetrahedra.size() << " = "
              << (100.0 * nnz / (m_indices_conduction_bands.size() * m_list_tetrahedra.size())) << "%\n";
    return acc;
}

// -------------------- Hole rates --------------------

RateValues ElectronPhonon::compute_hole_phonon_rate(int idx_n1, std::size_t idx_k1) {
    RateValues  rates_k1_n1;
    const auto& list_tetrahedra       = m_list_tetrahedra;
    const auto& indices_valence_bands = m_indices_valence_bands;

    const double  Ei_eV = m_list_vertices[idx_k1].get_energy_at_band(idx_n1);
    const vector3 k1    = m_list_vertices[idx_k1].get_position();

    constexpr double SMALL_OMEGA_CUTOFF = 1.0;

    for (int idx_n2 : indices_valence_bands) {
        for (const auto& tetra : list_tetrahedra) {
            const vector3 k2 = tetra.compute_barycenter();

            const double overlap  = hole_overlap_integral(idx_n1, k1, idx_n2, k2);
            const double overlap2 = overlap * overlap;

            vector3 q = k2 - k1;
            if (!is_inside_mesh_geometry(q)) q = retrieve_k_inside_mesh_geometry(q);
            if (!is_inside_mesh_geometry(q)) continue;

            const double q_norm = q.norm();

            for (int md = 0; md < 4; ++md) {
                const auto&           disp = m_phonon_dispersion[md];
                const PhononMode      mode = (md < 2) ? PhononMode::acoustic : PhononMode::optical;
                const PhononDirection dir  = (md & 1) ? PhononDirection::transverse : PhononDirection::longitudinal;

                const double omega = disp.omega_lookup(q_norm);
                if (omega <= SMALL_OMEGA_CUTOFF) continue;

                const double Eph_eV = EmpiricalPseudopotential::Constants::h_bar_eV * omega;
                const double N0     = bose_einstein_distribution(Eph_eV, m_temperature);

                const DeformationPotential& defpot = (mode == PhononMode::acoustic) ? m_ac_defpot_h : m_op_defpot_h;
                const double Delta_J = defpot.get_fischetti_deformation_potential(q, idx_n1) * EmpiricalPseudopotential::Constants::q_e;

                // Emission
                {
                    const double Ef_eV  = Ei_eV - Eph_eV;
                    const double dos_eV = tetra.interpolate_dos_at_energy_per_band(Ef_eV, static_cast<std::size_t>(idx_n2));
                    if (dos_eV > 0.0) {
                        const double dos_per_J = dos_eV / EmpiricalPseudopotential::Constants::q_e;
                        double rate_value = (EmpiricalPseudopotential::Constants::pi / (m_rho * omega)) * (Delta_J * Delta_J) * overlap2 *
                                            (N0 + 1.0) * dos_per_J;
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
                        const double dos_per_J = dos_eV / EmpiricalPseudopotential::Constants::q_e;
                        double       rate_value =
                            (EmpiricalPseudopotential::Constants::pi / (m_rho * omega)) * (Delta_J * Delta_J) * overlap2 * (N0)*dos_per_J;
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

// -------------------- Mesh sweeps & exporters --------------------

void ElectronPhonon::compute_electron_phonon_rates_over_mesh(double energy_max, bool irreducible_wedge_only, bool populate_nk_npkp) {
    auto indices_conduction_bands = m_indices_conduction_bands;
    auto min_idx_conduction_band  = *std::min_element(indices_conduction_bands.begin(), indices_conduction_bands.end());
    auto max_idx_conduction_band  = *std::max_element(indices_conduction_bands.begin(), indices_conduction_bands.end());

    std::cout << "Min index conduction band: " << min_idx_conduction_band << "\n";
    std::cout << "Max index conduction band: " << max_idx_conduction_band << "\n";
    std::cout << "Computing electron-phonon rates over mesh for " << m_list_vertices.size() << " k-points.\n";

    // --- Decide band set for rows/cols (compact conduction example) ---
    const std::vector<int>& rows_bands = m_indices_conduction_bands;  // (n,k)
    const std::vector<int>& cols_bands = m_indices_conduction_bands;  // (n',k')

    if (populate_nk_npkp) {
        const std::size_t Nk      = m_list_vertices.size();
        const std::size_t Nt      = m_list_tetrahedra.size();
        const int         Nb_rows = static_cast<int>(rows_bands.size());
        const int         Nb_cols = static_cast<int>(cols_bands.size());

        const Eigen::Index nrows = static_cast<Eigen::Index>(static_cast<std::size_t>(Nb_rows) * Nk);
        const Eigen::Index ncols = static_cast<Eigen::Index>(static_cast<std::size_t>(Nb_cols) * Nt);
        std::cout << "Matrix size: " << nrows << " x " << ncols << " (" << (nrows * ncols) << " elements)\n";
        constexpr double percentage_reserve = 0.01;  // 1% of elements
        std::cout << "Reserving " << (nrows * ncols * percentage_reserve) << " elements per matrix.\n";

        m_rates_nk_npkp.clear();
        m_rates_nk_npkp.resize(8);  // indices 0..7 must equal rate_index(m,d,e)

        Eigen::VectorXi reserve_vec(nrows);
        reserve_vec.setConstant(std::max<Eigen::Index>(1, static_cast<Eigen::Index>(percentage_reserve * Nb_cols * Nt)));

        for (int M = 0; M < 2; ++M) {
            const PhononMode m = (M == 0) ? PhononMode::acoustic : PhononMode::optical;
            for (int D = 0; D < 2; ++D) {
                const PhononDirection d = (D == 0) ? PhononDirection::longitudinal : PhononDirection::transverse;
                for (int E = 0; E < 2; ++E) {
                    const PhononEvent e   = (E == 0) ? PhononEvent::absorption : PhononEvent::emission;  // matches your enum
                    const int         idx = rate_index(m, d, e);                                         // (M<<2)|(D<<1)|E  → [0..7]

                    Rates_nk_npkp_ctor R;
                    R.mode      = m;
                    R.direction = d;
                    R.event     = e;
                    R.matrix.resize(nrows, ncols);
                    R.matrix.reserve(reserve_vec);

                    m_rates_nk_npkp[idx] = std::move(R);
                }
            }
        }
    } // if populate_nk_npkp

    std::cout << "Progress: 0%";
    std::atomic<std::size_t> counter{0};
    constexpr int            chunk_size = 32;

    // Create a shuffled list of indices to balance load when using irreducible wedge only
    std::vector<std::size_t> random_indices(m_list_vertices.size());
    for (std::size_t i = 0; i < m_list_vertices.size(); ++i)
        random_indices[i] = i;
    if (irreducible_wedge_only) {
        constexpr int seed = 0;
        std::mt19937  g(seed);
        std::shuffle(random_indices.begin(), random_indices.end(), g);
    }

#pragma omp parallel for schedule(dynamic, chunk_size)
    for (std::size_t idx = 0; idx < m_list_vertices.size(); ++idx) {
        // Use random index if irreducible wedge only
        const std::size_t idx_k1 = irreducible_wedge_only ? random_indices[idx] : idx;
        const bool        to_compute =
            (!irreducible_wedge_only) || (irreducible_wedge_only && is_irreducible_wedge(m_list_vertices[idx_k1].get_position()));

        const auto done = ++counter;
        if (omp_get_thread_num() == 0) {
            std::cout << "\rDone " << done << "/" << m_list_vertices.size() << " (" << std::fixed << std::setprecision(1)
                      << (100.0 * done / m_list_vertices.size()) << "%)" << std::flush;
        }

        // Valence bands (holes)
        for (std::size_t idx_n1 = 0; idx_n1 < static_cast<std::size_t>(min_idx_conduction_band); ++idx_n1) {
            if (!to_compute) continue;
            auto hole_rate = compute_hole_phonon_rate(static_cast<int>(idx_n1), idx_k1);
            (void)hole_rate;  // attach to vertices if/when needed
            // m_list_vertices[idx_k1].add_electron_phonon_rates(hole_rate.as_array());
        }

        // Conduction bands (electrons)
        for (std::size_t idx_n1 = static_cast<std::size_t>(min_idx_conduction_band);
             idx_n1 <= static_cast<std::size_t>(max_idx_conduction_band);
             ++idx_n1) {
            if (!to_compute) continue;

            if (m_list_vertices[idx_k1].get_energy_at_band(static_cast<int>(idx_n1)) > energy_max) {
                m_list_vertices[idx_k1].add_electron_phonon_rates(std::array<double, 8>{});
                continue;
            } else {
                auto rate = compute_electron_phonon_rate(static_cast<int>(idx_n1), idx_k1);
                m_list_vertices[idx_k1].add_electron_phonon_rates(rate.as_array());
            }
        }
    }

    std::cout << "\rComputed " << counter << " k-points out of " << m_list_vertices.size() << " (100%)\n";

    if (irreducible_wedge_only) {
        std::cout << "Set electron-phonon rates for all mesh vertices.\n";
#pragma omp parallel for schedule(dynamic)
        for (std::size_t idx_k1 = 0; idx_k1 < m_list_vertices.size(); ++idx_k1) {
            if ((idx_k1 % 1000) == 0 && omp_get_thread_num() == 0) {
                std::cout << "\rSetting rates for all k-points: " << idx_k1 << "/" << m_list_vertices.size() << " (" << std::fixed
                          << std::setprecision(1) << (100.0 * idx_k1 / m_list_vertices.size()) << "%)" << std::flush;
            }
            if (!is_irreducible_wedge(m_list_vertices[idx_k1].get_position())) {
                std::size_t idx_k1_symm = get_index_irreducible_wedge(m_list_vertices[idx_k1].get_position());
                for (std::size_t idx_n1 = 0; idx_n1 <= static_cast<std::size_t>(max_idx_conduction_band); ++idx_n1) {
                    auto rates_symm = m_list_vertices[idx_k1_symm].get_electron_phonon_rates(idx_n1);
                    m_list_vertices[idx_k1].add_electron_phonon_rates(rates_symm);
                }
            }
        }
        std::cout << "\rSet rates for all k-points: " << m_list_vertices.size() << "/" << m_list_vertices.size() << " (100%)\n";
    }
}

void ElectronPhonon::compute_electron_phonon_rates_over_mesh_nk_npkp(bool irreducible_wedge_only) {
    auto indices_conduction_bands = m_indices_conduction_bands;
    auto min_idx_conduction_band  = *std::min_element(indices_conduction_bands.begin(), indices_conduction_bands.end());

    std::cout << "Min index conduction band: " << min_idx_conduction_band << "\n";
    std::cout << "Computing electron-phonon rates over mesh for " << m_list_vertices.size() << " k-points.\n";

    std::cout << "Progress: 0%";
    std::atomic<std::size_t> counter{0};

#pragma omp parallel for schedule(dynamic)
    for (std::size_t idx_k1 = 0; idx_k1 < m_list_vertices.size(); ++idx_k1) {
        bool to_compute = is_irreducible_wedge(m_list_vertices[idx_k1].get_position()) && irreducible_wedge_only;

        auto done = ++counter;
        if ((done % 100) == 0 || (done == m_list_vertices.size() && omp_get_thread_num() == 0)) {
            std::cout << "\rDone " << done << "/" << m_list_vertices.size() << " (" << std::fixed << std::setprecision(1)
                      << (100.0 * done / m_list_vertices.size()) << "%)" << std::flush;
        }

        for (std::size_t idx_n1 = 0; idx_n1 < static_cast<std::size_t>(min_idx_conduction_band); ++idx_n1) {
            if (!to_compute) continue;
            (void)compute_hole_phonon_rate(static_cast<int>(idx_n1), idx_k1);
        }
    }
}

// -------------------- Sampling a final state --------------------

std::pair<int, std::size_t> ElectronPhonon::select_final_state(std::size_t     idx_band_initial,
                                                               std::size_t     idx_k_initial,
                                                               PhononMode      mode,
                                                               PhononDirection direction,
                                                               PhononEvent     event) const {
    using std::size_t;

    if (m_list_vertices.empty() || m_list_tetrahedra.empty()) throw std::runtime_error("select_final_state: empty mesh.");
    if (idx_k_initial >= m_list_vertices.size()) throw std::out_of_range("select_final_state: idx_k_initial OOB.");

    const int md = md_index(mode, direction);
    if (md < 0) throw std::runtime_error("select_final_state: invalid mode/direction.");
    const auto& disp = m_phonon_dispersion[md];  // ω(|q|) in s^-1

    // Initial state
    const double  Ei_eV = m_list_vertices[idx_k_initial].get_energy_at_band(static_cast<int>(idx_band_initial));
    const vector3 k1    = m_list_vertices[idx_k_initial].get_position();

    const double sign_ph = (event == PhononEvent::emission) ? -1.0 : +1.0;

    const size_t        nb_bands = m_list_vertices.front().get_number_bands();
    const size_t        nb_tetra = m_list_tetrahedra.size();
    std::vector<double> probs_flat(nb_bands * nb_tetra, 0.0);
    auto                P_ref = [&](int n2, size_t t) -> double& { return probs_flat[static_cast<size_t>(n2) * nb_tetra + t]; };

    const double pi      = EmpiricalPseudopotential::Constants::pi;
    const double qe      = EmpiricalPseudopotential::Constants::q_e;
    const double hbar_eV = EmpiricalPseudopotential::Constants::h_bar_eV;

    const double Eph_max_eV = get_max_phonon_energy();

    for (int n2 : m_indices_conduction_bands) {
        // Band window
        const double Ef_min = Ei_eV - Eph_max_eV;
        const double Ef_max = Ei_eV + Eph_max_eV;
        if (Ef_min > m_max_band[n2] || Ef_max < m_min_band[n2]) continue;

        for (size_t t = 0; t < nb_tetra; ++t) {
            const auto& tetra = m_list_tetrahedra[t];
            if (!tetra.does_intersect_band_energy_range(Ef_min, Ef_max, n2)) continue;

            const vector3 k2_bary = tetra.compute_barycenter();

            vector3 q = k2_bary - k1;
            if (!is_inside_mesh_geometry(q)) q = retrieve_k_inside_mesh_geometry(q);
            if (!is_inside_mesh_geometry(q)) continue;

            const double qn    = q.norm();
            const double omega = disp.omega_lookup(qn);
            if (!(omega > 0.0) || omega < 1e-12) continue;

            const double Eph_eV = hbar_eV * omega;
            const double N0     = bose_einstein_distribution(Eph_eV, m_temperature);
            const double bose   = (sign_ph < 0.0) ? (N0 + 1.0) : N0;

            const double Ef_eV  = Ei_eV + sign_ph * Eph_eV;
            const double dos_eV = tetra.interpolate_dos_at_energy_per_band(Ef_eV, static_cast<std::size_t>(n2));
            if (dos_eV <= 0.0) continue;

            const double dos_per_J = dos_eV / qe;

            const double I  = electron_overlap_integral(k1, k2_bary);
            const double I2 = I * I;

            const DeformationPotential& defpot  = (mode == PhononMode::acoustic) ? m_ac_defpot_e : m_op_defpot_e;
            const double                Delta_J = defpot.get_fischetti_deformation_potential(q, static_cast<int>(idx_band_initial)) * qe;

            double P = (pi / (m_rho * omega)) * (Delta_J * Delta_J) * I2 * bose * dos_per_J;
            P /= m_reduce_bz_factor;
            P *= m_spin_degeneracy;

            if (P > 0.0 && std::isfinite(P)) P_ref(n2, t) = P;
        }
    }

    double total = 0.0;
    for (double p : probs_flat)
        total += p;
    if (!(total > 0.0) || !std::isfinite(total))
        throw std::runtime_error("select_final_state: no admissible final states (total probability = 0).");

    // Thread-local RNG
    thread_local std::mt19937_64           rng([] {
        std::random_device rd;
        auto               s1 = static_cast<uint64_t>(rd());
        auto               s2 = static_cast<uint64_t>(rd());
        return (s1 << 32) ^ s2;
    }());
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
    // Fallback: last positive
    for (size_t flat = probs_flat.size(); flat-- > 0;) {
        if (probs_flat[flat] > 0.0) {
            const int    n2 = static_cast<int>(flat / nb_tetra);
            const size_t t  = static_cast<size_t>(flat % nb_tetra);
            return {n2, t};
        }
    }
    throw std::runtime_error("select_final_state: internal sampling error.");
}

// -------------------- Energy sweep exporter --------------------

void ElectronPhonon::compute_plot_electron_phonon_rates_vs_energy_over_mesh(int                nb_bands,
                                                                            double             max_energy,
                                                                            double             energy_step,
                                                                            const std::string& filename,
                                                                            bool /*irreducible_wedge_only*/) {
    if (energy_step <= 0.0) throw std::invalid_argument("energy_step must be > 0");
    if (max_energy < 0.0) throw std::invalid_argument("max_energy must be >= 0");

    if (m_list_vertices.empty()) throw std::runtime_error("No vertices in mesh.");
    const int total_bands = static_cast<int>(m_list_vertices.front().get_number_bands());
    if (total_bands <= 0) throw std::runtime_error("Mesh reports zero bands.");

    nb_bands = std::clamp(nb_bands, 0, total_bands);
    if (nb_bands == 0) throw std::runtime_error("nb_bands clamped to 0; nothing to process.");

    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open " + filename + " for writing.");

    static constexpr std::array<const char*, 8>
        kLabels{"ac_L_ab", "ac_T_ab", "op_L_ab", "op_T_ab", "ac_L_em", "ac_T_em", "op_L_em", "op_T_em"};

    out << "# E[eV],DOS(E)";
    for (auto* s : kLabels)
        out << ',' << s;
    out << '\n';

    const std::size_t n_steps = static_cast<std::size_t>(std::floor(max_energy / energy_step)) + 1;

    for (std::size_t istep = 0; istep < n_steps; ++istep) {
        const double E = std::min(max_energy, istep * energy_step);
        std::cout << "\rEnergy: " << E << " / " << max_energy << std::flush;

        double                dos_sum = 0.0;
        std::array<double, 8> num{};  // accumulators
        num.fill(0.0);

        for (const auto& tetra : m_list_tetrahedra) {
            for (int b = 0; b < nb_bands; ++b) {
                const double dos_t = tetra.compute_tetra_dos_energy_band(E, static_cast<std::size_t>(b));
                if (!std::isfinite(dos_t))
                    throw std::runtime_error("DOS is NaN/Inf at E=" + std::to_string(E) + " band=" + std::to_string(b));
                if (dos_t <= 0.0) continue;

                dos_sum += dos_t;

                const std::array<double, 8> rates = tetra.get_tetra_electron_phonon_rates(static_cast<std::size_t>(b));
                for (int i = 0; i < 8; ++i) {
                    num[i] += rates[i] * dos_t;
                }
            }
        }

        std::array<double, 8> mean{};
        if (dos_sum > 0.0) {
            const double inv_dos = 1.0 / dos_sum;
            for (int i = 0; i < 8; ++i)
                mean[i] = num[i] * inv_dos;
        } else {
            mean.fill(0.0);
        }

        out << std::scientific << std::setprecision(10) << E << ',' << dos_sum;
        for (double v : mean)
            out << ',' << v;
        out << '\n';
    }

    std::cout << std::endl;
    out.close();
}

/**
 * @brief Interpolate the phonon scattering rate at a given location for a given band.
 *
 * @param location
 * @param idx_band
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
    Rate8            rates;
    constexpr double inv_num_vertices = 1.0 / 4.0;
    for (std::size_t i = 0; i < vertex_indices.size(); ++i) {
        const auto& vertex = m_list_vertices[vertex_indices[i]];
        for (std::size_t idx_mode = 0; idx_mode < rates.size(); ++idx_mode) {
            rates[idx_mode] += m_list_phonon_scattering_rates[vertex.get_index()][idx_band][idx_mode];
        }
    }
    for (std::size_t idx_mode = 0; idx_mode < rates.size(); ++idx_mode) {
        rates[idx_mode] *= inv_num_vertices;
    }
    return rates;
}

// -------------------- Phonon dispersion dump --------------------

void ElectronPhonon::plot_phonon_dispersion(const std::string& filename) const {
    std::ofstream file(filename);
    for (auto&& vtx : m_list_vertices) {
        auto k = vtx.get_position();
        file << k.x() << " " << k.y() << " " << k.z() << " ";
        for (int md = 0; md < 4; ++md) {
            auto q = k;
            q /= m_material.get_fourier_factor();
            const auto&  disp = m_phonon_dispersion[md];
            const double e_ph = disp.omega_lookup(q.norm()) * EmpiricalPseudopotential::Constants::h_bar_eV;
            file << e_ph << " ";
        }
        file << '\n';
    }
}

// -------------------- Gmsh export --------------------

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

    auto indices_conduction_bands = m_indices_conduction_bands;
    auto max_idx_conduction_band  = *std::max_element(indices_conduction_bands.begin(), indices_conduction_bands.end());
    int  nb_bands                 = static_cast<int>(m_indices_conduction_bands.size() + m_indices_valence_bands.size());

    for (int idx_val_band = 0; idx_val_band < max_idx_conduction_band; ++idx_val_band) {
        std::vector<double> rates_ac_lo_em(m_list_vertices.size());
        std::vector<double> rates_ac_lo_ab(m_list_vertices.size());
        std::vector<double> rates_ac_tr_em(m_list_vertices.size());
        std::vector<double> rates_ac_tr_ab(m_list_vertices.size());
        std::vector<double> rates_opt_lo_em(m_list_vertices.size());
        std::vector<double> rates_opt_lo_ab(m_list_vertices.size());
        std::vector<double> rates_opt_tr_em(m_list_vertices.size());
        std::vector<double> rates_opt_tr_ab(m_list_vertices.size());

        for (std::size_t idx_k1 = 0; idx_k1 < m_list_vertices.size(); ++idx_k1) {
            auto rates              = m_list_vertices[idx_k1].get_electron_phonon_rates(static_cast<std::size_t>(idx_val_band));
            rates_ac_lo_em[idx_k1]  = rates[0];
            rates_ac_lo_ab[idx_k1]  = rates[1];
            rates_ac_tr_em[idx_k1]  = rates[2];
            rates_ac_tr_ab[idx_k1]  = rates[3];
            rates_opt_lo_em[idx_k1] = rates[4];
            rates_opt_lo_ab[idx_k1] = rates[5];
            rates_opt_tr_em[idx_k1] = rates[6];
            rates_opt_tr_ab[idx_k1] = rates[7];
        }

        std::string name_rate_ac_lo_em  = "ac_lo_em_" + std::to_string(idx_val_band);
        std::string name_rate_ac_lo_ab  = "ac_lo_ab_" + std::to_string(idx_val_band);
        std::string name_rate_ac_tr_em  = "ac_tr_em_" + std::to_string(idx_val_band);
        std::string name_rate_ac_tr_ab  = "ac_tr_ab_" + std::to_string(idx_val_band);
        std::string name_rate_opt_lo_em = "opt_lo_em_" + std::to_string(idx_val_band);
        std::string name_rate_opt_lo_ab = "opt_lo_ab_" + std::to_string(idx_val_band);
        std::string name_rate_opt_tr_em = "opt_tr_em_" + std::to_string(idx_val_band);
        std::string name_rate_opt_tr_ab = "opt_tr_ab_" + std::to_string(idx_val_band);

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

// -------------------- Load parameters (YAML) --------------------

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

    double radiusWS      = material["Radius-WS"].as<double>();
    m_radii_wigner_seitz = radiusWS;

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
            if (md < 0) throw std::runtime_error("load_phonon_parameters: bad mode/direction index.");
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
                if (mode == PhononMode::acoustic)
                    m_ac_defpot_e = deformationPotential;
                else
                    m_op_defpot_e = deformationPotential;
            } else {
                if (mode == PhononMode::acoustic)
                    m_ac_defpot_h = deformationPotential;
                else
                    m_op_defpot_h = deformationPotential;
            }
        }
    }
}

// -------------------- CSV import/export of rates --------------------

void ElectronPhonon::export_rate_values(const std::string& filename) const {
    std::ofstream file(filename);
    for (auto&& vertex : m_list_vertices) {
        std::vector<std::array<double, 8>> all_rates = vertex.get_electron_phonon_rates_all_bands();
        for (std::size_t idx_band = 0; idx_band < all_rates.size(); ++idx_band) {
            double energy = vertex.get_energy_at_band(idx_band);
            file << idx_band << "," << energy << ",";
            for (std::size_t idx_rate = 0; idx_rate < all_rates[idx_band].size(); ++idx_rate) {
                double rate = all_rates[idx_band][idx_rate];
                file << rate << ((idx_rate + 1 < all_rates[idx_band].size()) ? "," : "");
            }
            file << '\n';
        }
    }
    file.close();
}

void ElectronPhonon::read_phonon_scattering_rates_from_file(const std::filesystem::path& path) {
    std::cout << "Reading phonon scattering rates (CSV) from file " << path.string() << " ...\n";

    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("Could not open file " + path.string());

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

            for (char& c : line)
                if (c == ',') c = ' ';
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

// -------------------- Simple reductions --------------------

inline double ElectronPhonon::sum_modes(const Rate8& r) const noexcept {
    double s = 0.0;
    for (int i = 0; i < 8; ++i)
        s += r[i];
    return s;
}

double ElectronPhonon::compute_P_Gamma() const {
    double pgamma_max = 0.0;
    for (const auto& perVertex : m_list_phonon_scattering_rates) {
        for (const auto& rate8 : perVertex) {
            const double tot = sum_modes(rate8);
            if (tot > pgamma_max) pgamma_max = tot;
        }
    }
    return pgamma_max;
}

}  // namespace bz_mesh
