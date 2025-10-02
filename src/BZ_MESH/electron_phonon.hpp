/**
 * @file elelectron_phonon.hpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-02-09
 *
 * @copyright Copyright (c) 2024
 *
 * Phonon are always stored in the same order:
 * 0: acoustic longitudinal absorption
 * 1: acoustic transverse absorption
 * 2: optical longitudinal absorption
 * 3: optical transverse absorption
 *
 */

#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

#include "bz_states.hpp"
#include "yaml-cpp/yaml.h"

// Eigen space

namespace bz_mesh {

enum class PhononMode { acoustic, optical, none };
enum class PhononDirection { longitudinal, transverse, none };
enum class PhononEvent { absorption, emission, none };

// Order: ac/opt × L/T × ab/em  → indices 0..7
constexpr int rate_index(PhononMode m, PhononDirection d, PhononEvent e) {
    int M = (m == PhononMode::acoustic ? 0 : m == PhononMode::optical ? 1 : -1);
    int D = (d == PhononDirection::longitudinal ? 0 : d == PhononDirection::transverse ? 1 : -1);
    int E = (e == PhononEvent::absorption ? 0 : e == PhononEvent::emission ? 1 : -1);
    return (M < 0 || D < 0 || E < 0) ? -1 : ((M * 2 + D) * 2 + E);
}

// Named indices (handy for headers/export)
constexpr int IDX_AC_L_AB = rate_index(PhononMode::acoustic, PhononDirection::longitudinal, PhononEvent::absorption);  // 0
constexpr int IDX_AC_T_AB = rate_index(PhononMode::acoustic, PhononDirection::transverse, PhononEvent::absorption);    // 1
constexpr int IDX_OP_L_AB = rate_index(PhononMode::optical, PhononDirection::longitudinal, PhononEvent::absorption);   // 2
constexpr int IDX_OP_T_AB = rate_index(PhononMode::optical, PhononDirection::transverse, PhononEvent::absorption);     // 3
constexpr int IDX_AC_L_EM = rate_index(PhononMode::acoustic, PhononDirection::longitudinal, PhononEvent::emission);    // 4
constexpr int IDX_AC_T_EM = rate_index(PhononMode::acoustic, PhononDirection::transverse, PhononEvent::emission);      // 5
constexpr int IDX_OP_L_EM = rate_index(PhononMode::optical, PhononDirection::longitudinal, PhononEvent::emission);     // 6
constexpr int IDX_OP_T_EM = rate_index(PhononMode::optical, PhononDirection::transverse, PhononEvent::emission);       // 7

static_assert(IDX_AC_L_AB == 0 && IDX_OP_T_EM == 7, "Rate index mapping changed; update I/O accordingly.");

struct PhononDispersion {
    PhononMode      m_mode      = PhononMode::none;
    PhononDirection m_direction = PhononDirection::none;

    // Analytic model: ω(q) = m_omega + m_vs * q + m_c * q^2
    // q in 1/m, ω(q) in rad/s (keep exactly as used elsewhere in your code)
    double m_omega = 0.0;
    double m_vs    = 0.0;
    double m_c     = 0.0;

    // Lookup tables (uniform grid on q); names kept as-is
    std::vector<double> m_q_values;  // in 1/m
    std::vector<double> m_e_values;  // stores ω(q) samples (rad/s)

    PhononDispersion() = default;

    PhononDispersion(PhononMode mode, PhononDirection direction, double omega, double vs, double c)
        : m_mode(mode),
          m_direction(direction),
          m_omega(omega),
          m_vs(vs),
          m_c(c) {}

    PhononDispersion(PhononMode mode, PhononDirection direction)
        : m_mode(mode),
          m_direction(direction),
          m_omega(0.0),
          m_vs(0.0),
          m_c(0.0) {}

    // Analytic dispersion: ω(q) in rad/s
    double get_phonon_dispersion(double q) const {
        const double q_squared   = q * q;
        const double linear_part = std::fma(m_vs, q, m_omega);             // m_omega + m_vs*q
        const double omega_q     = std::fma(m_c, q_squared, linear_part);  // + m_c*q^2

        return omega_q;
    }

    // Build uniform lookup on [0, q_max]; stores ω(q) in m_e_values
    void fill_lookup_table(double q_max, std::size_t n_points) {
        if (q_max <= 0.0) {
            throw std::invalid_argument("q_max must be > 0");
        }
        if (n_points < 2) {
            throw std::invalid_argument("n_points must be >= 2");
        }

        m_q_values.resize(n_points);
        m_e_values.resize(n_points);

        const std::size_t grid_point_count = n_points;
        const double      grid_spacing     = q_max / static_cast<double>(grid_point_count - 1);

        for (std::size_t sample_index = 0; sample_index < grid_point_count; ++sample_index) {
            const double q_value     = static_cast<double>(sample_index) * grid_spacing;
            m_q_values[sample_index] = q_value;
            m_e_values[sample_index] = get_phonon_dispersion(q_value);  // ω(q) in rad/s
        }
    }

    // Fast O(1) lookup with linear interpolation; returns ω(q) in rad/s
    double get_phonon_dispersion_from_lookup(double q) const {
        if (m_q_values.empty() || m_e_values.empty()) {
            throw std::runtime_error("Lookup table is empty");
        }

        const std::size_t grid_point_count = m_q_values.size();
        const double      q_minimum        = m_q_values.front();
        const double      q_maximum        = m_q_values.back();

        if (q < q_minimum || q > q_maximum) {
            std::cerr << "Q value: " << q << " out of bounds [" << q_minimum << ", " << q_maximum << "]" << std::endl;
            throw std::runtime_error("Q value out of bounds of lookup table");
        }

        // Uniform grid → compute bin directly (no search)
        const double grid_spacing      = (q_maximum - q_minimum) / static_cast<double>(grid_point_count - 1);
        const double inverse_spacing   = 1.0 / grid_spacing;
        const double position_in_grid  = (q - q_minimum) * inverse_spacing;  // in “grid steps”
        std::size_t  left_sample_index = static_cast<std::size_t>(position_in_grid);
        if (left_sample_index >= grid_point_count - 1) {
            left_sample_index = grid_point_count - 2;  // guard top edge
        }

        const double fraction_right = position_in_grid - static_cast<double>(left_sample_index);
        const double left_value     = m_e_values[left_sample_index];
        const double right_value    = m_e_values[left_sample_index + 1];

        // Linear interpolation with FMA: left + t * (right - left)
        return std::fma(fraction_right, (right_value - left_value), left_value);
    }

    double get_max_energy_exchange() const {
        if (m_e_values.empty()) {
            throw std::runtime_error("Lookup table is empty");
        }
        return *std::max_element(m_e_values.begin(), m_e_values.end());
    }
};

struct DeformationPotential {
    PhononMode m_mode             = PhononMode::none;
    double     m_A                = 0.0;
    double     m_B                = 0.0;
    double     m_energy_threshold = 1e6;  // eV, default no threshold

    DeformationPotential() = default;
    DeformationPotential(PhononMode type, double A, double B, double energy_threshold)
        : m_mode(type),
          m_A(A),
          m_B(B),
          m_energy_threshold(energy_threshold) {}

    double get_deformation_potential(const vector3& q, double energy) const {
        double clamp_energy = std::min(energy, m_energy_threshold);
        if (m_mode == PhononMode::acoustic) {
            return std::sqrt(m_A + clamp_energy * m_B) * q.norm();
        } else {
            return std::sqrt(m_A + clamp_energy * m_B);
        }
    }

    double get_fischetti_deformation_potential(const vector3& q, int idx_band) const {
        constexpr double cm_to_m = 1e2;

        double boost_factor = 1.5;  // empirical boost to match bulk mobility

        if (m_mode == PhononMode::acoustic) {
            if (idx_band == 0) {
                return boost_factor * 1.2 * q.norm();  // in eV/m
            } else {
                return 1.0 * 1.7 * q.norm();  // in eV/m
            }
        } else {
            if (idx_band == 0) {
                return boost_factor * 1.75e8 * cm_to_m;  // in eV/m
            } else {
                return boost_factor * 2.10e8 * cm_to_m;  // in eV/m
            }
        }
    }

    void set_energy_threshold(double energy) { m_energy_threshold = energy; }

    void print() const {
        std::cout << "Deformation potential for mode "
                  << (m_mode == PhononMode::acoustic  ? "acoustic"
                      : m_mode == PhononMode::optical ? "optical"
                                                      : "none")
                  << ": A = " << m_A << ", B = " << m_B << ", energy threshold = " << m_energy_threshold << " eV" << std::endl;
    }
};

struct HoleOverlapIntParams {
    std::vector<std::pair<std::pair<int, int>, double>> A_params =
        {{{1, 1}, 1.0}, {{2, 2}, 1.0}, {{3, 3}, 5.0 / 8.0}, {{1, 2}, 3.0}, {{1, 3}, 3.0 / 8.0}, {{2, 3}, 3.0 / 8.0}};

    std::vector<std::pair<std::pair<int, int>, double>> B_params =
        {{{1, 1}, 3.0}, {{2, 2}, 3.0}, {{3, 3}, 0.0}, {{1, 2}, -3.0}, {{1, 3}, 0.0}, {{2, 3}, 0.0}};

    /**
     * @brief Return the parameters for the overlap integral, given the indices of the bands.
     *
     * @param n1
     * @param n2
     * @return std::array<double, 2>
     */
    std::array<double, 2> get_params(int n1, int n2) const {
        for (const auto& p : A_params) {
            if (p.first.first == n1 && p.first.second == n2) {
                return {p.second, 0.0};
            }
        }
        for (const auto& p : B_params) {
            if (p.first.first == n1 && p.first.second == n2) {
                return {0.0, p.second};
            }
        }
        return {0.0, 0.0};
    }

    HoleOverlapIntParams() = default;
};

struct RateValue {
    PhononMode      m_mode      = PhononMode::none;
    PhononDirection m_direction = PhononDirection::none;
    PhononEvent     m_type      = PhononEvent::none;
    double          m_value     = 0.0;

    RateValue(PhononMode mode, PhononDirection direction, PhononEvent type, double value)
        : m_mode(mode),
          m_direction(direction),
          m_type(type),
          m_value(value) {}
};

struct RateValues {
    std::array<double, 8> values{};  // zero-initialized

    void print_rates() const {
        std::cout << "Acoustic longitudinal absorption : " << values[IDX_AC_L_AB] << "\n"
                  << "Acoustic transverse absorption   : " << values[IDX_AC_T_AB] << "\n"
                  << "Optical longitudinal absorption  : " << values[IDX_OP_L_AB] << "\n"
                  << "Optical transverse absorption    : " << values[IDX_OP_T_AB] << "\n"
                  << "Acoustic longitudinal emission   : " << values[IDX_AC_L_EM] << "\n"
                  << "Acoustic transverse emission     : " << values[IDX_AC_T_EM] << "\n"
                  << "Optical longitudinal emission    : " << values[IDX_OP_L_EM] << "\n"
                  << "Optical transverse emission      : " << values[IDX_OP_T_EM] << std::endl;
    }

    const std::array<double, 8>& to_array() const { return values; }

    void add_rate(const RateValue& rate) {
        if (int idx = rate_index(rate.m_mode, rate.m_direction, rate.m_type); idx >= 0) values[idx] += rate.m_value;
    }

    // Optional: direct accessors if you like named use sites
    double&       at(PhononMode m, PhononDirection d, PhononEvent e) { return values[rate_index(m, d, e)]; }
    const double& at(PhononMode m, PhononDirection d, PhononEvent e) const { return values[rate_index(m, d, e)]; }

    inline static const char* rate_label(int i) {
        static const char* L[8] = {"ac_L_ab", "ac_T_ab", "op_L_ab", "op_T_ab", "ac_L_em", "ac_T_em", "op_L_em", "op_T_em"};
        return (i >= 0 && i < 8) ? L[i] : "invalid";
    }
};

typedef std::pair<PhononMode, PhononDirection> PhononModeDirection;
typedef Eigen::SparseMatrix<double>            EigenSparseMatrix;
typedef Eigen::Triplet<double>                 EigenTriplet;

class ElectronPhonon : public BZ_States {
 private:
    double               m_temperature        = 300.0;
    double               m_rho                = 2.329e3;
    double               m_radii_wigner_seitz = 0.0;
    HoleOverlapIntParams m_hole_overlap_int_params;
    DeformationPotential m_acoustic_deformation_potential_e;
    DeformationPotential m_optical_deformation_potential_e;
    DeformationPotential m_acoustic_deformation_potential_h;
    DeformationPotential m_optical_deformation_potential_h;

    std::map<PhononModeDirection, PhononDispersion> m_phonon_dispersion;

    /**
     * @brief Count the number of tetrahedra connected to each vertex in the el-ph computation.
     * m_count_weight_tetra_per_vertex[vertex].
     */
    std::vector<double> m_count_weight_tetra_per_vertex;

    /**
     * @brief Phonon rates.
     * m_phonon_rates[mode] gives the phonon dispersion for the given mode.
     * Each matrix stores all the ((n,k) → (n',k')) rates for electron-phonon scattering.
     * It is a sparse matrix, as most transitions are not allowed (energy/momentum conservation
     * rules).
     * M_i,j is the rate from state i to state j.
     * States are indexed as n * N_k + k, where n is the band index and k the k-point index.
     *
     *
     */
    std::vector<EigenSparseMatrix> m_phonon_nk_npkp_modes;

 public:
    explicit ElectronPhonon(const EmpiricalPseudopotential::Material& material) : BZ_States(material) {}

    void   load_phonon_parameters(const std::string& filename);
    void   plot_phonon_dispersion(const std::string& filename) const;
    double get_max_phonon_energy() const;

    inline double bose_einstein_distribution(double energy, double temperature);
    double        electron_overlap_integral(const vector3& k1, const vector3& k2);
    double        hole_overlap_integral(int n1, const vector3& k1, int n2, const vector3& k2);

    RateValues compute_electron_phonon_rate(int idx_n1, std::size_t idx_k1, bool populate_nk_npkp = false);
    RateValues compute_hole_phonon_rate(int idx_n1, std::size_t idx_k1);

    void set_temperature(double temperature) { m_temperature = temperature; }
    void set_density(double rho) { m_rho = rho; }

    void compute_electron_phonon_rates_over_mesh(bool irreducible_wedge_only = false);
    void add_electron_phonon_rates_to_mesh(const std::string& initial_filename, const std::string& final_filename);
    void compute_electron_phonon_rates_over_mesh_nk_npkp(bool irreducible_wedge_only = false);

    std::pair<int, std::size_t> select_final_state(std::size_t idx_band_initial,
                                                   std::size_t idx_k_initial,
                                                   PhononMode  mode,
                                                   PhononDirection direction,
                                                   PhononEvent event) const;

    void export_rate_values(const std::string& filename) const;

    void compute_plot_electron_phonon_rates_vs_energy_over_mesh(int                nb_bands,
                                                                double             max_energy,
                                                                double             energy_step,
                                                                const std::string& filename,
                                                                bool               irreducible_wedge_only = false);
};

}  // namespace bz_mesh