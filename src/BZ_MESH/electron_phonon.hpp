/**
 * @file electron_phonon.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2025-10-04
 *
 *
 */

#pragma once

#include <Eigen/Sparse>  // needed for member type in class; keep here
#include <algorithm>     // max_element
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "bz_states.hpp"  // for vector3, BZ_States, material

namespace bz_mesh {

// ---------- Discrete axes ----------
enum class PhononMode : uint8_t { acoustic = 0, optical = 1, none = 2 };
enum class PhononDirection : uint8_t { longitudinal = 0, transverse = 1, none = 2 };
enum class PhononEvent : uint8_t { absorption = 0, emission = 1, none = 2 };

// Compact index:  m∈{ac,op}, d∈{L,T}, e∈{ab,em}  → [0..7] = (m<<2)|(d<<1)|e
constexpr int rate_index(PhononMode m, PhononDirection d, PhononEvent e) noexcept {
    const int M = (m == PhononMode::acoustic) ? 0 : (m == PhononMode::optical) ? 1 : -1;
    const int D = (d == PhononDirection::longitudinal) ? 0 : (d == PhononDirection::transverse) ? 1 : -1;
    const int E = (e == PhononEvent::absorption) ? 0 : (e == PhononEvent::emission) ? 1 : -1;
    return (M | D | E) < 0 ? -1 : ((M << 2) | (D << 1) | E);
}

constexpr int IDX_AC_L_AB = 0;
constexpr int IDX_AC_T_AB = 1;
constexpr int IDX_OP_L_AB = 2;
constexpr int IDX_OP_T_AB = 3;
constexpr int IDX_AC_L_EM = 4;
constexpr int IDX_AC_T_EM = 5;
constexpr int IDX_OP_L_EM = 6;
constexpr int IDX_OP_T_EM = 7;

// ---------- Compact rate containers ----------
using Rate8           = std::array<double, 8>;
using BandRates       = std::vector<Rate8>;      // per band
using VertexBandRates = std::vector<BandRates>;  // per vertex

struct RateValues {
    Rate8 v{};
    void  add(PhononMode m, PhononDirection d, PhononEvent e, double x) noexcept {
        const int idx = rate_index(m, d, e);
        if (idx >= 0) v[static_cast<size_t>(idx)] += x;
    }
    double&            at(PhononMode m, PhononDirection d, PhononEvent e) { return v[rate_index(m, d, e)]; }
    const double&      at(PhononMode m, PhononDirection d, PhononEvent e) const { return v[rate_index(m, d, e)]; }
    const Rate8&       as_array() const noexcept { return v; }
    static const char* label(int i) {
        static const char* L[8] = {"ac_L_ab", "ac_T_ab", "op_L_ab", "op_T_ab", "ac_L_em", "ac_T_em", "op_L_em", "op_T_em"};
        return (i >= 0 && i < 8) ? L[i] : "invalid";
    }
};

// ---------- Phonon dispersion (uniform lookup + analytic fallback) ----------
struct PhononDispersion {
    PhononMode      mode      = PhononMode::none;
    PhononDirection direction = PhononDirection::none;

    // Analytic ω(q) = w0 + vs*q + c*q^2 (optional fast path)
    double w0 = 0.0, vs = 0.0, c = 0.0;

    // Uniform lookup on q∈[0,qmax]: store only ω; derive indices arithmetically
    double              q0     = 0.0;  // usually 0
    double              inv_dq = 0.0;  // 1/Δq
    double              qmax   = 0.0;
    uint32_t            N      = 0;     // number of samples
    std::vector<double> omega_samples;  // size N

    PhononDispersion() = default;
    PhononDispersion(PhononMode m, PhononDirection d, double w0_, double vs_, double c_) : mode(m), direction(d), w0(w0_), vs(vs_), c(c_) {}

    inline double omega_analytic(double q) const noexcept { return std::fma(c, q * q, std::fma(vs, q, w0)); }

    void build_lookup(double q_max, std::size_t n_points) {
        if (q_max <= 0.0 || n_points < 2) throw std::invalid_argument("bad lookup grid");
        q0              = 0.0;
        qmax            = q_max;
        N               = static_cast<uint32_t>(n_points);
        const double dq = (qmax - q0) / (N - 1);
        inv_dq          = 1.0 / dq;
        omega_samples.resize(N);
        for (uint32_t i = 0; i < N; ++i) {
            const double q   = q0 + i * dq;
            omega_samples[i] = omega_analytic(q);
        }
    }

    inline double omega_lookup(double q) const {
        if (omega_samples.empty()) throw std::runtime_error("phonon lookup empty");
        if (q <= q0) return omega_samples.front();
        if (q >= qmax) return omega_samples.back();
        const double t = (q - q0) * inv_dq;
        uint32_t     i = static_cast<uint32_t>(t);
        if (i >= N - 1) i = N - 2;
        const double frac = t - static_cast<double>(i);
        const double a = omega_samples[i], b = omega_samples[i + 1];
        return std::fma(frac, (b - a), a);
    }

    inline double max_omega() const {
        if (omega_samples.empty()) throw std::runtime_error("phonon lookup empty");
        return *std::max_element(omega_samples.begin(), omega_samples.end());
    }
};

// ---------- Deformation potentials (no i/o in header) ----------
struct DeformationPotential {
    PhononMode mode             = PhononMode::none;
    double     A                = 0.0;
    double     B                = 0.0;
    double     energy_threshold = 1e6;  // eV

    DeformationPotential() = default;
    DeformationPotential(PhononMode m, double A_, double B_, double thr) : mode(m), A(A_), B(B_), energy_threshold(thr) {}

    double get_deformation_potential(const vector3& q, double energy) const {
        const double Ee = (energy < energy_threshold ? energy : energy_threshold);
        return (mode == PhononMode::acoustic) ? std::sqrt(A + Ee * B) * q.norm() : std::sqrt(A + Ee * B);
    }

    double get_fischetti_deformation_potential(const vector3& q, int idx_band) const {
        constexpr double cm_to_m = 1e2;
        const double     boost   = 1.5;
        if (mode == PhononMode::acoustic) {
            return (idx_band == 0 ? boost * 1.2 : 1.0 * 1.7) * q.norm();
        } else {
            return (idx_band == 0 ? boost * 1.75e8 : boost * 2.10e8) * cm_to_m;
        }
    }
};

// ---------- Hole overlap params (fixed arrays instead of linear search) ----------
struct HoleOverlapIntParams {
    // Indices 1..3 used in your logic; we store 0..2 and subtract 1 on access
    std::array<std::array<double, 3>, 3> A{{{{1.0, 3.0, 3.0 / 8.0}}, {{3.0, 1.0, 3.0 / 8.0}}, {{3.0 / 8.0, 3.0 / 8.0, 5.0 / 8.0}}}};
    std::array<std::array<double, 3>, 3> B{{{{3.0, -3.0, 0.0}}, {{-3.0, 3.0, 0.0}}, {{0.0, 0.0, 0.0}}}};

    // returns {A,B} for pair (n1,n2) where n∈{1,2,3}; outside → {0,0}
    std::array<double, 2> get_params(int n1, int n2) const noexcept {
        if (n1 >= 1 && n1 <= 3 && n2 >= 1 && n2 <= 3) {
            return {A[n1 - 1][n2 - 1], B[n1 - 1][n2 - 1]};
        }
        return {0.0, 0.0};
    }
};

// Helper: compact container for 4 (mode×direction) dispersions
inline constexpr int md_index(PhononMode m, PhononDirection d) noexcept {
    const int M = (m == PhononMode::acoustic) ? 0 : (m == PhononMode::optical) ? 1 : -1;
    const int D = (d == PhononDirection::longitudinal) ? 0 : (d == PhononDirection::transverse) ? 1 : -1;
    return (M | D) < 0 ? -1 : ((M << 1) | D);  // 0..3
}

// ---------- Eigen aliases ----------
using EigenSparseMatrix = Eigen::SparseMatrix<double>;
using EigenTriplet      = Eigen::Triplet<double>;

struct Rates_nk_npkp_ctor {
    PhononMode      mode      = PhononMode::none;
    PhononDirection direction = PhononDirection::none;
    PhononEvent     event     = PhononEvent::none;
    /**
     * @brief Preliminary constructor for the rate matrix (n,k) → (n',k') for a given (m,d,e) triplet.
    // (n,k) → (n',k') transition rate matrix. /!\ k' are the barycenters of final tetrahedra, not the mesh vertices!
    // So the size is (nk, n'k') where n'k' are the number of final tetrahedra × number of bands.
     *
     */
    EigenSparseMatrix matrix;  // (nk, n'k')
};

class ElectronPhonon : public BZ_States {
 private:
    double m_temperature        = 300.0;
    double m_rho                = 2.329e3;
    double m_radii_wigner_seitz = 0.0;

    HoleOverlapIntParams m_hole_overlap_int_params;
    DeformationPotential m_ac_defpot_e, m_op_defpot_e;
    DeformationPotential m_ac_defpot_h, m_op_defpot_h;

    // 4 dispersions: (ac/op) × (L/T); index via md_index()
    std::array<PhononDispersion, 4> m_phonon_dispersion;

    // Precomputed rates on mesh: [vertex][band][8]
    std::vector<std::vector<Rate8>> m_list_phonon_scattering_rates;
    std::vector<double>             m_count_weight_tetra_per_vertex;

    // (n,k) → (n',k') transition rate matrices. /!\ k' are the barycenters of final tetrahedra, not the mesh vertices!
    // So the size is (nk, n'k') where n'k' are the number of final tetrahedra × number of bands.
    // One matrix per (m,d,e) triplet.
    std::vector<Rates_nk_npkp_ctor> m_rates_nk_npkp;

 public:
    explicit ElectronPhonon(const EmpiricalPseudopotential::Material& material) : BZ_States(material) {}

    // .cpp implements these (keep header slim—no YAML/iostream here)
    void   load_phonon_parameters(const std::string& filename);
    void   plot_phonon_dispersion(const std::string& filename) const;
    double get_max_phonon_energy() const;

    inline double bose_einstein_distribution(double energy_eV, double temperature_K) const;  // inline in .cpp
    double        electron_overlap_integral(const vector3& k1, const vector3& k2) const;
    double        hole_overlap_integral(int n1, const vector3& k1, int n2, const vector3& k2) const;

    Rate8      compute_transition_rates_pair(int idx_n1, std::size_t idx_k1, int idx_n2, std::size_t idx_tetra_final, bool push_nk_npkp);
    RateValues compute_electron_phonon_rate(int idx_n1, std::size_t idx_k1, bool populate_nk_npkp = false);
    RateValues compute_hole_phonon_rate(int idx_n1, std::size_t idx_k1);

    void set_temperature(double T) noexcept { m_temperature = T; }
    void set_density(double rho) noexcept { m_rho = rho; }

    void compute_electron_phonon_rates_over_mesh(double energy_max = 100.0, bool irreducible_wedge_only = false, bool populate_nk_npkp = false);
    void add_electron_phonon_rates_to_mesh(const std::string& initial_filename, const std::string& final_filename);
    void compute_electron_phonon_rates_over_mesh_nk_npkp(bool irreducible_wedge_only = false);

    std::pair<int, std::size_t> select_final_state(std::size_t     idx_band_initial,
                                                   std::size_t     idx_k_initial,
                                                   PhononMode      mode,
                                                   PhononDirection direction,
                                                   PhononEvent     event) const;

    void export_rate_values(const std::string& filename) const;

    void compute_plot_electron_phonon_rates_vs_energy_over_mesh(int                nb_bands,
                                                                double             max_energy,
                                                                double             energy_step,
                                                                const std::string& filename,
                                                                bool               irreducible_wedge_only = false);

    void          read_phonon_scattering_rates_from_file(const std::filesystem::path& path);
    Rate8         interpolate_phonon_scattering_rate_at_location(const vector3& location, const std::size_t& idx_band) const;
    inline double sum_modes(const Rate8& r) const noexcept;
    double        compute_P_Gamma() const;

    void compute_RTA_mobility();
};

}  // namespace bz_mesh
