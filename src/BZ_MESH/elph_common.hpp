/**
 * @file elph_common.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2025-10-10
 *
 *
 */

// phonon_axes.hpp
#pragma once
#include <array>
#include <cstdint>

namespace uepm::mesh_bz {

// ---------- Discrete axes ----------
enum class PhononMode : uint8_t { acoustic = 0, optical = 1, none = 2 };
enum class PhononDirection : uint8_t { longitudinal = 0, transverse = 1, none = 2 };
enum class PhononEvent : uint8_t { absorption = 0, emission = 1, none = 2 };

// ---------- Compact indices ----------
// Channel packing: (M,D,E) -> [0..7] with (m<<2)|(d<<1)|e
constexpr int rate_index(PhononMode m, PhononDirection d, PhononEvent e) noexcept {
    const int M = (m == PhononMode::acoustic) ? 0 : (m == PhononMode::optical) ? 1 : -1;
    const int D = (d == PhononDirection::longitudinal) ? 0 : (d == PhononDirection::transverse) ? 1 : -1;
    const int E = (e == PhononEvent::absorption) ? 0 : (e == PhononEvent::emission) ? 1 : -1;
    if ((M | D | E) < 0) return -1;
    return (E << 2) | (M << 1) | D;  // bit2=E, bit1=M, bit0=D
}

// Optional named indices (kept here so they can’t drift)
constexpr int IDX_AC_L_AB = 0;
constexpr int IDX_AC_T_AB = 1;
constexpr int IDX_OP_L_AB = 2;
constexpr int IDX_OP_T_AB = 3;
constexpr int IDX_AC_L_EM = 4;
constexpr int IDX_AC_T_EM = 5;
constexpr int IDX_OP_L_EM = 6;
constexpr int IDX_OP_T_EM = 7;

// Helper: 4-branch (mode×direction) index
inline constexpr int md_index(PhononMode m, PhononDirection d) noexcept {
    const int M = (m == PhononMode::acoustic) ? 0 : (m == PhononMode::optical) ? 1 : -1;
    const int D = (d == PhononDirection::longitudinal) ? 0 : (d == PhononDirection::transverse) ? 1 : -1;
    return (M | D) < 0 ? -1 : ((M << 1) | D);  // 0..3 or -1 if invalid
}

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
    double&       at(PhononMode m, PhononDirection d, PhononEvent e) { return v[rate_index(m, d, e)]; }
    const double& at(PhononMode m, PhononDirection d, PhononEvent e) const { return v[rate_index(m, d, e)]; }
    const Rate8&  as_array() const noexcept { return v; }

    static const char* label(int i) noexcept {
        static constexpr const char* L[8] = {"ac_L_ab", "ac_T_ab", "op_L_ab", "op_T_ab", "ac_L_em", "ac_T_em", "op_L_em", "op_T_em"};
        return (i >= 0 && i < 8) ? L[i] : "invalid";
    }
};

/**
 * @brief Compute the safe cosine of the angle between two vectors.
 *
 * @param v First vector.
 * @param w Second vector.
 * @return double Cosine of the angle between v and w, or 0 if either is nearly zero.
 */
[[nodiscard]] inline double cos_angle_safe(const vector3& v, const vector3& w) noexcept {
    const double nv2 = v.norm_squared();
    const double nw2 = w.norm_squared();
    if (nv2 <= 0.0 || nw2 <= 0.0) return 0.0;
    double c = v.dot(w) / std::sqrt(nv2 * nw2);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    return c;
}

inline double fermi_dirac_distribution(double energy_eV, double fermi_level_eV, double temperature_K) {
    const double kT = uepm::Constants::k_b_eV * temperature_K;

    // Handle T <= 0 K as the T→0 limit (Heaviside step).
    if (!(kT > 0.0)) {
        if (energy_eV < fermi_level_eV) return 1.0;
        if (energy_eV > fermi_level_eV) return 0.0;
        return 0.5;  // at E = Ef, take the symmetric limit
    }

    const double x = (energy_eV - fermi_level_eV) / kT;

    // In double precision, |x| >= 40 puts f within ~1e-17 of 0 or 1.
    constexpr double X_CUTOFF = 40.0;
    if (x >= X_CUTOFF) return 0.0;
    if (x <= -X_CUTOFF) return 1.0;

    // Stable logistic: avoid large exp(x) when x > 0.
    if (x > 0.0) {
        const double emx = std::exp(-x);
        return emx / (1.0 + emx);
    } else {
        const double ex = std::exp(x);
        return 1.0 / (1.0 + ex);
    }
}

inline double d_de_fermi_dirac_dE(double energy_eV, double fermi_level_eV, double temperature_K) {
    const double kT = uepm::Constants::k_b_eV * temperature_K;

    // T <= 0: derivative is a Dirac delta in theory; return 0 numerically.
    if (!(kT > 0.0)) return 0.0;

    const double x = (energy_eV - fermi_level_eV) / kT;

    constexpr double X_CUTOFF = 40.0;
    if (x >= X_CUTOFF || x <= -X_CUTOFF) return 0.0;

    // Compute f stably, then use df/dE = -(1/kT) * f * (1 - f).
    double f;
    if (x > 0.0) {
        const double emx = std::exp(-x);
        f                = emx / (1.0 + emx);
    } else {
        const double ex = std::exp(x);
        f               = 1.0 / (1.0 + ex);
    }
    return -(f * (1.0 - f)) / kT;
}

inline double bose_einstein_distribution(double energy_eV, double temperature_K)  {
    // N0 = 1 / (exp(E / kT) - 1)
    const double x = energy_eV / (uepm::Constants::k_b_eV * temperature_K);
    return 1.0 / std::expm1(x);  // stable for small x
}

/**
 * @brief Transport weight factor 1 - cos(θ) where θ is the angle between v0 and v1, the group velocities.
 *
 * @param v0
 * @param v1
 * @return double
 */
inline double transport_weight_RTA(const vector3& v0, const vector3& v1) {
    const double norm_v0 = v0.norm_squared();
    const double norm_v1 = v1.norm_squared();
    if (norm_v0 < 1e-24 || norm_v1 < 1e-24) return 1.0;  // degeneracy/edge guard
    double cos_theta = v0.dot(v1) / std::sqrt(norm_v0 * norm_v1);
    if (cos_theta > 1.0) cos_theta = 1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;
    return 1.0 - cos_theta;  // = 1 - cos θ
}

// ---- Optional compile-time sanity checks (won’t cost runtime) ----
static_assert(rate_index(PhononMode::acoustic, PhononDirection::longitudinal, PhononEvent::absorption) == 0);
static_assert(rate_index(PhononMode::acoustic, PhononDirection::transverse, PhononEvent::absorption) == 1);
static_assert(rate_index(PhononMode::optical, PhononDirection::longitudinal, PhononEvent::absorption) == 2);
static_assert(rate_index(PhononMode::optical, PhononDirection::transverse, PhononEvent::absorption) == 3);
static_assert(rate_index(PhononMode::acoustic, PhononDirection::longitudinal, PhononEvent::emission) == 4);
static_assert(rate_index(PhononMode::acoustic, PhononDirection::transverse, PhononEvent::emission) == 5);
static_assert(rate_index(PhononMode::optical, PhononDirection::longitudinal, PhononEvent::emission) == 6);
static_assert(rate_index(PhononMode::optical, PhononDirection::transverse, PhononEvent::emission) == 7);

}  // namespace uepm::mesh_bz
