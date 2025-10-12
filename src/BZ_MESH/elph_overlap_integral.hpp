// elph_overlap_integral.hpp
#pragma once

#include <array>
#include <cmath>

#include "elph_common.hpp"
#include "vector.hpp"

namespace uepm::mesh_bz {

// Parameters for hole-hole overlap integrals (HH/LH/SO style).
// Units: dimensionless (fits are empirical).
struct HoleOverlapIntParams {
    // Indices 1..3 used in legacy logic; stored 0..2.
    std::array<std::array<double, 3>, 3> A{std::array<double, 3>{1.0, 3.0, 3.0 / 8.0},
                                           std::array<double, 3>{3.0, 1.0, 3.0 / 8.0},
                                           std::array<double, 3>{3.0 / 8.0, 3.0 / 8.0, 5.0 / 8.0}};
    std::array<std::array<double, 3>, 3> B{std::array<double, 3>{3.0, -3.0, 0.0},
                                           std::array<double, 3>{-3.0, 3.0, 0.0},
                                           std::array<double, 3>{0.0, 0.0, 0.0}};

    // Return {A,B} for pair (n1,n2) where n ∈ {1,2,3}; outside → {0,0}
    [[nodiscard]] std::array<double, 2> get_params_one_based(int n1, int n2) const noexcept {
        if (n1 >= 1 && n1 <= 3 && n2 >= 1 && n2 <= 3) {
            return {A[static_cast<size_t>(n1 - 1)][static_cast<size_t>(n2 - 1)],
                    B[static_cast<size_t>(n1 - 1)][static_cast<size_t>(n2 - 1)]};
        }
        return {0.0, 0.0};
    }
};


/**
 * Electron overlap integral.
 * Model: spherical j1 kernel: 3 (sin x − x cos x) / x^3 with x = |k1−k2| * R_WS.
 * k1,k2 in SI [1/m] if R_WS is in [m]. Result: dimensionless.
 */
[[nodiscard]] inline double electron_overlap_integral(const vector3& k1, const vector3& k2, double R_Wigner_Seitz_m) noexcept {
    const double x = (k1 - k2).norm() * R_Wigner_Seitz_m;

    // Series for small x: 1 − x^2/10 + O(x^4). Use a slightly larger threshold than 1e-6 for safety.
    if (std::abs(x) < 1e-4) {
        const double x2 = x * x;
        return 1.0 - 0.1 * x2;
    }
    // Guard against overflow in denominator (extremely unlikely with physical k)
    const double x2 = x * x;
    const double x3 = x2 * x;
    return 3.0 * (std::sin(x) - x * std::cos(x)) / x3;
}

/**
 * Hole overlap integral between bands n1 and n2 (1..3) at k1 and k2.
 * Uses params A,B: I = 0.5 * sqrt( A + B * cos^2(theta) ), clamped to R>=0 for FP safety.
 * Returns dimensionless.
 */
[[nodiscard]] inline double hole_overlap_integral(int                         n1,
                                                  const vector3&              k1,
                                                  int                         n2,
                                                  const vector3&              k2,
                                                  const HoleOverlapIntParams& params) noexcept {
    const auto   AB  = params.get_params_one_based(n1, n2);
    const double c   = cos_angle_safe(k1, k2);
    const double c2  = c * c;
    double       rad = AB[0] + AB[1] * c2;
    if (rad < 0.0) rad = 0.0;  // clamp due to FP noise
    return 0.5 * std::sqrt(rad);
}

}  // namespace uepm::mesh_bz
