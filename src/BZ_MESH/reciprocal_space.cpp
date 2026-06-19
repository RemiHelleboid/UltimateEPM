/**
 * @file reciprocal_space.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-18
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "reciprocal_space.hpp"

#include <Eigen/LU>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "bz_mesh.hpp"

namespace uepm::mesh_bz {

namespace {

inline bool contains_bcc_reduced(const vector3& k, double h) noexcept {
    constexpr double eps = 1e-12;

    const double ax = std::abs(k.x());
    const double ay = std::abs(k.y());
    const double az = std::abs(k.z());

    return ax <= h + eps && ay <= h + eps && az <= h + eps && ax + ay + az <= 1.5 * h + eps;
}

inline double nearest_even(double x) noexcept { return 2.0 * std::nearbyint(0.5 * x); }

inline double nearest_odd(double x) noexcept { return 1.0 + 2.0 * std::nearbyint(0.5 * (x - 1.0)); }

inline vector3 fold_bcc_reduced_fast(const vector3& k, double h) noexcept {
    const double x = k.x() / h;
    const double y = k.y() / h;
    const double z = k.z() / h;

    const double ge_x = nearest_even(x);
    const double ge_y = nearest_even(y);
    const double ge_z = nearest_even(z);

    const double fe_x = x - ge_x;
    const double fe_y = y - ge_y;
    const double fe_z = z - ge_z;
    const double de2  = fe_x * fe_x + fe_y * fe_y + fe_z * fe_z;

    const double go_x = nearest_odd(x);
    const double go_y = nearest_odd(y);
    const double go_z = nearest_odd(z);

    const double fo_x = x - go_x;
    const double fo_y = y - go_y;
    const double fo_z = z - go_z;
    const double do2  = fo_x * fo_x + fo_y * fo_y + fo_z * fo_z;

    if (de2 <= do2) {
        return {h * fe_x, h * fe_y, h * fe_z};
    }

    return {h * fo_x, h * fo_y, h * fo_z};
}

}  // namespace

void ReciprocalSpace::initialize_basis(const Eigen::Vector3d& b1_SI,
                                       const Eigen::Vector3d& b2_SI,
                                       const Eigen::Vector3d& b3_SI,
                                       double                 halfwidth_reduced,
                                       double                 si_to_reduced) {
    if (!(si_to_reduced > 0.0)) {
        throw std::invalid_argument("SI-to-reduced reciprocal-space scale must be positive");
    }
    m_basis.col(0)  = b1_SI;
    m_basis.col(1)  = b2_SI;
    m_basis.col(2)  = b3_SI;
    m_inverse_basis = m_basis.inverse();
    m_halfwidth     = halfwidth_reduced;
    m_si_to_reduced = si_to_reduced;
}

void ReciprocalSpace::precompute_bcc_shifts(double si_to_reduced, int maximum_shell) {
    if (!(si_to_reduced > 0.0) || maximum_shell < 1) {
        throw std::invalid_argument("Invalid reciprocal-space shift configuration");
    }

    m_reciprocal_shifts.clear();
    m_reciprocal_shifts.push_back({0.0, 0.0, 0.0});

    vector3 b1{-1.0, 1.0, 1.0};
    vector3 b2{1.0, -1.0, 1.0};
    vector3 b3{1.0, 1.0, -1.0};
    b1 /= si_to_reduced;
    b2 /= si_to_reduced;
    b3 /= si_to_reduced;

    for (int shell = 1; shell <= maximum_shell; ++shell) {
        for (int n1 = -shell; n1 <= shell; ++n1) {
            for (int n2 = -shell; n2 <= shell; ++n2) {
                for (int n3 = -shell; n3 <= shell; ++n3) {
                    if (std::abs(n1) + std::abs(n2) + std::abs(n3) == shell) {
                        m_reciprocal_shifts.push_back(n1 * b1 + n2 * b2 + n3 * b3);
                    }
                }
            }
        }
    }
}

bool ReciprocalSpace::contains_bcc(const vector3& k_SI, double si_to_reduced) const noexcept {
    constexpr double eps = 1e-12;
    const double     ax  = std::abs(k_SI.x() * si_to_reduced);
    const double     ay  = std::abs(k_SI.y() * si_to_reduced);
    const double     az  = std::abs(k_SI.z() * si_to_reduced);

    return ax <= m_halfwidth + eps && ay <= m_halfwidth + eps && az <= m_halfwidth + eps &&
           ax + ay + az <= 1.5 * m_halfwidth + eps;
}

// vector3 ReciprocalSpace::retrieve_bcc_image(const vector3& k_SI, double si_to_reduced) const {
//     for (const vector3& shift : m_reciprocal_shifts) {
//         const vector3 candidate = k_SI + shift;
//         if (contains_bcc(candidate, si_to_reduced)) {
//             return candidate;
//         }
//     }

//     const vector3 folded = fold_wigner_seitz(k_SI);
//     if (contains_bcc(folded, si_to_reduced)) {
//         return folded;
//     }
//     throw std::runtime_error("No reciprocal-lattice image lies inside the Brillouin zone");
// }

vector3 ReciprocalSpace::fold_bcc_fast_SI(const vector3& k_SI, double si_to_reduced) const {
    if (!(si_to_reduced > 0.0)) {
        throw std::invalid_argument("fold_bcc_fast_SI: si_to_reduced must be positive");
    }

    const vector3 k_reduced{k_SI.x() * si_to_reduced, k_SI.y() * si_to_reduced, k_SI.z() * si_to_reduced};
    if (contains_bcc_reduced(k_reduced, m_halfwidth)) {
        return k_SI;
    }

    const vector3 folded_reduced = fold_bcc_reduced_fast(k_reduced, m_halfwidth);
    assert(contains_bcc_reduced(folded_reduced, m_halfwidth));
    return {folded_reduced.x() / si_to_reduced, folded_reduced.y() / si_to_reduced, folded_reduced.z() / si_to_reduced};
}

vector3 ReciprocalSpace::fold_bcc_fast_SI(const vector3& k_SI) const {
    return fold_bcc_fast_SI(k_SI, m_si_to_reduced);
}

vector3 ReciprocalSpace::retrieve_bcc_image(const vector3& k_SI, double si_to_reduced) const {
    return fold_bcc_fast_SI(k_SI, si_to_reduced);
}

vector3 ReciprocalSpace::fold_wigner_seitz(const vector3& k_SI) const noexcept {
    const Eigen::Vector3d k(k_SI.x(), k_SI.y(), k_SI.z());
    const Eigen::Vector3d reduced = m_inverse_basis * k;
    const Eigen::Vector3d nearest = reduced.array().round().matrix();

    Eigen::Vector3d best              = k - m_basis * nearest;
    double          best_norm_squared = best.squaredNorm();

    for (int i = -2; i <= 2; ++i) {
        for (int j = -2; j <= 2; ++j) {
            for (int l = -2; l <= 2; ++l) {
                const Eigen::Vector3d lattice_coordinates =
                    nearest + Eigen::Vector3d(static_cast<double>(i), static_cast<double>(j), static_cast<double>(l));
                const Eigen::Vector3d candidate = k - m_basis * lattice_coordinates;
                const double          norm      = candidate.squaredNorm();
                if (norm < best_norm_squared) {
                    best_norm_squared = norm;
                    best              = candidate;
                }
            }
        }
    }
    return {best.x(), best.y(), best.z()};
}

bool ReciprocalSpace::inside_wigner_seitz(const vector3& k_SI, double si_to_reduced) const noexcept {
    return contains_bcc(k_SI, si_to_reduced);
}

bool ReciprocalSpace::inside_irreducible_wedge(const vector3& k_SI, double si_to_reduced) const noexcept {
    const double     x   = k_SI.x() * si_to_reduced;
    const double     y   = k_SI.y() * si_to_reduced;
    const double     z   = k_SI.z() * si_to_reduced;
    constexpr double eps = 1e-12;
    return z >= -eps && z <= y + eps && y <= x + eps && x <= m_halfwidth + eps && x + y + z <= 1.5 * m_halfwidth + eps;
}

bool MeshBZ::is_inside_mesh_geometry(const vector3& k) const {
    return m_reciprocal_space.contains_bcc(k, si_to_reduced_scale());
}

void MeshBZ::precompute_G_shifts() { m_reciprocal_space.precompute_bcc_shifts(si_to_reduced_scale()); }

vector3 MeshBZ::retrieve_k_inside_mesh_geometry(const vector3& k) const {
    return m_reciprocal_space.fold_bcc_fast_SI(k);
}

void MeshBZ::init_reciprocal_basis(const Eigen::Vector3d& b1_SI,
                                   const Eigen::Vector3d& b2_SI,
                                   const Eigen::Vector3d& b3_SI,
                                   double                 halfwidth_reduced,
                                   double                 si_to_reduced) {
    if (!(si_to_reduced > 0.0)) {
        throw std::invalid_argument("SI-to-reduced reciprocal-space scale must be positive");
    }
    m_reciprocal_space.initialize_basis(b1_SI, b2_SI, b3_SI, halfwidth_reduced, si_to_reduced);
}

vector3 MeshBZ::fold_ws_bcc(const vector3& k_SI) const noexcept { return m_reciprocal_space.fold_wigner_seitz(k_SI); }

bool MeshBZ::inside_ws_bcc(const vector3& k_SI) const noexcept {
    return m_reciprocal_space.inside_wigner_seitz(k_SI, si_to_reduced_scale());
}

bool MeshBZ::is_irreducible_wedge(const vector3& k_SI) const noexcept {
    return m_reciprocal_space.inside_irreducible_wedge(k_SI, si_to_reduced_scale());
}

}  // namespace uepm::mesh_bz
