/**
 * @file elph_deformation_potential.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2025-10-10
 *
 *
 */

#pragma once

#include "elph_common.hpp"

namespace uepm::mesh_bz {

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
        const double     boost   = 1.6;
        if (mode == PhononMode::acoustic) {
            return (idx_band == 0 ? boost * 1.2 : boost * 1.0 * 1.7) * q.norm();
        } else {
            return (idx_band == 0 ? boost * 1.75e8 : boost * 2.10e8) * cm_to_m;
        }
    }
};

}  // namespace uepm::mesh_bz