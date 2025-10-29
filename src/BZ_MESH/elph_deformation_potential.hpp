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

#include <cmath>

#include "vector.hpp"
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

    /**
     * @brief Get the fischetti deformation potential object
     *
     * Values of Fischetti et al. for Si :
     *  if (mode == PhononMode::acoustic) {
     *      return (idx_band == 0 ? boost_acc * 1.2 : boost_acc * 1.0 * 1.7) * norm_q;
     *   } else {
     *       return (idx_band == 0 ? boost_opt * 1.75e8 : boost_opt * 2.10e8) * cm_to_m;
     *   }
     * Values from XANG : . Appl. Phys. 73, 3339–3347 (1993) https://doi.org/10.1063/1.352959
     *
     *  if (mode == PhononMode::acoustic) {
            return (idx_band == 0 ? boost_acc * 1.8 : boost_acc * 1.0 * 2.5) * norm_q;
        } else {
            return (idx_band == 0 ? boost_opt * 3.4e8 : boost_opt * 4.10e8) * cm_to_m;
        }
     *
     * @param norm_q
     * @param idx_band
     * @return double
     */
    double get_fischetti_deformation_potential(double norm_q, int idx_band) const {
        constexpr double cm_to_m = 1e2;
        const double     boost_acc   = 1.25;
        const double     boost_opt   = 1.25;
        if (mode == PhononMode::acoustic) {
            return (idx_band == 0 ? boost_acc * 1.8 : boost_acc * 2.5) * norm_q;
        } else {
            return (idx_band == 0 ? boost_opt * 3.4e8 : boost_opt * 4.10e8) * cm_to_m;
        }
    }
};

}  // namespace uepm::mesh_bz