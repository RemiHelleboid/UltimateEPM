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

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "elph_common.hpp"
#include "vector_bz.hpp"

namespace uepm::mesh_bz {

struct DeformationPotential {
    PhononMode mode             = PhononMode::none;
    double     A                = 0.0;
    double     B                = 0.0;
    double     energy_threshold = 1e6;  // eV

    DeformationPotential() = default;
    DeformationPotential(PhononMode m, double A_, double B_, double thr)
        : mode(m),
          A(A_),
          B(B_),
          energy_threshold(thr) {}

    double get_deformation_potential(const vector3& q, double energy_eV) const {
        if (!std::isfinite(energy_eV)) {
            throw std::invalid_argument("Deformation-potential energy must be finite");
        }
        if (mode != PhononMode::acoustic && mode != PhononMode::optical) {
            throw std::logic_error("Deformation-potential phonon mode is not configured");
        }

        const double effective_energy_eV = std::min(energy_eV, energy_threshold);
        const double squared_magnitude   = A + effective_energy_eV * B;
        if (!(squared_magnitude >= 0.0) || !std::isfinite(squared_magnitude)) {
            throw std::domain_error("Deformation-potential A + B*E is negative or non-finite");
        }

        const double magnitude = std::sqrt(squared_magnitude);
        return mode == PhononMode::acoustic ? magnitude * q.norm() : magnitude;
    }
};

}  // namespace uepm::mesh_bz
