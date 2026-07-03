/**
 * @file screened_coulomb_kernel.hpp
 * @brief Screened Coulomb interaction helper for ab initio impact ionization.
 */

#pragma once

#include <complex>

#include "dielectric_mesh.hpp"

namespace uepm::mesh_bz {

struct ScreenedCoulombKernelConfig {
    double m_min_q_norm_SI       = 1e6;    // 1/m, avoids the singular q = 0 Fourier component.
    double m_min_abs_epsilon     = 1e-12;  // dimensionless
    bool   m_throw_on_singular_q = true;
};

/**
 * @brief Evaluates e^2 / (eps0 * epsilon(q, E) * |q|^2).
 *
 * The returned interaction does not include the plane-wave normalization volume.
 * Its unit is eV*m^3, so multiplying by 1/Omega gives an energy in eV.
 */
class ScreenedCoulombKernel {
 private:
    const DielectricMesh*        m_dielectric_mesh = nullptr;
    ScreenedCoulombKernelConfig  m_config{};

 public:
    explicit ScreenedCoulombKernel(const DielectricMesh& dielectric_mesh,
                                   ScreenedCoulombKernelConfig config = {});

    const ScreenedCoulombKernelConfig& config() const noexcept { return m_config; }

    complex_d dielectric_function(const vector3& q_SI, double energy_eV) const;
    complex_d interaction_eV_m3(const vector3& q_SI, double energy_eV) const;
};

}  // namespace uepm::mesh_bz
