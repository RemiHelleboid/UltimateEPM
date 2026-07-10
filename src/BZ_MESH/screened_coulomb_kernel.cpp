/**
 * @file screened_coulomb_kernel.cpp
 * @brief Screened Coulomb interaction helper for ab initio impact ionization.
 */

#include "screened_coulomb_kernel.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "physical_constants.hpp"

namespace uepm::mesh_bz {

namespace {

void require_finite_vector(const vector3& q_SI) {
    if (!std::isfinite(q_SI.x()) || !std::isfinite(q_SI.y()) || !std::isfinite(q_SI.z())) {
        throw std::invalid_argument("ScreenedCoulombKernel: q vector must be finite.");
    }
}

}  // namespace

ScreenedCoulombKernel::ScreenedCoulombKernel(const DielectricMesh& dielectric_mesh, ScreenedCoulombKernelConfig config)
    : m_dielectric_mesh(&dielectric_mesh),
      m_config(config) {
    if (!(m_config.m_min_q_norm_SI >= 0.0) || !std::isfinite(m_config.m_min_q_norm_SI)) {
        throw std::invalid_argument("ScreenedCoulombKernel: minimum q norm must be finite and non-negative.");
    }
    if (!(m_config.m_min_abs_epsilon > 0.0) || !std::isfinite(m_config.m_min_abs_epsilon)) {
        throw std::invalid_argument("ScreenedCoulombKernel: minimum epsilon magnitude must be finite and positive.");
    }
}

complex_d ScreenedCoulombKernel::dielectric_function(const vector3& q_SI, double energy_eV) const {
    if (!m_dielectric_mesh) {
        throw std::logic_error("ScreenedCoulombKernel: missing dielectric mesh.");
    }
    require_finite_vector(q_SI);
    if (!std::isfinite(energy_eV)) {
        throw std::invalid_argument("ScreenedCoulombKernel: energy must be finite.");
    }
    const double q_norm = q_SI.norm();
    if (q_norm < m_config.m_min_q_norm_SI && m_config.m_throw_on_singular_q) {
        throw std::domain_error("ScreenedCoulombKernel: q is too close to the Coulomb singularity.");
    }

    complex_d epsilon = m_dielectric_mesh->interpolate_dielectric_function(q_SI, energy_eV);
    if (!std::isfinite(epsilon.real()) || !std::isfinite(epsilon.imag())) {
        throw std::runtime_error("ScreenedCoulombKernel: dielectric interpolation returned a non-finite value.");
    }
    if (std::abs(epsilon) < m_config.m_min_abs_epsilon) {
        throw std::domain_error("ScreenedCoulombKernel: dielectric function magnitude is too small.");
    }
    return epsilon;
}

complex_d ScreenedCoulombKernel::interaction_eV_m3(const vector3& q_SI, double energy_eV) const {
    require_finite_vector(q_SI);
    const double q_norm_squared = q_SI.norm_squared();
    const double min_q_squared  = m_config.m_min_q_norm_SI * m_config.m_min_q_norm_SI;
    if (q_norm_squared < min_q_squared) {
        if (m_config.m_throw_on_singular_q) {
            throw std::domain_error("ScreenedCoulombKernel: q is too close to the Coulomb singularity.");
        }
    }
    const double regularized_q_squared = std::max(q_norm_squared, min_q_squared);
    const auto   epsilon               = dielectric_function(q_SI, energy_eV);
    const double prefactor_eV_m3       = (uepm::constants::q_e * uepm::constants::q_e) /
                                   (uepm::constants::eps_0 * regularized_q_squared * uepm::constants::q_e);
    return prefactor_eV_m3 / epsilon;
}

}  // namespace uepm::mesh_bz
