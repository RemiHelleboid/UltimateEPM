/**
 * @file valley_model.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-03-27
 *
 *
 */

/**
 * @file valley_model.hpp
 * @brief Analytical anisotropic Kane valley model for electron transport.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

#include "physical_constants.hpp"
#include "vector.hpp"

namespace uepm::amc {

class valley_model {
 public:
    using vector3 = uepm::mesh_bz::vector3;
    using mat3    = std::array<std::array<double, 3>, 3>;
    enum class band_type { parabolic, kane };

    struct parameters {
        double      transverse_effective_mass   = 0.0;
        double      longitudinal_effective_mass = 0.0;
        double      non_parabolicity            = 0.0;  // [1 / eV]
        double      energy_offset               = 0.0;  // [eV]
        double      phonon_reference_energy     = 0.0;  // [eV], optional
        std::size_t degeneracy                  = 1;
        band_type   dispersion                  = band_type::kane;
        mat3        rotation                    = identity_matrix();
        std::string name                        = "unnamed_valley";
    };

 private:
    double      m_transverse_effective_mass   = 0.0;
    double      m_longitudinal_effective_mass = 0.0;
    double      m_non_parabolicity            = 0.0;
    double      m_energy_offset               = 0.0;
    double      m_phonon_reference_energy     = 0.0;
    std::size_t m_degeneracy                  = 1;
    band_type   m_dispersion                  = band_type::kane;
    mat3        m_rotation                    = identity_matrix();
    std::string m_name                        = "unnamed_valley";


 public:
    valley_model()                                   = default;
    valley_model(const valley_model&)                = default;
    valley_model(valley_model&&) noexcept            = default;
    valley_model& operator=(const valley_model&)     = default;
    valley_model& operator=(valley_model&&) noexcept = default;
    ~valley_model()                                  = default;

    explicit valley_model(const parameters& p)
        : m_transverse_effective_mass(p.transverse_effective_mass),
          m_longitudinal_effective_mass(p.longitudinal_effective_mass),
          m_non_parabolicity(p.non_parabolicity),
          m_energy_offset(p.energy_offset),
          m_phonon_reference_energy(p.phonon_reference_energy),
          m_degeneracy(p.degeneracy),
          m_dispersion(p.dispersion),
          m_rotation(p.rotation),
          m_name(p.name) {
        validate();
    }

    const std::string& name() const noexcept { return m_name; }

    double transverse_effective_mass() const noexcept { return m_transverse_effective_mass; }

    double longitudinal_effective_mass() const noexcept { return m_longitudinal_effective_mass; }

    double non_parabolicity() const noexcept { return m_non_parabolicity; }

    double energy_offset() const noexcept { return m_energy_offset; }

    double phonon_reference_energy() const noexcept { return m_phonon_reference_energy; }

    std::size_t degeneracy() const noexcept { return m_degeneracy; }

    band_type dispersion() const noexcept { return m_dispersion; }

    const mat3& rotation() const noexcept { return m_rotation; }

    bool is_isotropic(double tolerance = 1e-18) const noexcept {
        return std::abs(m_longitudinal_effective_mass - m_transverse_effective_mass) < tolerance;
    }

    double conductivity_effective_mass() const noexcept {
        return 3.0 / (1.0 / m_longitudinal_effective_mass + 2.0 / m_transverse_effective_mass);
    }

    double density_of_states_effective_mass() const noexcept {
        return std::cbrt(m_longitudinal_effective_mass * m_transverse_effective_mass * m_transverse_effective_mass);
    }

    vector3 to_valley_frame(const vector3& k_global) const noexcept { return multiply(transpose(m_rotation), k_global); }

    vector3 to_global_frame(const vector3& k_valley) const noexcept { return multiply(m_rotation, k_valley); }

    /**
     * @brief Quadratic form:
     *        gamma(k) = (hbar^2 / 2) * (kt1^2 / mt + kt2^2 / mt + kl^2 / ml)
     *        where k is expressed in the valley principal frame.
     */
    double gamma_from_k_valley(const vector3& k_valley) const noexcept {
        const double kt1 = k_valley.x();
        const double kt2 = k_valley.y();
        const double kl  = k_valley.z();

        const double gamma_joule = 0.5 * uepm::constants::h_bar * uepm::constants::h_bar *
                                   ((kt1 * kt1) / m_transverse_effective_mass + (kt2 * kt2) / m_transverse_effective_mass +
                                    (kl * kl) / m_longitudinal_effective_mass);

        return gamma_joule / uepm::constants::eV_to_J;
    }

    double gamma_from_k(const vector3& k_global) const noexcept { return gamma_from_k_valley(to_valley_frame(k_global)); }

    /**
     * @brief Kinetic energy from gamma using:
     *        gamma = E                    for parabolic
     *        gamma = E * (1 + alpha * E)  for Kane
     */
    double kinetic_energy_from_gamma(double gamma) const {
        if (gamma < 0.0) {
            throw std::invalid_argument("gamma must be non-negative");
        }

        if (m_dispersion == band_type::parabolic || m_non_parabolicity == 0.0) {
            return gamma;
        }

        const double alpha = m_non_parabolicity;
        const double disc  = 1.0 + 4.0 * alpha * gamma;
        if (disc < 0.0) {
            throw std::runtime_error("negative discriminant in Kane energy inversion");
        }

        return (std::sqrt(disc) - 1.0) / (2.0 * alpha);
    }

    double gamma_from_kinetic_energy(double energy) const {
        if (energy < 0.0) {
            throw std::invalid_argument("kinetic energy must be non-negative");
        }

        if (m_dispersion == band_type::parabolic || m_non_parabolicity == 0.0) {
            return energy;
        }

        return energy * (1.0 + m_non_parabolicity * energy);
    }

    double kinetic_energy_from_k(const vector3& k_global) const { return kinetic_energy_from_gamma(gamma_from_k(k_global)); }

    double total_energy_from_k(const vector3& k_global) const { return m_energy_offset + kinetic_energy_from_k(k_global); }

    /**
     * @brief Group velocity in valley frame:
     *        v_i = (1 / hbar) * dE/dk_i
     *            = hbar * k_i / (m_i * (1 + 2 alpha E))
     */
    vector3 velocity_from_k_valley(const vector3& k_valley) const {
        const double energy = kinetic_energy_from_gamma(gamma_from_k_valley(k_valley));
        const double factor = 1.0 + 2.0 * m_non_parabolicity * energy;

        return vector3{uepm::constants::h_bar * k_valley.x() / (m_transverse_effective_mass * factor),
                       uepm::constants::h_bar * k_valley.y() / (m_transverse_effective_mass * factor),
                       uepm::constants::h_bar * k_valley.z() / (m_longitudinal_effective_mass * factor)};
    }

    vector3 velocity_from_k(const vector3& k_global) const { return to_global_frame(velocity_from_k_valley(to_valley_frame(k_global))); }

    /**
     * @brief Build a valley-frame wave-vector from a kinetic energy and a direction
     *        in the scaled ellipsoidal coordinate system.
     *
     * direction must be normalized.
     */
    vector3 k_valley_from_energy_direction(double energy, const vector3& direction) const {
        if (energy < 0.0) {
            throw std::invalid_argument("kinetic energy must be non-negative");
        }

        const double norm = direction.norm();
        if (norm <= std::numeric_limits<double>::epsilon()) {
            throw std::invalid_argument("direction must be non-zero");
        }

        const vector3 u           = direction / norm;
        const double  gamma_eV    = gamma_from_kinetic_energy(energy);
        const double  gamma_joule = gamma_eV * uepm::constants::eV_to_J;
        const double  pref        = std::sqrt(2.0 * gamma_joule) / uepm::constants::h_bar;

        return vector3{pref * std::sqrt(m_transverse_effective_mass) * u.x(),
                       pref * std::sqrt(m_transverse_effective_mass) * u.y(),
                       pref * std::sqrt(m_longitudinal_effective_mass) * u.z()};
    }

    vector3 k_global_from_energy_direction(double energy, const vector3& direction_valley_frame) const {
        return to_global_frame(k_valley_from_energy_direction(energy, direction_valley_frame));
    }

    /**
     * @brief Non-parabolic density-of-states correction factor often appearing in rates.
     *        For Kane:
     *        g(E) ~ sqrt(E (1 + alpha E)) * (1 + 2 alpha E)
     */
    double density_of_states_jacobian(double energy) const {
        if (energy < 0.0) {
            throw std::invalid_argument("kinetic energy must be non-negative");
        }

        if (m_dispersion == band_type::parabolic || m_non_parabolicity == 0.0) {
            return std::sqrt(energy);
        }

        return std::sqrt(energy * (1.0 + m_non_parabolicity * energy)) * (1.0 + 2.0 * m_non_parabolicity * energy);
    }

    void validate() const {
        if (m_transverse_effective_mass <= 0.0) {
            throw std::invalid_argument("transverse effective mass must be > 0");
        }

        if (m_longitudinal_effective_mass <= 0.0) {
            throw std::invalid_argument("longitudinal effective mass must be > 0");
        }

        if (m_non_parabolicity < 0.0) {
            throw std::invalid_argument("non-parabolicity must be >= 0");
        }

        if (m_degeneracy == 0) {
            throw std::invalid_argument("degeneracy must be >= 1");
        }

        validate_rotation_matrix();
    }

    static mat3 identity_matrix() noexcept { return {{{{1.0, 0.0, 0.0}}, {{0.0, 1.0, 0.0}}, {{0.0, 0.0, 1.0}}}}; }

 private:
    static mat3 transpose(const mat3& m) noexcept {
        return {{{{m[0][0], m[1][0], m[2][0]}}, {{m[0][1], m[1][1], m[2][1]}}, {{m[0][2], m[1][2], m[2][2]}}}};
    }

    static vector3 multiply(const mat3& m, const vector3& v) noexcept {
        return vector3{m[0][0] * v.x() + m[0][1] * v.y() + m[0][2] * v.z(),
                       m[1][0] * v.x() + m[1][1] * v.y() + m[1][2] * v.z(),
                       m[2][0] * v.x() + m[2][1] * v.y() + m[2][2] * v.z()};
    }

    void validate_rotation_matrix() const {
        const mat3 rt                 = transpose(m_rotation);
        const mat3 should_be_identity = multiply(rt, m_rotation);

        constexpr double tol = 1e-10;
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                const double expected = (i == j) ? 1.0 : 0.0;
                if (std::abs(should_be_identity[i][j] - expected) > tol) {
                    throw std::invalid_argument("rotation matrix must be orthonormal");
                }
            }
        }
    }

    static mat3 multiply(const mat3& a, const mat3& b) noexcept {
        mat3 out{};

        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                out[i][j] = 0.0;
                for (std::size_t k = 0; k < 3; ++k) {
                    out[i][j] += a[i][k] * b[k][j];
                }
            }
        }

        return out;
    }
};

}  // namespace uepm::amc