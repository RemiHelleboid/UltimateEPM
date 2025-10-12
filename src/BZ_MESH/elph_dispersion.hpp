/**
 * @file elph_dispersion.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2025-10-10
 * 
 * 
 */

#pragma once

#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "elph_common.hpp"

namespace uepm::mesh_bz {


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

}  // namespace uepm::mesh_bz