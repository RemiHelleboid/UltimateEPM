/**
 * @file amc_quench_detector.hpp
 * @brief Successful-quench detection for AMC device simulations.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace uepm::amc {

struct successful_quench_detection_state {
    bool                  m_detected = false;
    std::optional<double> m_time_s{};
    std::optional<double> m_quiet_period_start_s{};
};

class successful_quench_detector {
 public:
    explicit successful_quench_detector(double quiet_time_s) : m_quiet_time_s(quiet_time_s) {
        if (!std::isfinite(m_quiet_time_s) || m_quiet_time_s <= 0.0) {
            throw std::invalid_argument("Successful-quench quiet time must be positive and finite.");
        }
    }

    void update(bool   avalanche_detected,
                double time_s,
                bool   has_high_field_particle,
                bool   had_impact_ionization_event) {
        if (m_state.m_detected || !avalanche_detected) {
            return;
        }
        if (!std::isfinite(time_s) || time_s < 0.0) {
            throw std::invalid_argument("Successful-quench detection time must be non-negative and finite.");
        }

        if (has_high_field_particle || had_impact_ionization_event) {
            m_state.m_quiet_period_start_s.reset();
            return;
        }

        if (!m_state.m_quiet_period_start_s.has_value()) {
            m_state.m_quiet_period_start_s = time_s;
            return;
        }

        const double elapsed_s = time_s - *m_state.m_quiet_period_start_s;
        const double comparison_scale =
            std::max({std::abs(time_s), std::abs(*m_state.m_quiet_period_start_s), m_quiet_time_s});
        const double tolerance_s = 16.0 * std::numeric_limits<double>::epsilon() * comparison_scale;
        if (elapsed_s + tolerance_s >= m_quiet_time_s) {
            m_state.m_detected = true;
            m_state.m_time_s   = time_s;
        }
    }

    [[nodiscard]] double                                   quiet_time_s() const { return m_quiet_time_s; }
    [[nodiscard]] const successful_quench_detection_state& state() const { return m_state; }

 private:
    double                            m_quiet_time_s;
    successful_quench_detection_state m_state{};
};

}  // namespace uepm::amc
