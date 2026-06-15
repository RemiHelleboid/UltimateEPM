/**
 * @file pbmc_avalanche_detector.hpp
 * @brief Voltage-drop based avalanche detection for circuit-coupled PBMC simulations.
 */

#pragma once

#include <cmath>
#include <optional>
#include <stdexcept>

namespace uepm::PBMC {

struct avalanche_detection_state {
    bool                  m_detected = false;
    std::optional<double> m_time_s{};
    std::optional<double> m_voltage_drop_V{};
};

class voltage_drop_avalanche_detector {
 public:
    explicit voltage_drop_avalanche_detector(double threshold_V) : m_threshold_V(threshold_V) {
        if (!std::isfinite(m_threshold_V) || m_threshold_V <= 0.0) {
            throw std::invalid_argument("Avalanche voltage-drop threshold must be positive and finite.");
        }
    }

    void update(bool circuit_enabled, double time_s, double voltage_drop_V) {
        if (m_state.m_detected || !circuit_enabled) {
            return;
        }
        if (!std::isfinite(time_s) || time_s < 0.0) {
            throw std::invalid_argument("Avalanche detection time must be non-negative and finite.");
        }
        if (!std::isfinite(voltage_drop_V)) {
            throw std::invalid_argument("Avalanche voltage drop must be finite.");
        }
        if (std::abs(voltage_drop_V) >= m_threshold_V) {
            m_state.m_detected       = true;
            m_state.m_time_s         = time_s;
            m_state.m_voltage_drop_V = voltage_drop_V;
        }
    }

    [[nodiscard]] double                           threshold_V() const { return m_threshold_V; }
    [[nodiscard]] const avalanche_detection_state& state() const { return m_state; }

 private:
    double                    m_threshold_V;
    avalanche_detection_state m_state{};
};

}  // namespace uepm::PBMC
