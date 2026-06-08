/**
 * @file amc_quench_circuit.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-06-08
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#include "amc_quench_circuit.hpp"

#include <cmath>
#include <stdexcept>


namespace uepm::amc {

passive_quench_circuit::passive_quench_circuit(passive_quench_circuit_options options) : m_options(options) {
    validate_options();
    reset();
}

void passive_quench_circuit::set_options(const passive_quench_circuit_options& options) {
    m_options = options;
    validate_options();
    reset();
}

void passive_quench_circuit::reset() {
    m_state                    = {};
    m_state.m_device_voltage_V = m_options.m_initial_device_voltage_V;
    m_state.m_voltage_drop_V   = m_options.m_bias_voltage_V - m_state.m_device_voltage_V;

    if (m_options.m_resistance_ohm > 0.0) {
        m_state.m_resistor_current_A = m_state.m_voltage_drop_V / m_options.m_resistance_ohm;
    }
}

bool passive_quench_circuit::is_enabled() const { return m_options.m_enabled; }

const passive_quench_circuit_options& passive_quench_circuit::options() const { return m_options; }

const passive_quench_circuit_state& passive_quench_circuit::state() const { return m_state; }

double passive_quench_circuit::device_voltage_V() const { return m_state.m_device_voltage_V; }

double passive_quench_circuit::device_current_A() const { return m_state.m_device_current_A; }

double passive_quench_circuit::resistor_current_A() const { return m_state.m_resistor_current_A; }

double passive_quench_circuit::voltage_drop_V() const { return m_state.m_voltage_drop_V; }

passive_quench_circuit_state passive_quench_circuit::advance(double device_current_A, double dt_s) {
    if (!std::isfinite(device_current_A)) {
        throw std::invalid_argument("quench circuit device current must be finite");
    }

    if (!std::isfinite(dt_s) || dt_s < 0.0) {
        throw std::invalid_argument("quench circuit timestep must be finite and non-negative");
    }

    if (!m_options.m_enabled) {
        return m_state;
    }

    const double resistance  = m_options.m_resistance_ohm;
    const double capacitance = m_options.m_capacitance_F;
    const double tau_s       = resistance * capacitance;

    const double steady_voltage_V = m_options.m_bias_voltage_V - resistance * device_current_A;

    if (dt_s > 0.0) {
        const double decay = std::exp(-dt_s / tau_s);

        m_state.m_device_voltage_V = steady_voltage_V + (m_state.m_device_voltage_V - steady_voltage_V) * decay;

        m_state.m_time_s += dt_s;
    }

    m_state.m_device_current_A   = device_current_A;
    m_state.m_voltage_drop_V     = m_options.m_bias_voltage_V - m_state.m_device_voltage_V;
    m_state.m_resistor_current_A = m_state.m_voltage_drop_V / resistance;

    return m_state;
}

void passive_quench_circuit::validate_options() const {
    if (!std::isfinite(m_options.m_bias_voltage_V)) {
        throw std::invalid_argument("quench circuit bias voltage must be finite");
    }

    if (!std::isfinite(m_options.m_initial_device_voltage_V)) {
        throw std::invalid_argument("quench circuit initial device voltage must be finite");
    }

    if (!std::isfinite(m_options.m_resistance_ohm) || m_options.m_resistance_ohm <= 0.0) {
        throw std::invalid_argument("quench circuit resistance must be positive");
    }

    if (!std::isfinite(m_options.m_capacitance_F) || m_options.m_capacitance_F <= 0.0) {
        throw std::invalid_argument("quench circuit capacitance must be positive");
    }
}

}  // namespace uepm::amc