/**
 * @file amc_quench_circuit.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-06-08
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#pragma once

namespace uepm::amc {

struct passive_quench_circuit_options {
    bool   m_enabled                  = false;
    double m_bias_voltage_V           = 0.0;
    double m_initial_device_voltage_V = 0.0;
    double m_resistance_ohm           = 1.0e5;
    double m_capacitance_F            = 1.0e-15;
};

struct passive_quench_circuit_state {
    double m_time_s             = 0.0;
    double m_device_voltage_V   = 0.0;
    double m_device_current_A   = 0.0;
    double m_resistor_current_A = 0.0;
    double m_voltage_drop_V     = 0.0;
};

class passive_quench_circuit {
 public:
    explicit passive_quench_circuit(passive_quench_circuit_options options = {});

    void set_options(const passive_quench_circuit_options& options);
    void reset();

    [[nodiscard]] bool is_enabled() const;

    [[nodiscard]] const passive_quench_circuit_options& options() const;
    [[nodiscard]] const passive_quench_circuit_state&   state() const;

    [[nodiscard]] double device_voltage_V() const;
    [[nodiscard]] double device_current_A() const;
    [[nodiscard]] double resistor_current_A() const;
    [[nodiscard]] double voltage_drop_V() const;

    passive_quench_circuit_state advance(double device_current_A, double dt_s);

 private:
    passive_quench_circuit_options m_options;
    passive_quench_circuit_state   m_state;

    void validate_options() const;
};

}  // namespace uepm::amc