/**
 * @file pbmc_device_history.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-26
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once
#include <algorithm>
#include <cmath>
#include <deque>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "export_vector_to_csv.hpp"
#include "vector.hpp"

namespace uepm::PBMC {

/**
 * @brief Structure contatining the history of the device_simulation states.
 *
 */
struct history_device_PBMC {
    mesh::vector3                    m_last_impact_ionization_position{};
    std::vector<double>              m_list_times{};
    std::vector<std::size_t>         m_list_nb_electrons{};
    std::vector<std::size_t>         m_list_nb_holes{};
    std::vector<double>              m_list_mean_electron_kinetic_energy_eV{};
    std::vector<double>              m_list_mean_hole_kinetic_energy_eV{};
    std::vector<double>              m_list_mean_particle_kinetic_energy_eV{};
    std::vector<std::size_t>         m_list_nb_impact_ionization{};
    std::vector<double>              m_list_ramo_current_electron{};
    std::vector<double>              m_list_ramo_current_hole{};
    std::vector<double>              m_list_ramo_current{};
    std::vector<double>              m_list_probe_ramo_current_electron{};
    std::vector<double>              m_list_probe_ramo_current_hole{};
    std::vector<double>              m_list_probe_ramo_current{};
    std::vector<mesh::vector3>       m_impact_ionization_positions{};
    std::vector<double>              m_max_electric_field{};
    std::vector<double>              m_list_ramo_electrode_voltage_V{};
    std::vector<double>              m_list_reference_electrode_voltage_V{};
    std::vector<double>              m_list_quench_bias_voltage_V{};
    std::vector<double>              m_list_quench_device_current_A{};
    std::vector<double>              m_list_quench_resistor_current_A{};
    std::vector<double>              m_list_quench_voltage_drop_V{};
    std::vector<std::string>         m_contact_voltage_names{};
    std::vector<std::vector<double>> m_list_contact_voltages_V{};
    std::vector<std::string>         m_contact_flow_names{};
    std::vector<std::vector<double>> m_list_collected_current_electron_A{};
    std::vector<std::vector<double>> m_list_collected_current_hole_A{};
    std::vector<std::vector<double>> m_list_collected_current_A{};
    std::vector<std::vector<double>> m_list_cumulative_collected_charge_C{};
    std::vector<std::vector<double>> m_list_injected_current_A{};
    std::vector<std::vector<double>> m_list_net_contact_current_A{};
    std::vector<std::vector<double>> m_list_cumulative_injected_charge_C{};
    std::vector<std::vector<double>> m_list_cumulative_net_contact_charge_C{};
    std::size_t                      m_initial_seed_rng{0};

    struct contact_current_window_sample {
        double                      start_time_s = 0.0;
        double                      end_time_s   = 0.0;
        std::vector<double>         collected_electron_charge_C{};
        std::vector<double>         collected_hole_charge_C{};
        std::vector<double>         injected_charge_C{};
    };

    struct windowed_contact_currents {
        std::vector<double> collected_electron_A{};
        std::vector<double> collected_hole_A{};
        std::vector<double> collected_total_A{};
        std::vector<double> injected_A{};
        std::vector<double> net_A{};
    };

    double                                    m_contact_current_window_s = 0.0;
    std::deque<contact_current_window_sample> m_contact_current_window_samples{};
    std::vector<double>                       m_window_collected_electron_charge_C{};
    std::vector<double>                       m_window_collected_hole_charge_C{};
    std::vector<double>                       m_window_injected_charge_C{};
    double                                    m_contact_current_origin_time_s = 0.0;
    bool                                      m_has_contact_current_origin    = false;

    history_device_PBMC() = default;

    void reserve_memory(std::size_t size) {
        m_list_times.reserve(size);
        m_list_nb_electrons.reserve(size);
        m_list_nb_holes.reserve(size);
        m_list_mean_electron_kinetic_energy_eV.reserve(size);
        m_list_mean_hole_kinetic_energy_eV.reserve(size);
        m_list_mean_particle_kinetic_energy_eV.reserve(size);
        m_list_nb_impact_ionization.reserve(size);
        m_list_ramo_current_electron.reserve(size);
        m_list_ramo_current_hole.reserve(size);
        m_list_ramo_current.reserve(size);
        m_list_probe_ramo_current_electron.reserve(size);
        m_list_probe_ramo_current_hole.reserve(size);
        m_list_probe_ramo_current.reserve(size);
        m_max_electric_field.reserve(size);
        m_list_ramo_electrode_voltage_V.reserve(size);
        m_list_reference_electrode_voltage_V.reserve(size);
        m_list_quench_bias_voltage_V.reserve(size);
        m_list_quench_device_current_A.reserve(size);
        m_list_quench_resistor_current_A.reserve(size);
        m_list_quench_voltage_drop_V.reserve(size);
        m_list_contact_voltages_V.reserve(size);
        m_list_collected_current_electron_A.reserve(size);
        m_list_collected_current_hole_A.reserve(size);
        m_list_collected_current_A.reserve(size);
        m_list_cumulative_collected_charge_C.reserve(size);
        m_list_injected_current_A.reserve(size);
        m_list_net_contact_current_A.reserve(size);
        m_list_cumulative_injected_charge_C.reserve(size);
        m_list_cumulative_net_contact_charge_C.reserve(size);
    }

    void set_contact_voltage_names(const std::vector<std::string>& contact_names) {
        m_contact_voltage_names = contact_names;
    }

    void reset_contact_current_window() {
        m_contact_current_window_samples.clear();
        m_window_collected_electron_charge_C.assign(m_contact_flow_names.size(), 0.0);
        m_window_collected_hole_charge_C.assign(m_contact_flow_names.size(), 0.0);
        m_window_injected_charge_C.assign(m_contact_flow_names.size(), 0.0);
        m_contact_current_origin_time_s = 0.0;
        m_has_contact_current_origin    = false;
    }

    void set_contact_flow_names(const std::vector<std::string>& contact_names) {
        m_contact_flow_names = contact_names;
        reset_contact_current_window();
    }

    void set_contact_current_window_s(double window_s) {
        if (!std::isfinite(window_s) || window_s < 0.0) {
            throw std::invalid_argument("Contact-current averaging window must be finite and non-negative.");
        }
        m_contact_current_window_s = window_s;
        reset_contact_current_window();
    }

    windowed_contact_currents average_contact_currents(
        double                     time_s,
        const std::vector<double>& collected_current_electron_A,
        const std::vector<double>& collected_current_hole_A,
        const std::vector<double>& injected_current_A) {
        const std::size_t contact_count = m_contact_flow_names.size();
        if (collected_current_electron_A.size() != contact_count ||
            collected_current_hole_A.size() != contact_count ||
            injected_current_A.size() != contact_count) {
            throw std::invalid_argument("Contact-current vector size does not match the configured contact count.");
        }

        const double previous_time_s = m_list_times.empty() ? 0.0 : m_list_times.back();
        const double step_duration_s = time_s - previous_time_s;
        if (!(step_duration_s > 0.0) || !std::isfinite(step_duration_s)) {
            throw std::invalid_argument("Contact-current history times must be finite and strictly increasing.");
        }
        if (!m_has_contact_current_origin) {
            m_contact_current_origin_time_s = previous_time_s;
            m_has_contact_current_origin    = true;
        }

        contact_current_window_sample sample;
        sample.start_time_s = previous_time_s;
        sample.end_time_s   = time_s;
        sample.collected_electron_charge_C.resize(contact_count);
        sample.collected_hole_charge_C.resize(contact_count);
        sample.injected_charge_C.resize(contact_count);
        for (std::size_t index = 0; index < contact_count; ++index) {
            sample.collected_electron_charge_C[index] = collected_current_electron_A[index] * step_duration_s;
            sample.collected_hole_charge_C[index]     = collected_current_hole_A[index] * step_duration_s;
            sample.injected_charge_C[index]           = injected_current_A[index] * step_duration_s;
            m_window_collected_electron_charge_C[index] += sample.collected_electron_charge_C[index];
            m_window_collected_hole_charge_C[index] += sample.collected_hole_charge_C[index];
            m_window_injected_charge_C[index] += sample.injected_charge_C[index];
        }
        m_contact_current_window_samples.push_back(std::move(sample));

        const double window_start_s =
            std::max(m_contact_current_origin_time_s, time_s - m_contact_current_window_s);
        while (!m_contact_current_window_samples.empty() &&
               m_contact_current_window_samples.front().end_time_s <= window_start_s) {
            const auto& expired = m_contact_current_window_samples.front();
            for (std::size_t index = 0; index < contact_count; ++index) {
                m_window_collected_electron_charge_C[index] -= expired.collected_electron_charge_C[index];
                m_window_collected_hole_charge_C[index] -= expired.collected_hole_charge_C[index];
                m_window_injected_charge_C[index] -= expired.injected_charge_C[index];
            }
            m_contact_current_window_samples.pop_front();
        }

        auto effective_electron_charge_C = m_window_collected_electron_charge_C;
        auto effective_hole_charge_C     = m_window_collected_hole_charge_C;
        auto effective_injected_charge_C = m_window_injected_charge_C;
        if (!m_contact_current_window_samples.empty()) {
            const auto& oldest = m_contact_current_window_samples.front();
            if (oldest.start_time_s < window_start_s) {
                const double excluded_fraction =
                    (window_start_s - oldest.start_time_s) / (oldest.end_time_s - oldest.start_time_s);
                for (std::size_t index = 0; index < contact_count; ++index) {
                    effective_electron_charge_C[index] -=
                        excluded_fraction * oldest.collected_electron_charge_C[index];
                    effective_hole_charge_C[index] -= excluded_fraction * oldest.collected_hole_charge_C[index];
                    effective_injected_charge_C[index] -= excluded_fraction * oldest.injected_charge_C[index];
                }
            }
        }

        const double averaging_duration_s = time_s - window_start_s;
        windowed_contact_currents result;
        result.collected_electron_A.resize(contact_count);
        result.collected_hole_A.resize(contact_count);
        result.collected_total_A.resize(contact_count);
        result.injected_A.resize(contact_count);
        result.net_A.resize(contact_count);
        for (std::size_t index = 0; index < contact_count; ++index) {
            result.collected_electron_A[index] = effective_electron_charge_C[index] / averaging_duration_s;
            result.collected_hole_A[index]     = effective_hole_charge_C[index] / averaging_duration_s;
            result.collected_total_A[index] =
                result.collected_electron_A[index] + result.collected_hole_A[index];
            result.injected_A[index] = effective_injected_charge_C[index] / averaging_duration_s;
            result.net_A[index]      = result.collected_total_A[index] - result.injected_A[index];
        }
        return result;
    }

    std::vector<double> contact_voltage_values_from_map(const std::map<std::string, double>& contact_voltages_V) const {
        std::vector<double> values;
        values.reserve(m_contact_voltage_names.size());
        for (const auto& contact_name : m_contact_voltage_names) {
            const auto voltage_it = contact_voltages_V.find(contact_name);
            values.push_back(voltage_it == contact_voltages_V.end() ? 0.0 : voltage_it->second);
        }
        return values;
    }

    void add_data_to_history(double                     time,
                             std::size_t                nb_electrons,
                             std::size_t                nb_holes,
                             double                     mean_electron_kinetic_energy_eV,
                             double                     mean_hole_kinetic_energy_eV,
                             double                     mean_particle_kinetic_energy_eV,
                             std::size_t                nb_impact_ionization,
                             double                     ramo_current_electron,
                             double                     ramo_current_hole,
                             double                     ramo_current,
                             double                     probe_ramo_current_electron,
                             double                     probe_ramo_current_hole,
                             double                     probe_ramo_current,
                             double                     max_electric_field,
                             double                     ramo_electrode_voltage_V,
                             double                     reference_electrode_voltage_V,
                             double                     quench_bias_voltage_V,
                             double                     quench_device_current_A,
                             double                     quench_resistor_current_A,
                             double                     quench_voltage_drop_V,
                             const std::vector<double>& contact_voltages_V            = {},
                             const std::vector<double>& collected_current_electron_A  = {},
                             const std::vector<double>& collected_current_hole_A      = {},
                             const std::vector<double>& collected_current_A           = {},
                             const std::vector<double>& cumulative_collected_charge_C = {},
                             const std::vector<double>& injected_current_A = {},
                             const std::vector<double>& net_contact_current_A = {},
                             const std::vector<double>& cumulative_injected_charge_C = {},
                             const std::vector<double>& cumulative_net_contact_charge_C = {}) {
        const bool average_contact_flow = m_contact_current_window_s > 0.0 && !m_contact_flow_names.empty();
        windowed_contact_currents averaged_contact_flow;
        if (average_contact_flow) {
            averaged_contact_flow = average_contact_currents(
                time, collected_current_electron_A, collected_current_hole_A, injected_current_A);
        }

        m_list_times.push_back(time);
        m_list_nb_electrons.push_back(nb_electrons);
        m_list_nb_holes.push_back(nb_holes);
        m_list_mean_electron_kinetic_energy_eV.push_back(mean_electron_kinetic_energy_eV);
        m_list_mean_hole_kinetic_energy_eV.push_back(mean_hole_kinetic_energy_eV);
        m_list_mean_particle_kinetic_energy_eV.push_back(mean_particle_kinetic_energy_eV);
        m_list_nb_impact_ionization.push_back(nb_impact_ionization);
        m_list_ramo_current_electron.push_back(ramo_current_electron);
        m_list_ramo_current_hole.push_back(ramo_current_hole);
        m_list_ramo_current.push_back(ramo_current);
        m_list_probe_ramo_current_electron.push_back(probe_ramo_current_electron);
        m_list_probe_ramo_current_hole.push_back(probe_ramo_current_hole);
        m_list_probe_ramo_current.push_back(probe_ramo_current);
        m_max_electric_field.push_back(max_electric_field);
        m_list_ramo_electrode_voltage_V.push_back(ramo_electrode_voltage_V);
        m_list_reference_electrode_voltage_V.push_back(reference_electrode_voltage_V);
        m_list_quench_bias_voltage_V.push_back(quench_bias_voltage_V);
        m_list_quench_device_current_A.push_back(quench_device_current_A);
        m_list_quench_resistor_current_A.push_back(quench_resistor_current_A);
        m_list_quench_voltage_drop_V.push_back(quench_voltage_drop_V);
        m_list_contact_voltages_V.push_back(contact_voltages_V);
        m_list_collected_current_electron_A.push_back(
            average_contact_flow ? averaged_contact_flow.collected_electron_A : collected_current_electron_A);
        m_list_collected_current_hole_A.push_back(
            average_contact_flow ? averaged_contact_flow.collected_hole_A : collected_current_hole_A);
        m_list_collected_current_A.push_back(
            average_contact_flow ? averaged_contact_flow.collected_total_A : collected_current_A);
        m_list_cumulative_collected_charge_C.push_back(cumulative_collected_charge_C);
        m_list_injected_current_A.push_back(
            average_contact_flow ? averaged_contact_flow.injected_A : injected_current_A);
        m_list_net_contact_current_A.push_back(
            average_contact_flow ? averaged_contact_flow.net_A : net_contact_current_A);
        m_list_cumulative_injected_charge_C.push_back(cumulative_injected_charge_C);
        m_list_cumulative_net_contact_charge_C.push_back(cumulative_net_contact_charge_C);
    }

    void set_last_contact_voltages(const std::vector<double>& contact_voltages_V) {
        if (!m_list_contact_voltages_V.empty()) {
            m_list_contact_voltages_V.back() = contact_voltages_V;
        }
    }

    void print_header_csv(const std::string& filename) {
        std::ofstream file(filename);
        file << "time,nb_electrons,nb_holes,mean_electron_kinetic_energy_eV,mean_hole_kinetic_energy_eV,"
                "mean_particle_kinetic_energy_eV,nb_impact_ionization,ramo_current_electron,ramo_current_hole,ramo_"
                "current,"
                "probe_ramo_current_electron,probe_ramo_current_hole,probe_ramo_current,"
                "max_electric_field,ramo_electrode_voltage_V,reference_electrode_voltage_V,quench_bias_voltage_V,"
                "quench_device_current_A,"
                "quench_resistor_current_A,quench_voltage_drop_V";
        for (const auto& contact_name : m_contact_voltage_names) {
            file << ",V_" << contact_name;
        }
        for (const auto& contact_name : m_contact_flow_names) {
            file << ",collected_current_electron_" << contact_name << "_A"
                 << ",collected_current_hole_" << contact_name << "_A"
                 << ",collected_current_" << contact_name << "_A"
                 << ",cumulative_collected_charge_" << contact_name << "_C"
                 << ",injected_current_" << contact_name << "_A"
                 << ",net_contact_current_" << contact_name << "_A"
                 << ",cumulative_injected_charge_" << contact_name << "_C"
                 << ",cumulative_net_contact_charge_" << contact_name << "_C";
        }
        file << '\n';
        file.close();
    }

    void append_last_iter_to_csv(std::fstream& file) {
        file << m_list_times.back() << ',' << m_list_nb_electrons.back() << ',' << m_list_nb_holes.back() << ','
             << m_list_mean_electron_kinetic_energy_eV.back() << ',' << m_list_mean_hole_kinetic_energy_eV.back() << ','
             << m_list_mean_particle_kinetic_energy_eV.back() << ',' << m_list_nb_impact_ionization.back() << ','
             << m_list_ramo_current_electron.back() << ',' << m_list_ramo_current_hole.back() << ','
             << m_list_ramo_current.back() << ',' << m_list_probe_ramo_current_electron.back() << ','
             << m_list_probe_ramo_current_hole.back() << ',' << m_list_probe_ramo_current.back() << ','
             << m_max_electric_field.back() << ',' << m_list_ramo_electrode_voltage_V.back() << ','
             << m_list_reference_electrode_voltage_V.back() << ',' << m_list_quench_bias_voltage_V.back() << ','
             << m_list_quench_device_current_A.back() << ',' << m_list_quench_resistor_current_A.back() << ','
             << m_list_quench_voltage_drop_V.back();
        const auto& contact_voltages_V =
            m_list_contact_voltages_V.empty() ? std::vector<double>{} : m_list_contact_voltages_V.back();
        for (std::size_t index = 0; index < m_contact_voltage_names.size(); ++index) {
            file << ',' << (index < contact_voltages_V.size() ? contact_voltages_V[index] : 0.0);
        }
        const auto append_contact_flow_values = [&](const std::vector<std::vector<double>>& rows,
                                                    std::size_t                             contact_index) {
            if (rows.empty() || contact_index >= rows.back().size()) {
                file << ",0";
            } else {
                file << ',' << rows.back()[contact_index];
            }
        };
        for (std::size_t index = 0; index < m_contact_flow_names.size(); ++index) {
            append_contact_flow_values(m_list_collected_current_electron_A, index);
            append_contact_flow_values(m_list_collected_current_hole_A, index);
            append_contact_flow_values(m_list_collected_current_A, index);
            append_contact_flow_values(m_list_cumulative_collected_charge_C, index);
            append_contact_flow_values(m_list_injected_current_A, index);
            append_contact_flow_values(m_list_net_contact_current_A, index);
            append_contact_flow_values(m_list_cumulative_injected_charge_C, index);
            append_contact_flow_values(m_list_cumulative_net_contact_charge_C, index);
        }
        file << '\n';
    }

    void export_to_csv(const std::string& filename, std::size_t frequency = 1) {
        if (frequency == 0) {
            throw std::invalid_argument("history export frequency must be greater than zero");
        }

        if (m_list_times.empty()) {
            print_header_csv(filename);
            return;
        }

        std::vector<double>              double_list_time;
        std::vector<double>              double_list_nb_electrons;
        std::vector<double>              double_list_nb_hole;
        std::vector<double>              double_list_mean_electron_kinetic_energy_eV;
        std::vector<double>              double_list_mean_hole_kinetic_energy_eV;
        std::vector<double>              double_list_mean_particle_kinetic_energy_eV;
        std::vector<double>              double_list_nb_impact_ionization;
        std::vector<double>              double_list_ramo_current_electron;
        std::vector<double>              double_list_ramo_current_hole;
        std::vector<double>              double_list_ramo_current;
        std::vector<double>              double_list_probe_ramo_current_electron;
        std::vector<double>              double_list_probe_ramo_current_hole;
        std::vector<double>              double_list_probe_ramo_current;
        std::vector<double>              double_list_ramo_electrode_voltage_V;
        std::vector<double>              double_list_reference_electrode_voltage_V;
        std::vector<double>              double_list_quench_bias_voltage_V;
        std::vector<double>              double_list_quench_device_current_A;
        std::vector<double>              double_list_quench_resistor_current_A;
        std::vector<double>              double_list_quench_voltage_drop_V;
        std::vector<std::vector<double>> double_list_contact_voltages_V(m_contact_voltage_names.size());
        std::vector<std::vector<double>> double_list_collected_current_electron_A(m_contact_flow_names.size());
        std::vector<std::vector<double>> double_list_collected_current_hole_A(m_contact_flow_names.size());
        std::vector<std::vector<double>> double_list_collected_current_A(m_contact_flow_names.size());
        std::vector<std::vector<double>> double_list_cumulative_collected_charge_C(m_contact_flow_names.size());
        std::vector<std::vector<double>> double_list_injected_current_A(m_contact_flow_names.size());
        std::vector<std::vector<double>> double_list_net_contact_current_A(m_contact_flow_names.size());
        std::vector<std::vector<double>> double_list_cumulative_injected_charge_C(m_contact_flow_names.size());
        std::vector<std::vector<double>> double_list_cumulative_net_contact_charge_C(m_contact_flow_names.size());

        const auto append_contact_flow_row = [&](std::size_t row_index) {
            for (std::size_t contact_index = 0; contact_index < m_contact_flow_names.size(); ++contact_index) {
                const auto value_at = [&](const std::vector<std::vector<double>>& rows) {
                    return row_index < rows.size() && contact_index < rows[row_index].size()
                               ? rows[row_index][contact_index]
                               : 0.0;
                };
                double_list_collected_current_electron_A[contact_index].push_back(
                    value_at(m_list_collected_current_electron_A));
                double_list_collected_current_hole_A[contact_index].push_back(
                    value_at(m_list_collected_current_hole_A));
                double_list_collected_current_A[contact_index].push_back(value_at(m_list_collected_current_A));
                double_list_cumulative_collected_charge_C[contact_index].push_back(
                    value_at(m_list_cumulative_collected_charge_C));
                double_list_injected_current_A[contact_index].push_back(value_at(m_list_injected_current_A));
                double_list_net_contact_current_A[contact_index].push_back(value_at(m_list_net_contact_current_A));
                double_list_cumulative_injected_charge_C[contact_index].push_back(
                    value_at(m_list_cumulative_injected_charge_C));
                double_list_cumulative_net_contact_charge_C[contact_index].push_back(
                    value_at(m_list_cumulative_net_contact_charge_C));
            }
        };

        std::vector<double> double_list_max_electric_field;
        for (std::size_t iter_nb = 0; iter_nb < m_list_times.size() - 1; iter_nb += frequency) {
            double_list_time.push_back(m_list_times[iter_nb]);
            double_list_nb_electrons.push_back(m_list_nb_electrons[iter_nb]);
            double_list_nb_hole.push_back(m_list_nb_holes[iter_nb]);
            double_list_mean_electron_kinetic_energy_eV.push_back(m_list_mean_electron_kinetic_energy_eV[iter_nb]);
            double_list_mean_hole_kinetic_energy_eV.push_back(m_list_mean_hole_kinetic_energy_eV[iter_nb]);
            double_list_mean_particle_kinetic_energy_eV.push_back(m_list_mean_particle_kinetic_energy_eV[iter_nb]);
            double_list_nb_impact_ionization.push_back(m_list_nb_impact_ionization[iter_nb]);
            double_list_ramo_current_electron.push_back(m_list_ramo_current_electron[iter_nb]);
            double_list_ramo_current_hole.push_back(m_list_ramo_current_hole[iter_nb]);
            double_list_ramo_current.push_back(m_list_ramo_current[iter_nb]);
            double_list_probe_ramo_current_electron.push_back(m_list_probe_ramo_current_electron[iter_nb]);
            double_list_probe_ramo_current_hole.push_back(m_list_probe_ramo_current_hole[iter_nb]);
            double_list_probe_ramo_current.push_back(m_list_probe_ramo_current[iter_nb]);
            double_list_max_electric_field.push_back(m_max_electric_field[iter_nb]);
            double_list_ramo_electrode_voltage_V.push_back(m_list_ramo_electrode_voltage_V[iter_nb]);
            double_list_reference_electrode_voltage_V.push_back(m_list_reference_electrode_voltage_V[iter_nb]);
            double_list_quench_bias_voltage_V.push_back(m_list_quench_bias_voltage_V[iter_nb]);
            double_list_quench_device_current_A.push_back(m_list_quench_device_current_A[iter_nb]);
            double_list_quench_resistor_current_A.push_back(m_list_quench_resistor_current_A[iter_nb]);
            double_list_quench_voltage_drop_V.push_back(m_list_quench_voltage_drop_V[iter_nb]);
            for (std::size_t contact_index = 0; contact_index < m_contact_voltage_names.size(); ++contact_index) {
                const auto& contact_voltages_V = m_list_contact_voltages_V[iter_nb];
                double_list_contact_voltages_V[contact_index].push_back(
                    contact_index < contact_voltages_V.size() ? contact_voltages_V[contact_index] : 0.0);
            }
            append_contact_flow_row(iter_nb);
        }
        // Always add the last iteration
        double_list_time.push_back(m_list_times[m_list_nb_electrons.size() - 1]);
        double_list_nb_electrons.push_back(m_list_nb_electrons[m_list_nb_electrons.size() - 1]);
        double_list_nb_hole.push_back(m_list_nb_holes[m_list_nb_electrons.size() - 1]);
        double_list_mean_electron_kinetic_energy_eV.push_back(
            m_list_mean_electron_kinetic_energy_eV[m_list_nb_electrons.size() - 1]);
        double_list_mean_hole_kinetic_energy_eV.push_back(
            m_list_mean_hole_kinetic_energy_eV[m_list_nb_electrons.size() - 1]);
        double_list_mean_particle_kinetic_energy_eV.push_back(
            m_list_mean_particle_kinetic_energy_eV[m_list_nb_electrons.size() - 1]);
        double_list_nb_impact_ionization.push_back(m_list_nb_impact_ionization[m_list_nb_electrons.size() - 1]);
        double_list_ramo_current_electron.push_back(m_list_ramo_current_electron[m_list_nb_electrons.size() - 1]);
        double_list_ramo_current_hole.push_back(m_list_ramo_current_hole[m_list_nb_electrons.size() - 1]);
        double_list_ramo_current.push_back(m_list_ramo_current[m_list_nb_electrons.size() - 1]);
        double_list_probe_ramo_current_electron.push_back(
            m_list_probe_ramo_current_electron[m_list_nb_electrons.size() - 1]);
        double_list_probe_ramo_current_hole.push_back(m_list_probe_ramo_current_hole[m_list_nb_electrons.size() - 1]);
        double_list_probe_ramo_current.push_back(m_list_probe_ramo_current[m_list_nb_electrons.size() - 1]);
        double_list_max_electric_field.push_back(m_max_electric_field[m_list_nb_electrons.size() - 1]);
        double_list_ramo_electrode_voltage_V.push_back(m_list_ramo_electrode_voltage_V[m_list_nb_electrons.size() - 1]);
        double_list_reference_electrode_voltage_V.push_back(
            m_list_reference_electrode_voltage_V[m_list_nb_electrons.size() - 1]);
        double_list_quench_bias_voltage_V.push_back(m_list_quench_bias_voltage_V[m_list_nb_electrons.size() - 1]);
        double_list_quench_device_current_A.push_back(m_list_quench_device_current_A[m_list_nb_electrons.size() - 1]);
        double_list_quench_resistor_current_A.push_back(
            m_list_quench_resistor_current_A[m_list_nb_electrons.size() - 1]);
        double_list_quench_voltage_drop_V.push_back(m_list_quench_voltage_drop_V[m_list_nb_electrons.size() - 1]);
        for (std::size_t contact_index = 0; contact_index < m_contact_voltage_names.size(); ++contact_index) {
            const auto& contact_voltages_V = m_list_contact_voltages_V[m_list_nb_electrons.size() - 1];
            double_list_contact_voltages_V[contact_index].push_back(
                contact_index < contact_voltages_V.size() ? contact_voltages_V[contact_index] : 0.0);
        }
        append_contact_flow_row(m_list_nb_electrons.size() - 1);

        std::vector<std::string> header_csv = {"time",
                                               "nb_electrons",
                                               "nb_holes",
                                               "mean_electron_kinetic_energy_eV",
                                               "mean_hole_kinetic_energy_eV",
                                               "mean_particle_kinetic_energy_eV",
                                               "nb_impact_ionization",
                                               "ramo_current_electron",
                                               "ramo_current_hole",
                                               "ramo_current",
                                               "probe_ramo_current_electron",
                                               "probe_ramo_current_hole",
                                               "probe_ramo_current",
                                               "max_electric_field",
                                               "ramo_electrode_voltage_V",
                                               "reference_electrode_voltage_V",
                                               "quench_bias_voltage_V",
                                               "quench_device_current_A",
                                               "quench_resistor_current_A",
                                               "quench_voltage_drop_V"};
        for (const auto& contact_name : m_contact_voltage_names) {
            header_csv.push_back("V_" + contact_name);
        }
        for (const auto& contact_name : m_contact_flow_names) {
            header_csv.push_back("collected_current_electron_" + contact_name + "_A");
            header_csv.push_back("collected_current_hole_" + contact_name + "_A");
            header_csv.push_back("collected_current_" + contact_name + "_A");
            header_csv.push_back("cumulative_collected_charge_" + contact_name + "_C");
            header_csv.push_back("injected_current_" + contact_name + "_A");
            header_csv.push_back("net_contact_current_" + contact_name + "_A");
            header_csv.push_back("cumulative_injected_charge_" + contact_name + "_C");
            header_csv.push_back("cumulative_net_contact_charge_" + contact_name + "_C");
        }
        std::vector<std::vector<double>> columns = {double_list_time,
                                                    double_list_nb_electrons,
                                                    double_list_nb_hole,
                                                    double_list_mean_electron_kinetic_energy_eV,
                                                    double_list_mean_hole_kinetic_energy_eV,
                                                    double_list_mean_particle_kinetic_energy_eV,
                                                    double_list_nb_impact_ionization,
                                                    double_list_ramo_current_electron,
                                                    double_list_ramo_current_hole,
                                                    double_list_ramo_current,
                                                    double_list_probe_ramo_current_electron,
                                                    double_list_probe_ramo_current_hole,
                                                    double_list_probe_ramo_current,
                                                    double_list_max_electric_field,
                                                    double_list_ramo_electrode_voltage_V,
                                                    double_list_reference_electrode_voltage_V,
                                                    double_list_quench_bias_voltage_V,
                                                    double_list_quench_device_current_A,
                                                    double_list_quench_resistor_current_A,
                                                    double_list_quench_voltage_drop_V};
        for (const auto& contact_column : double_list_contact_voltages_V) {
            columns.push_back(contact_column);
        }
        for (std::size_t contact_index = 0; contact_index < m_contact_flow_names.size(); ++contact_index) {
            columns.push_back(double_list_collected_current_electron_A[contact_index]);
            columns.push_back(double_list_collected_current_hole_A[contact_index]);
            columns.push_back(double_list_collected_current_A[contact_index]);
            columns.push_back(double_list_cumulative_collected_charge_C[contact_index]);
            columns.push_back(double_list_injected_current_A[contact_index]);
            columns.push_back(double_list_net_contact_current_A[contact_index]);
            columns.push_back(double_list_cumulative_injected_charge_C[contact_index]);
            columns.push_back(double_list_cumulative_net_contact_charge_C[contact_index]);
        }
        utils::export_multiple_vector_to_csv(filename, header_csv, columns);
    }

    std::vector<std::size_t> get_history_total_nb_particles() const {
        std::vector<std::size_t> history_total_number_particles(m_list_nb_electrons.size());
        for (std::size_t index_iteration = 0; index_iteration < m_list_nb_electrons.size(); ++index_iteration) {
            history_total_number_particles[index_iteration] =
                m_list_nb_electrons[index_iteration] + m_list_nb_holes[index_iteration];
        }
        return history_total_number_particles;
    }

    double extract_final_current(double time_window_s) const {
        if (m_list_times.empty()) {
            return 0.0;
        }

        const double final_time_s  = m_list_times.back();
        double       current_sum_A = 0.0;
        std::size_t  count         = 0;

        for (std::size_t index_iteration = m_list_times.size(); index_iteration-- > 0;) {
            if (final_time_s - m_list_times[index_iteration] <= time_window_s) {
                current_sum_A += m_list_ramo_current[index_iteration];
                ++count;
            } else {
                break;
            }
        }

        return count > 0 ? current_sum_A / static_cast<double>(count) : 0.0;
    }
};

}  // namespace uepm::PBMC
