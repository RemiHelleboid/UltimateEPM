/**
 * @file amc_device_history.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-05-26
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#pragma once
#include <vector>

#include "vector.hpp"
#include "export_vector_to_csv.hpp"

namespace uepm::amc {

/**
 * @brief Structure contatining the history of the device_simulation states.
 *
 */
struct history_device_amc {
    mesh::vector3              m_last_impact_ionization_position{};
    std::vector<double>        m_list_times{};
    std::vector<std::size_t>   m_list_nb_electrons{};
    std::vector<std::size_t>   m_list_nb_holes{};
    std::vector<std::size_t>   m_list_nb_impact_ionization{};
    std::vector<double>        m_list_anode_current{};
    std::vector<double>        m_list_cathode_current{};
    std::vector<mesh::vector3> m_impact_ionization_positions{};
    std::vector<double>        m_max_electric_field{};
    std::size_t                m_initial_seed_rng{0};

    history_device_amc() = default;

    void reserve_memory(std::size_t size) {
        m_list_times.reserve(size);
        m_list_nb_electrons.reserve(size);
        m_list_nb_holes.reserve(size);
        m_list_nb_impact_ionization.reserve(size);
        m_list_anode_current.reserve(size);
        m_list_cathode_current.reserve(size);
        m_max_electric_field.reserve(size);
    }

    void add_data_to_history(double      time,
                             std::size_t nb_electrons,
                             std::size_t nb_holes,
                             std::size_t nb_impact_ionization,
                             double      anode_current,
                             double      cathode_current,
                             double      max_electric_field) {
        m_list_times.push_back(time);
        m_list_nb_electrons.push_back(nb_electrons);
        m_list_nb_holes.push_back(nb_holes);
        m_list_nb_impact_ionization.push_back(nb_impact_ionization);
        m_list_anode_current.push_back(anode_current);
        m_list_cathode_current.push_back(cathode_current);
    }

    void print_header_csv(const std::string &filename) {
        std::ofstream file(filename);
        file << "time,nb_electrons,nb_holes,nb_impact_ionization,anode_current,cathode_current,max_electric_field\n";
        file.close();
    }

    void append_last_iter_to_csv(std::fstream &file) {
        file << m_list_times.back() << ',' << m_list_nb_electrons.back() << ',' << m_list_nb_holes.back() << ','
             << m_list_nb_impact_ionization.back() << ',' << m_list_anode_current.back() << ',' << m_list_cathode_current.back() << ','
             << m_max_electric_field.back() << '\n';
    }

    void export_to_csv(const std::string &filename, std::size_t frequency = 1) {
        std::vector<double> double_list_time;
        std::vector<double> double_list_nb_electrons;
        std::vector<double> double_list_nb_hole;
        std::vector<double> double_list_nb_impact_ionization;
        std::vector<double> double_list_anode_current;
        std::vector<double> double_list_cathode_current;
        std::vector<double> double_list_max_electric_field;
        for (std::size_t iter_nb = 0; iter_nb < m_list_times.size() - 1; iter_nb++) {
            double_list_time.push_back(m_list_times[iter_nb]);
            double_list_nb_electrons.push_back(m_list_nb_electrons[iter_nb]);
            double_list_nb_hole.push_back(m_list_nb_holes[iter_nb]);
            double_list_nb_impact_ionization.push_back(m_list_nb_impact_ionization[iter_nb]);
            double_list_anode_current.push_back(m_list_anode_current[iter_nb]);
            double_list_cathode_current.push_back(m_list_cathode_current[iter_nb]);
            double_list_max_electric_field.push_back(m_max_electric_field[iter_nb]);
        }
        // Always add the last iteration
        double_list_time.push_back(m_list_times[m_list_nb_electrons.size() - 1]);
        double_list_nb_electrons.push_back(m_list_nb_electrons[m_list_nb_electrons.size() - 1]);
        double_list_nb_hole.push_back(m_list_nb_holes[m_list_nb_electrons.size() - 1]);
        double_list_nb_impact_ionization.push_back(m_list_nb_impact_ionization[m_list_nb_electrons.size() - 1]);
        double_list_anode_current.push_back(m_list_anode_current[m_list_nb_electrons.size() - 1]);
        double_list_cathode_current.push_back(m_list_cathode_current[m_list_nb_electrons.size() - 1]);
        double_list_max_electric_field.push_back(m_max_electric_field[m_list_nb_electrons.size() - 1]);

        std::vector<std::string> header_csv =
            {"time", "nb_electrons", "nb_holes", "nb_impact_ionization", "anode_current", "cathode_current", "max_electric_field"};
        utils::export_multiple_vector_to_csv(filename,
                                             header_csv,
                                             {double_list_time,
                                              double_list_nb_electrons,
                                              double_list_nb_hole,
                                              double_list_nb_impact_ionization,
                                              double_list_anode_current,
                                              double_list_cathode_current,
                                              double_list_max_electric_field});
    }

    std::vector<std::size_t> get_history_total_nb_particles() const {
        std::vector<std::size_t> history_total_number_particles(m_list_nb_electrons.size());
        for (std::size_t index_iteration = 0; index_iteration < m_list_nb_electrons.size(); ++index_iteration) {
            history_total_number_particles[index_iteration] = m_list_nb_electrons[index_iteration] + m_list_nb_holes[index_iteration];
        }
        return history_total_number_particles;
    }

};

}