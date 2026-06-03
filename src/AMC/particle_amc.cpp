/**
 * @file particle_amc.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-04-10
 *
 *
 */

#include "particle_amc.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

namespace uepm::amc {

void particle_amc::print_info() const {
    std::cout << "Particle index: " << m_index << "\n";
    std::cout << "Type: " << (m_type == particle_type::electron ? "Electron" : "Hole") << "\n";
    std::cout << "State:\n";
    std::cout << "  Time: " << m_state.time << " s\n";
    std::cout << "  Position: (" << m_state.position.x() << ", " << m_state.position.y() << ", " << m_state.position.z()
              << ") m\n";
    std::cout << "  Local k: (" << m_state.local_k.x() << ", " << m_state.local_k.y() << ", " << m_state.local_k.z()
              << ") 1/m\n";
    std::cout << "  Velocity: (" << m_state.velocity.x() << ", " << m_state.velocity.y() << ", " << m_state.velocity.z()
              << ") m/s\n";
    std::cout << "  Kinetic energy: " << m_state.kinetic_energy << " eV\n";
    std::cout << "  Gamma: " << m_state.gamma << " eV\n";
    std::cout << "  Valley index: " << m_state.valley_index << "\n";
}

void particle_amc::set_data_from_device(int m_dimension) {
    mesh::vector3 interp_position = m_state.position;
    if (m_dimension == 2) {
        interp_position.to_2d_inplace();
    }
    if (m_state.m_containing_element != nullptr) {
        m_state.electric_field = m_state.m_containing_element->interpolate_electric_field_at_location(interp_position);
        m_state.doping_concentration_cm_3 =
            m_state.m_containing_element->interpolate_doping_at_location(interp_position);
        m_state.impurity_concentration_cm_3 = std::abs(m_state.doping_concentration_cm_3);
    } else {
        std::cout << "Error no element at particle position." << std::endl;
    }
}

double particle_amc::compute_raw_impact_ionization_coefficient() const {
    // NB II / drift distance
    double drift_distance_m = m_state.position.x();
    return 0.01 *
           static_cast<double>(
               history().scattering_events()[static_cast<std::size_t>(scattering_event::impact_ionization)]) /
           drift_distance_m;
}

void particle_history::export_trajectory_as_csv(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file for writing trajectory data.");
    }

    // Write header
    file << "time,X,Y,Z,electric_field,kinetic_energy,gamma,valley_index\n";

    // Write data
    for (const auto& snapshot : m_snapshots) {
        file << snapshot.time << "," << snapshot.position.x() << "," << snapshot.position.y() << ","
             << snapshot.position.z() << "," << snapshot.electric_field_norm << "," << snapshot.kinetic_energy << ","
             << snapshot.gamma << "," << snapshot.valley_index << "\n";
    }

    file.close();
}

}  // namespace uepm::amc
