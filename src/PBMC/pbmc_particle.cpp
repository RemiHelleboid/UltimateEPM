/**
 * @file pbmc_particle.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-04-10
 *
 *
 */

#include "pbmc_particle.hpp"

#include <fmt/core.h>
#include <fmt/format.h>

#include <cmath>
#include <fstream>
#include <iostream>

#include "physical_constants.hpp"
#include "unit_conversion.hpp"

namespace uepm::PBMC {

std::string_view carrier_type_to_string(particle_type type) {
    switch (type) {
        case particle_type::electron:
            return "electron";
        case particle_type::hole:
            return "hole";
        default:
            return "unknown";
    }
}

double signed_charge_C(particle_type type) {
    constexpr double q = uepm::constants::q_e;
    switch (type) {
        case particle_type::electron:
            return -q;
        case particle_type::hole:
            return q;
        default:
            throw std::invalid_argument("invalid particle type");
    }
}

void pbmc_particle::print_info() const {
    std::cout << "Particle index: " << m_index << "\n";
    std::cout << "Type: " << carrier_type_to_string(m_type) << "\n";
    std::cout << "State:\n";
    std::cout << "  Time: " << m_state.time << " s\n";
    std::cout << "  Position: (" << m_state.position.x() << ", " << m_state.position.y() << ", " << m_state.position.z()
              << ") \n";
    std::cout << "  Local k: (" << m_state.local_k.x() << ", " << m_state.local_k.y() << ", " << m_state.local_k.z()
              << ") 1/m\n";
    std::cout << "  Velocity: (" << m_state.velocity.x() << ", " << m_state.velocity.y() << ", " << m_state.velocity.z()
              << ") m/s\n";
    std::cout << "  Kinetic energy: " << m_state.kinetic_energy << " eV\n";
    std::cout << "  Gamma: " << m_state.gamma << " eV\n";
    std::cout << "  Valley index: " << m_state.valley_index << "\n";
}

void pbmc_particle::set_data_from_device(int m_dimension) {
    mesh::vector3 interp_position = m_state.position;
    if (m_dimension == 2) {
        interp_position.to_2d_inplace();
    }
    if (m_state.m_containing_element != nullptr) {
        m_state.electric_field = m_state.m_containing_element->interpolate_electric_field_at_location(interp_position);
        m_state.doping_concentration_cm_3 =
            m_state.m_containing_element->interpolate_doping_at_location(interp_position);
        m_state.impurity_concentration_cm_3 = std::abs(m_state.doping_concentration_cm_3);
        m_state.lattice_temperature_K =
            m_state.m_containing_element->interpolate_temperature_at_location(interp_position);
    } else {
        std::cout << "Error no element at particle position." << std::endl;
    }
}

double pbmc_particle::compute_raw_impact_ionization_coefficient() const {
    // NB II / drift distance
    double drift_distance_m = std::abs(m_state.position.x()) *
                              uepm::units::micron_to_meter;  // assuming drift along x and position in microns
    // fmt::print("Drift distance for impact ionization coefficient: {:.3e} m\n", drift_distance_m);
    double ii_coef = 1e-2 *
                     static_cast<double>(
                         history().scattering_events()[static_cast<std::size_t>(scattering_event::impact_ionization)]) /
                     drift_distance_m;
    return ii_coef;
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

}  // namespace uepm::PBMC
