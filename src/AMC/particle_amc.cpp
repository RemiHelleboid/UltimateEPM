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

#include <iostream>

namespace uepm::amc {

void particle_amc::print_info() const {
    std::cout << "Particle index: " << m_index << "\n";
    std::cout << "Type: " << (m_type == particle_type::electron ? "Electron" : "Hole") << "\n";
    std::cout << "State:\n";
    std::cout << "  Time: " << m_state.time << " s\n";
    std::cout << "  Position: (" << m_state.position.x() << ", " << m_state.position.y() << ", " << m_state.position.z()
              << ") m\n";
    std::cout << "  Local k: (" << m_state.local_k.x() << ", " << m_state.local_k.y() << ", "
              << m_state.local_k.z() << ") 1/m\n";
    std::cout << "  Velocity: (" << m_state.velocity.x() << ", " << m_state.velocity.y() << ", "
              << m_state.velocity.z() << ") m/s\n";
    std::cout << "  Kinetic energy: " << m_state.kinetic_energy << " eV\n";
    std::cout << "  Gamma: " << m_state.gamma << " eV\n";
    std::cout << "  Valley index: " << m_state.valley_index << "\n";
}

}  // namespace uepm::amc
