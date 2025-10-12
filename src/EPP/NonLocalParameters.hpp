/**
 * @file NonLocalParameters.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-07-26
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <iostream>
#include <string>
#include <system_error>

#include "Constants.hpp"
#include "yaml-cpp/yaml.h"

namespace uepm::pseudopotential {

enum class non_local_well_type { square, gaussian, unknown };

/**
 * @brief Struct to store the non-local parameters.
 *
 * The formalism is taken from: 1. Pötz, W. & Vogl, P. Theory of optical-phonon deformation potentials in tetrahedral semiconductors.
 *  Phys. Rev. B 24, 2025–2037 (1981).
 *
 */
struct NonLocalParameters {
    double              m_alpha_0_cation;
    double              m_beta_0_cation;
    double              m_A2_cation;
    double              m_alpha_0_anion;
    double              m_beta_0_anion;
    double              m_A2_anion;
    double              m_R0_cation;
    double              m_R2_cation;
    double              m_R0_anion;
    double              m_R2_anion;
    non_local_well_type m_well_type;

    NonLocalParameters()
        : m_alpha_0_cation(0.0),
          m_beta_0_cation(0.0),
          m_A2_cation(0.0),
          m_alpha_0_anion(0.0),
          m_beta_0_anion(0.0),
          m_A2_anion(0.0),
          m_R0_cation(0.0),
          m_R2_cation(0.0),
          m_R0_anion(0.0),
          m_R2_anion(0.0),
          m_well_type(non_local_well_type::unknown) {}

    NonLocalParameters(double              alpha_0_cation,
                       double              beta_0_cation,
                       double              A2_cation,
                       double              alpha_0_anion,
                       double              beta_0_anion,
                       double              A2_anion,
                       double              R0_cation,
                       double              R2_cation,
                       double              R0_anion,
                       double              R2_anion,
                       non_local_well_type well_type)
        : m_alpha_0_cation(alpha_0_cation),
          m_beta_0_cation(beta_0_cation),
          m_A2_cation(A2_cation),
          m_alpha_0_anion(alpha_0_anion),
          m_beta_0_anion(beta_0_anion),
          m_A2_anion(A2_anion),
          m_R0_cation(R0_cation),
          m_R2_cation(R2_cation),
          m_R0_anion(R0_anion),
          m_R2_anion(R2_anion),
          m_well_type(well_type) {
        m_alpha_0_cation *= Constants::Ryd_to_eV;
        m_alpha_0_anion *= Constants::Ryd_to_eV;
        m_A2_cation *= Constants::Ryd_to_eV;
        m_A2_anion *= Constants::Ryd_to_eV;

        m_R0_anion *= Constants::angstrom_to_m;
        m_R2_anion *= Constants::angstrom_to_m;
        m_R0_cation *= Constants::angstrom_to_m;
        m_R2_cation *= Constants::angstrom_to_m;
    }

    void populate_non_local_parameters(const YAML::Node& node) {
        m_alpha_0_cation                        = Constants::Ryd_to_eV * node["alpha_0_cation"].as<double>();
        m_beta_0_cation                         = node["beta_0_cation"].as<double>();
        m_A2_cation                             = Constants::Ryd_to_eV * node["A2_cation"].as<double>();
        m_alpha_0_anion                         = Constants::Ryd_to_eV * node["alpha_0_anion"].as<double>();
        m_beta_0_anion                          = node["beta_0_anion"].as<double>();
        m_A2_anion                              = Constants::Ryd_to_eV * node["A2_anion"].as<double>();
        m_R0_cation                             = node["R0_cation"].as<double>();
        m_R2_cation                             = node["R2_cation"].as<double>();
        m_R0_anion                              = node["R0_anion"].as<double>();
        m_R2_anion                              = node["R2_anion"].as<double>();
        const std::string   well_type           = node["well_type"].as<std::string>();
        non_local_well_type non_local_well_type = non_local_well_type::unknown;
        if (well_type == "square") {
            non_local_well_type = non_local_well_type::square;
        } else if (well_type == "gaussian") {
            non_local_well_type = non_local_well_type::gaussian;
        } else {
            throw std::runtime_error("Unknown non-local well type");
        }
        m_well_type = non_local_well_type;
        m_R0_anion *= Constants::angstrom_to_m;
        m_R2_anion *= Constants::angstrom_to_m;
        m_R0_cation *= Constants::angstrom_to_m;
        m_R2_cation *= Constants::angstrom_to_m;
    }

    void print_non_local_parameters() const {
        std::cout << "alpha_0_cation: " << m_alpha_0_cation << std::endl;
        std::cout << "beta_0_cation: " << m_beta_0_cation << std::endl;
        std::cout << "A2_cation: " << m_A2_cation << std::endl;
        std::cout << "alpha_0_anion: " << m_alpha_0_anion << std::endl;
        std::cout << "beta_0_anion: " << m_beta_0_anion << std::endl;
        std::cout << "A2_anion: " << m_A2_anion << std::endl;
        std::cout << "R0_cation: " << m_R0_cation << std::endl;
        std::cout << "R2_cation: " << m_R2_cation << std::endl;
        std::cout << "R0_anion: " << m_R0_anion << std::endl;
        std::cout << "R2_anion: " << m_R2_anion << std::endl;
        std::cout << "well_type: " << static_cast<int>(m_well_type) << std::endl;
    }
};
}  // namespace uepm::pseudopotential