/**
 * @file SpinOrbitParameters.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-09-30
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

namespace EmpiricalPseudopotential {

/**
 * @brief Struct to store the Spin-Orbit Coupling (SOC) parameters.
 *
 */
struct SpinOrbitParameters {
    /**
     * @brief Period of the anion (line number in the periodic table).
     *
     */
    unsigned int m_period_anion = 0;

    /**
     * @brief Period of the cation (line number in the periodic table).
     *
     */
    unsigned int m_period_cation = 0;

    /**
     * @brief Length scale of the radial wave function for the cation.
     *
     */
    double m_radial_extent_cation = 0;

    /**
     * @brief Length scale of the radial wave function for the anion.
     *
     */
    double m_radial_extent_anion = 0;

    /**
     * @brief Fitting parameter alpha.
     *
     */
    double m_alpha = 0;

    /**
     * @brief Fitting parameter beta.
     *
     */
    double m_mu = 0;

    /**
     * @brief Default constructor.
     *
     */
    SpinOrbitParameters() = default;

    // Copy constructor
    SpinOrbitParameters(const SpinOrbitParameters& other) = default;
    SpinOrbitParameters& operator=(const SpinOrbitParameters& other) = default;
    SpinOrbitParameters(SpinOrbitParameters&& other)                 = default;

    /**
     * @brief Populate the spin-orbit parameters from a YAML node.
     *
     * @param node
     */
    void populate_from_yaml(const YAML::Node& node) {
        m_period_anion         = node["period_anion"].as<unsigned int>();
        m_period_cation        = node["period_cation"].as<unsigned int>();
        m_radial_extent_anion  = node["radial_extent_anion"].as<double>();
        m_radial_extent_cation = node["radial_extent_cation"].as<double>();
        m_alpha                = node["alpha_soc"].as<double>();
        m_mu                   = 1.0 * Constants::Ryd_to_eV * node["mu_soc"].as<double>();
        // print_parameters();
    }

    /**
     * @brief Display the spin-orbit parameters.
     *
     */
    void print_parameters() const {
        std::cout << "Spin-Orbit Coupling parameters:" << std::endl;
        std::cout << "Period anion: " << m_period_anion << std::endl;
        std::cout << "Period cation: " << m_period_cation << std::endl;
        std::cout << "Radial extent anion: " << m_radial_extent_anion << std::endl;
        std::cout << "Radial extent cation: " << m_radial_extent_cation << std::endl;
        std::cout << "Alpha: " << m_alpha << std::endl;
        std::cout << "Mu: " << m_mu << std::endl;
    }
};

}  // namespace EmpiricalPseudopotential