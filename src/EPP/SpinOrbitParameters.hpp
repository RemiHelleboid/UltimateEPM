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


struct SpinOrbitParameters {
    /**
     * @brief Length scale of the radial wave function for the cation.
     * 
     */
    double m_zeta_cation;

    /**
     * @brief Length scale of the radial wave function for the anion.
     * 
     */
    double m_zeta_anion;

    /**
     * @brief Fitting parameter alpha.
     * 
     */
    double m_alpha;

    /**
     * @brief Fitting parameter beta.
     * 
     */
    double m_mu;
    
    /**
     * @brief Default constructor.
     * 
     */
    SpinOrbitParameters() = default;

    /**
     * @brief Populate the spin-orbit parameters from a YAML node.
     * 
     * @param node 
     */
    void populate_from_yaml(const YAML::Node& node) {
        m_zeta_cation = node["zeta_cation"].as<double>();
        m_zeta_anion  = node["zeta_anion"].as<double>();
        m_alpha       = node["alpha"].as<double>();
        m_mu          = node["mu"].as<double>();
    }

    /**
     * @brief Display the spin-orbit parameters.
     * 
     */
    void print() const {
        std::cout << "zeta_cation = " << m_zeta_cation << std::endl;
        std::cout << "zeta_anion = " << m_zeta_anion << std::endl;
        std::cout << "alpha = " << m_alpha << std::endl;
        std::cout << "mu = " << m_mu << std::endl;
    }
};