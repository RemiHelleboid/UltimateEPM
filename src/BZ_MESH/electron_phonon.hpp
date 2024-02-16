/**
 * @file elelectron_phonon.hpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-02-09
 *
 * @copyright Copyright (c) 2024
 *
 * Phonon are always stored in the same order:
 * 0: acoustic longitudinal absorption
 * 1: acoustic transverse absorption
 * 2: optical longitudinal absorption
 * 3: optical transverse absorption
 *
 */

#pragma once

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

#include "bz_states.hpp"
#include "yaml-cpp/yaml.h"

namespace bz_mesh {

enum class PhononMode { acoustic, optical, none };
enum class PhononDirection { longitudinal, transverse, none };
enum class PhononEvent { absorption, emission, none };

struct PhononDispersion {
    PhononMode      m_mode      = PhononMode::none;
    PhononDirection m_direction = PhononDirection::none;
    double          m_omega     = 0.0;
    double          m_vs        = 0.0;
    double          m_c         = 0.0;

    PhononDispersion() = default;

    PhononDispersion(PhononMode mode, PhononDirection direction, double omega, double vs, double c)
        : m_mode(mode),
          m_direction(direction),
          m_omega(omega),
          m_vs(vs),
          m_c(c) {}
    PhononDispersion(PhononMode mode, PhononDirection direction)
        : m_mode(mode),
          m_direction(direction),
          m_omega(0.0),
          m_vs(0.0),
          m_c(0.0) {}

    double get_phonon_dispersion(double q) const {
        double e_ph = m_omega + m_vs * q + m_c * q * q;
        if (e_ph < 0.0) {
            std::cout << "ERROR: negative phonon energy: " << e_ph << std::endl;
            std::cout << "Q: " << q << std::endl;
            std::cout << "OMEGA: " << m_omega << std::endl;
            std::cout << "VS: " << m_vs << std::endl;
            std::cout << "C: " << m_c << std::endl;
        }
        return e_ph;
    }
};

struct DeformationPotential {
    PhononMode m_mode             = PhononMode::none;
    double     m_A                = 0.0;
    double     m_B                = 0.0;
    double     m_energy_threshold = 0.0;

    DeformationPotential() = default;
    DeformationPotential(PhononMode type, double A, double B) : m_mode(type), m_A(A), m_B(B) {}

    double get_deformation_potential(Vector3D<double> q, double energy) const {
        double clamp_energy = std::min(energy, m_energy_threshold);
        if (m_mode == PhononMode::acoustic) {
            return std::sqrt(m_A + clamp_energy * m_B) * q.Length();
        } else {
            return std::sqrt(m_A + clamp_energy * m_B);
        }
    }
};

struct HoleOverlapIntParams {
    std::vector<std::pair<std::pair<int, int>, double>> A_params =
        {{{1, 1}, 1.0}, {{2, 2}, 1.0}, {{3, 3}, 5.0 / 8.0}, {{1, 2}, 3.0}, {{1, 3}, 3.0 / 8.0}, {{2, 3}, 3.0 / 8.0}};

    std::vector<std::pair<std::pair<int, int>, double>> B_params =
        {{{1, 1}, 3.0}, {{2, 2}, 3.0}, {{3, 3}, 0.0}, {{1, 2}, -3.0}, {{1, 3}, 0.0}, {{2, 3}, 0.0}};

    /**
     * @brief Return the parameters for the overlap integral, given the indices of the bands.
     *
     * @param n1
     * @param n2
     * @return std::array<double, 2>
     */
    std::array<double, 2> get_params(int n1, int n2) const {
        for (const auto& p : A_params) {
            if (p.first.first == n1 && p.first.second == n2) {
                return {p.second, 0.0};
            }
        }
        for (const auto& p : B_params) {
            if (p.first.first == n1 && p.first.second == n2) {
                return {0.0, p.second};
            }
        }
        return {0.0, 0.0};
    }

    HoleOverlapIntParams() = default;
};

struct RateValue {
    PhononMode      m_mode      = PhononMode::none;
    PhononDirection m_direction = PhononDirection::none;
    PhononEvent     m_type      = PhononEvent::none;
    double          m_value     = 0.0;

    RateValue(PhononMode mode, PhononDirection direction, PhononEvent type, double value)
        : m_mode(mode),
          m_direction(direction),
          m_type(type),
          m_value(value) {}
};

struct RateValues {
    RateValue m_ac_long_absorption   = RateValue(PhononMode::acoustic, PhononDirection::longitudinal, PhononEvent::absorption, 0.0);
    RateValue m_ac_trans_absorption  = RateValue(PhononMode::acoustic, PhononDirection::transverse, PhononEvent::absorption, 0.0);
    RateValue m_opt_long_absorption  = RateValue(PhononMode::optical, PhononDirection::longitudinal, PhononEvent::absorption, 0.0);
    RateValue m_opt_trans_absorption = RateValue(PhononMode::optical, PhononDirection::transverse, PhononEvent::absorption, 0.0);

    RateValue m_ac_long_emission   = RateValue(PhononMode::acoustic, PhononDirection::longitudinal, PhononEvent::emission, 0.0);
    RateValue m_ac_trans_emission  = RateValue(PhononMode::acoustic, PhononDirection::transverse, PhononEvent::emission, 0.0);
    RateValue m_opt_long_emission  = RateValue(PhononMode::optical, PhononDirection::longitudinal, PhononEvent::emission, 0.0);
    RateValue m_opt_trans_emission = RateValue(PhononMode::optical, PhononDirection::transverse, PhononEvent::emission, 0.0);

    /**
     * @brief Print the rates to the standard output.
     *
     */
    void print_rates() const {
        std::cout << "Acoustic longitudinal absorption: " << m_ac_long_absorption.m_value << std::endl;
        std::cout << "Acoustic transverse absorption: " << m_ac_trans_absorption.m_value << std::endl;
        std::cout << "Optical longitudinal absorption: " << m_opt_long_absorption.m_value << std::endl;
        std::cout << "Optical transverse absorption: " << m_opt_trans_absorption.m_value << std::endl;
        std::cout << "Acoustic longitudinal emission: " << m_ac_long_emission.m_value << std::endl;
        std::cout << "Acoustic transverse emission: " << m_ac_trans_emission.m_value << std::endl;
        std::cout << "optical longitudinal emission: " << m_opt_long_emission.m_value << std::endl;
        std::cout << "Optical transverse emission: " << m_opt_trans_emission.m_value << std::endl;
    }

    std::array<double, 8> to_array() const {
        return {m_ac_long_absorption.m_value,
                m_ac_trans_absorption.m_value,
                m_opt_long_absorption.m_value,
                m_opt_trans_absorption.m_value,
                m_ac_long_emission.m_value,
                m_ac_trans_emission.m_value,
                m_opt_long_emission.m_value,
                m_opt_trans_emission.m_value};
    }

    void add_rate(RateValue rate) {
        if (rate.m_mode == PhononMode::acoustic) {
            if (rate.m_direction == PhononDirection::longitudinal) {
                if (rate.m_type == PhononEvent::absorption) {
                    m_ac_long_absorption.m_value += rate.m_value;
                } else {
                    m_ac_long_emission.m_value += rate.m_value;
                }
            } else {
                if (rate.m_type == PhononEvent::absorption) {
                    m_ac_trans_absorption.m_value += rate.m_value;
                } else {
                    m_ac_trans_emission.m_value += rate.m_value;
                }
            }
        } else {
            if (rate.m_direction == PhononDirection::longitudinal) {
                if (rate.m_type == PhononEvent::absorption) {
                    m_opt_long_absorption.m_value += rate.m_value;
                } else {
                    m_opt_long_emission.m_value += rate.m_value;
                }
            } else {
                if (rate.m_type == PhononEvent::absorption) {
                    m_opt_trans_absorption.m_value += rate.m_value;
                } else {
                    m_opt_trans_emission.m_value += rate.m_value;
                }
            }
        }
    }
};

typedef std::pair<PhononMode, PhononDirection> PhononModeDirection;

class ElectronPhonon : public BZ_States {
 private:
    double               m_temperature        = 300.0;
    double               m_rho                = 2.329e3;
    double               m_radii_wigner_seitz = 0.0;
    HoleOverlapIntParams m_hole_overlap_int_params;
    DeformationPotential m_acoustic_deformation_potential_e;
    DeformationPotential m_optical_deformation_potential_e;
    DeformationPotential m_acoustic_deformation_potential_h;
    DeformationPotential m_optical_deformation_potential_h;

    std::map<PhononModeDirection, PhononDispersion> m_phonon_dispersion;

 public:
    explicit ElectronPhonon(const EmpiricalPseudopotential::Material& material) : BZ_States(material) {}

    void load_phonon_parameters(const std::string& filename);
    void plot_phonon_dispersion(const std::string& filename) const;

    double bose_einstein_distribution(double energy, double temperature);
    double electron_overlap_integral(const Vector3D<double>& k1, const Vector3D<double>& k2);
    double hole_overlap_integral(int n1, const Vector3D<double>& k1, int n2, const Vector3D<double>& k2);

    RateValues compute_electron_phonon_rate(int idx_n1, std::size_t idx_k1);
    RateValues compute_hole_phonon_rate(int idx_n1, std::size_t idx_k1);

    void set_temperature(double temperature) { m_temperature = temperature; }
    void set_density(double rho) { m_rho = rho; }

    void compute_electron_phonon_rates_over_mesh();

    void export_rate_values(const std::string& filename) const;

    void compute_plot_electron_phonon_rates_vs_energy_over_mesh(int                nb_bands,
                                                                double             max_energy,
                                                                double             energy_step,
                                                                const std::string& filename);
};

}  // namespace bz_mesh