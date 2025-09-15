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

// Order: ac/opt × L/T × ab/em  → indices 0..7
constexpr int rate_index(PhononMode m, PhononDirection d, PhononEvent e) {
    int M = (m == PhononMode::acoustic ? 0 : m == PhononMode::optical ? 1 : -1);
    int D = (d == PhononDirection::longitudinal ? 0 : d == PhononDirection::transverse ? 1 : -1);
    int E = (e == PhononEvent::absorption ? 0 : e == PhononEvent::emission ? 1 : -1);
    return (M < 0 || D < 0 || E < 0) ? -1 : ((M * 2 + D) * 2 + E);
}

// Named indices (handy for headers/export)
constexpr int IDX_AC_L_AB = rate_index(PhononMode::acoustic, PhononDirection::longitudinal, PhononEvent::absorption);  // 0
constexpr int IDX_AC_T_AB = rate_index(PhononMode::acoustic, PhononDirection::transverse, PhononEvent::absorption);    // 1
constexpr int IDX_OP_L_AB = rate_index(PhononMode::optical, PhononDirection::longitudinal, PhononEvent::absorption);   // 2
constexpr int IDX_OP_T_AB = rate_index(PhononMode::optical, PhononDirection::transverse, PhononEvent::absorption);     // 3
constexpr int IDX_AC_L_EM = rate_index(PhononMode::acoustic, PhononDirection::longitudinal, PhononEvent::emission);    // 4
constexpr int IDX_AC_T_EM = rate_index(PhononMode::acoustic, PhononDirection::transverse, PhononEvent::emission);      // 5
constexpr int IDX_OP_L_EM = rate_index(PhononMode::optical, PhononDirection::longitudinal, PhononEvent::emission);     // 6
constexpr int IDX_OP_T_EM = rate_index(PhononMode::optical, PhononDirection::transverse, PhononEvent::emission);       // 7

static_assert(IDX_AC_L_AB == 0 && IDX_OP_T_EM == 7, "Rate index mapping changed; update I/O accordingly.");

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
            std::cout << "vs*q: " << m_vs * q << std::endl;
            std::cout << "c*q^2: " << m_c * q * q << std::endl;
            throw std::runtime_error("Negative phonon energy");
        }
        return e_ph;
    }
};

struct DeformationPotential {
    PhononMode m_mode             = PhononMode::none;
    double     m_A                = 0.0;
    double     m_B                = 0.0;
    double     m_energy_threshold = 1e6;  // eV, default no threshold

    DeformationPotential() = default;
    DeformationPotential(PhononMode type, double A, double B) : m_mode(type), m_A(A), m_B(B), m_energy_threshold(1e6) {}

    double get_deformation_potential(const vector3& q, double energy) const {
        double clamp_energy = std::min(energy, m_energy_threshold);
        if (m_mode == PhononMode::acoustic) {
            return std::sqrt(m_A + clamp_energy * m_B) * q.norm();
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
    std::array<double, 8> values{};  // zero-initialized

    void print_rates() const {
        std::cout << "Acoustic longitudinal absorption : " << values[IDX_AC_L_AB] << "\n"
                  << "Acoustic transverse absorption   : " << values[IDX_AC_T_AB] << "\n"
                  << "Optical longitudinal absorption  : " << values[IDX_OP_L_AB] << "\n"
                  << "Optical transverse absorption    : " << values[IDX_OP_T_AB] << "\n"
                  << "Acoustic longitudinal emission   : " << values[IDX_AC_L_EM] << "\n"
                  << "Acoustic transverse emission     : " << values[IDX_AC_T_EM] << "\n"
                  << "Optical longitudinal emission    : " << values[IDX_OP_L_EM] << "\n"
                  << "Optical transverse emission      : " << values[IDX_OP_T_EM] << std::endl;
    }

    [[nodiscard]] std::array<double, 8> to_array() const { return values; }

    void add_rate(const RateValue& rate) {
        if (int idx = rate_index(rate.m_mode, rate.m_direction, rate.m_type); idx >= 0) values[idx] += rate.m_value;
    }

    // Optional: direct accessors if you like named use sites
    double&       at(PhononMode m, PhononDirection d, PhononEvent e) { return values[rate_index(m, d, e)]; }
    const double& at(PhononMode m, PhononDirection d, PhononEvent e) const { return values[rate_index(m, d, e)]; }

    inline static const char* rate_label(int i) {
        static const char* L[8] = {"ac_L_ab", "ac_T_ab", "op_L_ab", "op_T_ab", "ac_L_em", "ac_T_em", "op_L_em", "op_T_em"};
        return (i >= 0 && i < 8) ? L[i] : "invalid";
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

    inline double bose_einstein_distribution(double energy, double temperature);
    double        electron_overlap_integral(const vector3& k1, const vector3& k2);
    double        hole_overlap_integral(int n1, const vector3& k1, int n2, const vector3& k2);

    RateValues compute_electron_phonon_rate(int idx_n1, std::size_t idx_k1);
    RateValues compute_hole_phonon_rate(int idx_n1, std::size_t idx_k1);

    void set_temperature(double temperature) { m_temperature = temperature; }
    void set_density(double rho) { m_rho = rho; }

    void compute_electron_phonon_rates_over_mesh();
    void add_electron_phonon_rates_to_mesh(const std::string& initial_filename, const std::string& final_filename);

    void export_rate_values(const std::string& filename) const;

    void compute_plot_electron_phonon_rates_vs_energy_over_mesh(int                nb_bands,
                                                                double             max_energy,
                                                                double             energy_step,
                                                                const std::string& filename);
};

}  // namespace bz_mesh