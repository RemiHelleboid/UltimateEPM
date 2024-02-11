/**
 * @file elelectron_phonon.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-02-09
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "electron_phonon.hpp"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <iostream>
#include <tuple>

#include "Constants.hpp"
#include "Vector3D.h"
#include "bz_states.hpp"

namespace bz_mesh {

double ElectronPhonon::bose_einstein_distribution(double energy, double temperature) {
    double f = 1.0 / (std::expm1(energy / (EmpiricalPseudopotential::Constants::k_b_eV * temperature)));
    return f;
}

double ElectronPhonon::electron_overlap_integral(const Vector3D<double>& k1, const Vector3D<double>& k2) {
    const double R_Wigner_Seitz = 2.122e-10;
    const double qRws           = (k1 - k2).Length() * R_Wigner_Seitz;
    double       integral       = 3.0 * (std::sin(qRws) - qRws * std::cos(qRws)) / (qRws * qRws * qRws);
    return integral;
}

double ElectronPhonon::hole_overlap_integral(int n1, const Vector3D<double>& k1, int n2, const Vector3D<double>& k2) {
    const double cos_angle_k1_k2_2 = compte_cos_angle(k1, k2) * compte_cos_angle(k1, k2);
    auto         A_B_params        = m_hole_overlap_int_params.get_params(n1, n2);
    double       integral          = (1.0 / 2.0) * std::sqrt(A_B_params[0] + A_B_params[1] * cos_angle_k1_k2_2);
    return integral;
}

RateValues ElectronPhonon::compute_electron_phonon_rate(int idx_n1, std::size_t idx_k1) {
    RateValues                rates_k1_n1;
    const std::vector<Tetra>& list_tetrahedra          = m_list_tetrahedra;
    auto                      indices_conduction_bands = m_indices_conduction_bands;
    double                    energy_n1_k1             = m_list_vertices[idx_k1].get_energy_at_band(idx_n1);

    for (auto&& idx_n2 : indices_conduction_bands) {
        for (auto&& tetra : list_tetrahedra) {
            auto             k_1 = m_list_vertices[idx_k1].get_position();
            auto             k_2 = tetra.compute_barycenter();
            Vector3D<double> k1{k_1.x(), k_1.y(), k_1.z()};
            Vector3D<double> k2{k_2.x(), k_2.y(), k_2.z()};

            double           overlap_integral = electron_overlap_integral(k1, k2);
            auto             q                = k1 - k2;
            Vector3D<double> q_ph{q.X, q.X, q.X};
            // for (auto&& mode_phonon : m_phonon_dispersion) {
            //     double e_ph = mode_phonon.get_phonon_dispersion(q_ph.Length());
            //     for (auto sign_phonon : {-1.0, 1.0}) {
            //         double               bose_part    = (bose_einstein_distribution(e_ph, m_temperature) + sign_phonon);
            //         double               energy_final = energy_n1_k1 + e_ph * sign_phonon;
            //         double               dos_tetra    = tetra.compute_tetra_dos_energy_band(energy_final, idx_n2);
            //         DeformationPotential deformation_potential =
            //             (mode_phonon.m_mode == PhononMode::acoustic) ? m_acoustic_deformation_potential :
            //             m_optical_deformation_potential;
            //         double          deformation_potential_value = deformation_potential.get_deformation_potential(q_ph, energy_final);
            //         PhononMode      phonon_mode                 = mode_phonon.m_mode;
            //         PhononDirection phonon_direction            = mode_phonon.m_direction;
            //         PhononEvent     phonon_event                = (sign_phonon == 1.0) ? PhononEvent::absorption : PhononEvent::emission;
            //         double          rate_value = (EmpiricalPseudopotential::Constants::pi / (m_rho * e_ph)) * deformation_potential_value
            //         *
            //                             deformation_potential_value * overlap_integral * overlap_integral * bose_part * dos_tetra;
            //         RateValue rate(phonon_mode, phonon_direction, phonon_event, rate_value);
            //         rates_k1_n1.add_rate(rate);
            //     }
            // }
            for (auto&& mode_phonon : m_phonon_dispersion) {
                PhononModeDirection mode_direction = mode_phonon.first;
                double              e_ph           = mode_phonon.second.get_phonon_dispersion(q_ph.Length());
                for (auto sign_phonon : {-1.0, 1.0}) {
                    double               bose_part    = (bose_einstein_distribution(e_ph, m_temperature) + sign_phonon);
                    double               energy_final = energy_n1_k1 + e_ph * sign_phonon;
                    double               dos_tetra    = tetra.compute_tetra_dos_energy_band(energy_final, idx_n2);
                    DeformationPotential deformation_potential =
                        (mode_direction.first == PhononMode::acoustic) ? m_acoustic_deformation_potential : m_optical_deformation_potential;
                    double          deformation_potential_value = deformation_potential.get_deformation_potential(q_ph, energy_final);
                    PhononMode      phonon_mode                 = mode_direction.first;
                    PhononDirection phonon_direction            = mode_direction.second;
                    PhononEvent     phonon_event                = (sign_phonon == 1.0) ? PhononEvent::absorption : PhononEvent::emission;
                    double          rate_value = (EmpiricalPseudopotential::Constants::pi / (m_rho * e_ph)) * deformation_potential_value *
                                        deformation_potential_value * overlap_integral * overlap_integral * bose_part * dos_tetra;
                    RateValue rate(phonon_mode, phonon_direction, phonon_event, rate_value);
                    rates_k1_n1.add_rate(rate);
                }
            }
        }
    }
    rates_k1_n1.print_rates();
    return rates_k1_n1;
}

void ElectronPhonon::compute_electron_phonon_rates_over_mesh() {
    for (std::size_t idx_k1 = 0; idx_k1 < m_list_vertices.size(); ++idx_k1) {
        std::cout << "Computing rates for k-point " << idx_k1 << std::endl;
        for (std::size_t idx_n1 = 0; idx_n1 < m_indices_valence_bands.size(); ++idx_n1) {
            auto rate = compute_electron_phonon_rate(idx_n1, idx_k1);
            m_list_vertices[idx_k1].add_electron_phonon_rates(rate.to_array());
        }
    }

    // Add rate to the mesh
}

void ElectronPhonon::compute_plot_electron_phonon_rates_vs_energy_over_mesh(int                nb_bands,
                                                                            double             max_energy,
                                                                            double             energy_step,
                                                                            const std::string& filename) {
    std::vector<double>              energies;
    std::vector<std::vector<double>> rates;
    for (double energy = 0.0; energy < max_energy; energy += energy_step) {
        energies.push_back(energy);
        for (int idx_band = 0; idx_band < nb_bands; ++idx_band) {
            double rate = 0.0;
            for (auto&& vertex : m_list_vertices) {
                for (std::size_t idx_band = 0; idx_band < vertex.get_number_bands(); ++idx_band) {
                    rate += vertex.get_energy_at_band(idx_band);
                }
            }
            // rates[idx_band].push_back(rate);
        }
    }
}

void ElectronPhonon::load_phonon_parameters(const std::string& filename) {
    YAML::Node config = YAML::LoadFile(filename);

    // Check if the file is empty
    if (config.IsNull()) {
        throw std::runtime_error("File " + filename + " is empty");
    }
    // Print the file
    std::cout << "File " << filename << " contains: " << std::endl;
    std::cout << config << std::endl;

    auto list_materials = config["materials"];

    const std::string& my_material = m_material.get_name();

    auto same_material = [my_material](const YAML::Node& node) { return node["name"].as<std::string>() == my_material; };
    auto it_material   = std::find_if(list_materials.begin(), list_materials.end(), same_material);
    


    if (it_material == list_materials.end()) {
        throw std::runtime_error("Material " + my_material + " not found in file " + filename);
    }
    auto material = *it_material;

    double radiusWS      = material["Radius-WS"].as<double>();
    m_radii_wigner_seitz = radiusWS;

    // Example for reading nested dispersion properties
    auto dispersion = material["dispersion"];
    for (const auto& type : {"longitudinal", "transverse"}) {
        auto dispType = dispersion[type];
        for (const auto& wave : {"acoustic", "optic"}) {
            auto   waveType = dispType[wave];
            double w0       = waveType["w0"].as<double>();
            double vs       = waveType["vs"].as<double>();
            double c        = waveType["c"].as<double>();
            std::cout << "w0: " << w0 << " vs: " << vs << " c: " << c << std::endl;

            PhononDirection  direction = (type == "longitudinal") ? PhononDirection::longitudinal : PhononDirection::transverse;
            PhononMode       mode      = (wave == "acoustic") ? PhononMode::acoustic : PhononMode::optical;
            PhononDispersion phononDispersion(mode, direction, w0, vs, c);
            m_phonon_dispersion[std::make_pair(mode, direction)] = phononDispersion;
        }
    }

    // Example for reading deformation potential
    auto node_deformationPotential = material["deformation-potential"];
    for (const auto& carrierType : {"electron", "hole"}) {
        auto carrier = node_deformationPotential[carrierType];
        for (const auto& wave : {"acoustic", "optic"}) {
            auto   waveType = carrier[wave];
            double A        = waveType["A"].as<double>();
            double B        = waveType["B"].as<double>();
            std::cout << "A: " << A << " B: " << B << std::endl;

            PhononMode           mode = (wave == "acoustic") ? PhononMode::acoustic : PhononMode::optical;
            DeformationPotential deformationPotential(mode, A, B);
            if (mode == PhononMode::acoustic) {
                m_acoustic_deformation_potential = deformationPotential;
            } else {
                m_optical_deformation_potential = deformationPotential;
            }
        }
    }
}

}  // namespace bz_mesh