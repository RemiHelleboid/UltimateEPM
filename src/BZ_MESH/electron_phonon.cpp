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
#include <fstream>
#include <iostream>
#include <random>
#include <tuple>

#include "Constants.hpp"
#include "Vector3D.h"
#include "bz_states.hpp"

#include "omp.h"

namespace bz_mesh {

double ElectronPhonon::bose_einstein_distribution(double energy, double temperature) {
    if (energy < 1e-9) {
        return 0.0;
    }
    double f = 1.0 / (std::expm1(energy / (EmpiricalPseudopotential::Constants::k_b_eV * temperature)));
    if (f < 0 || std::isnan(f) || std::isinf(f)) {
        std::cout << "Energy: " << energy << " Temperature: " << temperature << " f: " << f << std::endl;
        throw std::runtime_error("Bose-Einstein distribution is negative or NaN");
    }
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
            auto             q                = k2 - k1;
            auto             initial_q        = q;
            Vector3D<double> q_ph{q.X, q.X, q.X};
            bool             is_in_bz = is_inside_mesh_geometry(q_ph / m_material.get_fourier_factor());

            // No Umklapp process for now.
            if (!is_in_bz) {
                auto new_q = retrieve_k_inside_mesh_geometry(q_ph / m_material.get_fourier_factor()) * m_material.get_fourier_factor();
                q_ph       = Vector3D<double>{new_q.x(), new_q.y(), new_q.z()};
                // std::cout << "q: " << q_ph / m_material.get_fourier_factor() << " new_q: " << new_q / m_material.get_fourier_factor() <<
                // std::endl;
            }
            is_in_bz = is_inside_mesh_geometry(q_ph / m_material.get_fourier_factor());
            if (!is_in_bz) {
                throw std::runtime_error("Q is not inside the BZ");
            }

            for (auto&& mode_phonon : m_phonon_dispersion) {
                PhononModeDirection mode_direction = mode_phonon.first;
                double e_ph = mode_phonon.second.get_phonon_dispersion(q_ph.Length()) * EmpiricalPseudopotential::Constants::h_bar_eV;
                if (e_ph > 100e-3) {
                    std::cout << "Energy phonon: " << e_ph << std::endl << std::endl;
                    throw std::runtime_error("Energy phonon too high");
                }
                if (e_ph <= 0.0) {
                    // std::cout << "Energy phonon: " << e_ph << std::endl << std::endl;
                    // std::cout << "Q: " << q_ph << std::endl << std::endl;
                    // std::cout << "Q: " << q_ph.Length() / m_material.get_fourier_factor() << std::endl << std::endl;
                    // std::cout << "Energy (n1, k1): " << energy_n1_k1 << std::endl << std::endl;
                    bool is_in_bz = is_inside_mesh_geometry(q_ph / m_material.get_fourier_factor());
                    // std::cout << "Is in BZ: " << is_in_bz << std::endl << std::endl;
                    // throw std::runtime_error("Energy phonon negative");
                    continue;
                }
                // std::cout << "Energy phonon: " << e_ph << std::endl;
                for (auto sign_phonon : {-1.0, 1.0}) {
                    double bose_part    = (bose_einstein_distribution(e_ph, m_temperature) + 0.5 + 0.5 * std::pow(-1.0, sign_phonon));
                    double energy_final = energy_n1_k1 + e_ph * sign_phonon;
                    // std::cout << "Energy (n1, k1): " << energy_n1_k1 << " Energy phonon: " << e_ph << " Energy final: " << energy_final
                    //           << std::endl;
                    if (!tetra.is_energy_inside_band(energy_final, idx_n2)) {
                        continue;
                    }
                    // std::cout << "Volume: " << this->get_volume() << std::endl;
                    double dos_tetra = tetra.compute_tetra_dos_energy_band(energy_final, idx_n2);
                    // std::cout << "Energy final: " << energy_final << " dos tetra: " << dos_tetra
                    //           << std::endl;
                    DeformationPotential deformation_potential       = (mode_direction.first == PhononMode::acoustic)
                                                                           ? m_acoustic_deformation_potential_e
                                                                           : m_optical_deformation_potential_e;
                    double               deformation_potential_value = deformation_potential.get_deformation_potential(q_ph, energy_final);
                    PhononMode           phonon_mode                 = mode_direction.first;
                    PhononDirection      phonon_direction            = mode_direction.second;
                    PhononEvent          phonon_event = (sign_phonon == 1.0) ? PhononEvent::absorption : PhononEvent::emission;
                    double rate_value = (EmpiricalPseudopotential::Constants::pi / (m_rho * e_ph)) * deformation_potential_value *
                                        deformation_potential_value * overlap_integral * overlap_integral * bose_part * dos_tetra;
                    // std::cout << this->get_volume() << std::endl;
                    rate_value /= this->get_volume();
                    rate_value *= EmpiricalPseudopotential::Constants::q;

                    if (rate_value < 0.0 || rate_value > 1e50 || std::isnan(rate_value) || std::isinf(rate_value)) {
                        std::cout << "Rate value: " << rate_value << std::endl;
                        std::cout << "Overlap integral: " << overlap_integral << std::endl;
                        std::cout << "Deformation potential: " << deformation_potential_value << std::endl;
                        std::cout << "Bose part: " << bose_part << std::endl;
                        std::cout << "DOS tetra: " << dos_tetra << std::endl;
                        std::cout << "Energy final: " << energy_final << std::endl;
                        std::cout << "q: " << q_ph << std::endl;
                        std::cout << "initial q: " << initial_q / m_material.get_fourier_factor() << std::endl;
                        std::cout << "Energy phonon: " << e_ph << std::endl;
                        std::cout << "Energy (n1, k1): " << energy_n1_k1 << std::endl;
                        std::cout << "--------------------------------------------------------------------------------" << std::endl;
                    }
                    // std::cout << "Rate value: " << (deformation_potential_value) << std::endl;
                    RateValue rate(phonon_mode, phonon_direction, phonon_event, rate_value);
                    rates_k1_n1.add_rate(rate);
                }
            }
        }
    }
    // rates_k1_n1.print_rates();
    // std::cout << "*************************************************************************************************************" <<
    // std::endl;
    return rates_k1_n1;
}

RateValues ElectronPhonon::compute_hole_phonon_rate(int idx_n1, std::size_t idx_k1) {
    RateValues                rates_k1_n1;
    const std::vector<Tetra>& list_tetrahedra       = m_list_tetrahedra;
    auto                      indices_valence_bands = m_indices_valence_bands;
    double                    energy_n1_k1          = m_list_vertices[idx_k1].get_energy_at_band(idx_n1);

    for (auto&& idx_n2 : indices_valence_bands) {
        for (auto&& tetra : list_tetrahedra) {
            auto             k_1 = m_list_vertices[idx_k1].get_position();
            auto             k_2 = tetra.compute_barycenter();
            Vector3D<double> k1{k_1.x(), k_1.y(), k_1.z()};
            Vector3D<double> k2{k_2.x(), k_2.y(), k_2.z()};

            double           overlap_integral = hole_overlap_integral(idx_n1, k1, idx_n2, k2);
            auto             q                = k2 - k1;
            auto             initial_q        = q;
            Vector3D<double> q_ph{q.X, q.X, q.X};
            bool             is_in_bz = is_inside_mesh_geometry(q_ph / m_material.get_fourier_factor());

            // No Umklapp process for now.
            if (!is_in_bz) {
                auto new_q = retrieve_k_inside_mesh_geometry(q_ph / m_material.get_fourier_factor()) * m_material.get_fourier_factor();
                q_ph       = Vector3D<double>{new_q.x(), new_q.y(), new_q.z()};
                // std::cout << "q: " << q_ph / m_material.get_fourier_factor() << " new_q: " << new_q / m_material.get_fourier_factor() <<
                // std::endl;
            }
            is_in_bz = is_inside_mesh_geometry(q_ph / m_material.get_fourier_factor());
            if (!is_in_bz) {
                throw std::runtime_error("Q is not inside the BZ");
            }

            for (auto&& mode_phonon : m_phonon_dispersion) {
                PhononModeDirection mode_direction = mode_phonon.first;
                double e_ph = mode_phonon.second.get_phonon_dispersion(q_ph.Length()) * EmpiricalPseudopotential::Constants::h_bar_eV;
                if (e_ph > 100e-3) {
                    std::cout << "Energy phonon: " << e_ph << std::endl << std::endl;
                    throw std::runtime_error("Energy phonon too high");
                }
                if (e_ph <= 0.0) {
                    // std::cout << "Energy phonon: " << e_ph << std::endl << std::endl;
                    // std::cout << "Q: " << q_ph << std::endl << std::endl;
                    // std::cout << "Q: " << q_ph.Length() / m_material.get_fourier_factor() << std::endl << std::endl;
                    // std::cout << "Energy (n1, k1): " << energy_n1_k1 << std::endl << std::endl;
                    // bool is_in_bz = is_inside_mesh_geometry(q_ph / m_material.get_fourier_factor());
                    // std::cout << "Is in BZ: " << is_in_bz << std::endl << std::endl;
                    // throw std::runtime_error("Energy phonon negative");
                    continue;
                }
                // std::cout << "Energy phonon: " << e_ph << std::endl;
                for (auto sign_phonon : {-1.0, 1.0}) {
                    double bose_part    = (bose_einstein_distribution(e_ph, m_temperature) + 0.5 + 0.5 * std::pow(-1.0, sign_phonon));
                    double energy_final = energy_n1_k1 + e_ph * sign_phonon;
                    // std::cout << "Energy (n1, k1): " << energy_n1_k1 << " Energy phonon: " << e_ph << " Energy final: " << energy_final
                    //           << std::endl;
                    if (!tetra.is_energy_inside_band(energy_final, idx_n2)) {
                        continue;
                    }
                    // std::cout << "Volume: " << this->get_volume() << std::endl;
                    double dos_tetra = tetra.compute_tetra_dos_energy_band(energy_final, idx_n2);
                    // std::cout << "Energy final: " << energy_final << " dos tetra: " << dos_tetra
                    //           << std::endl;
                    DeformationPotential deformation_potential       = (mode_direction.first == PhononMode::acoustic)
                                                                           ? m_acoustic_deformation_potential_h
                                                                           : m_optical_deformation_potential_h;
                    double               deformation_potential_value = deformation_potential.get_deformation_potential(q_ph, energy_final);
                    PhononMode           phonon_mode                 = mode_direction.first;
                    PhononDirection      phonon_direction            = mode_direction.second;
                    PhononEvent          phonon_event = (sign_phonon == 1.0) ? PhononEvent::absorption : PhononEvent::emission;
                    double rate_value = (EmpiricalPseudopotential::Constants::pi / (m_rho * e_ph)) * deformation_potential_value *
                                        deformation_potential_value * overlap_integral * overlap_integral * bose_part * dos_tetra;
                    // std::cout << this->get_volume() << std::endl;
                    rate_value /= this->get_volume();
                    rate_value *= EmpiricalPseudopotential::Constants::q;

                    if (rate_value < 0.0 || rate_value > 1e50 || std::isnan(rate_value) || std::isinf(rate_value)) {
                        std::cout << "Rate value: " << rate_value << std::endl;
                        std::cout << "Overlap integral: " << overlap_integral << std::endl;
                        std::cout << "Deformation potential: " << deformation_potential_value << std::endl;
                        std::cout << "Bose part: " << bose_part << std::endl;
                        std::cout << "DOS tetra: " << dos_tetra << std::endl;
                        std::cout << "Energy final: " << energy_final << std::endl;
                        std::cout << "q: " << q_ph << std::endl;
                        std::cout << "initial q: " << initial_q / m_material.get_fourier_factor() << std::endl;
                        std::cout << "Energy phonon: " << e_ph << std::endl;
                        std::cout << "Energy (n1, k1): " << energy_n1_k1 << std::endl;
                        std::cout << "--------------------------------------------------------------------------------" << std::endl;
                    }
                    // std::cout << "Rate value: " << (deformation_potential_value) << std::endl;
                    RateValue rate(phonon_mode, phonon_direction, phonon_event, rate_value);
                    rates_k1_n1.add_rate(rate);
                }
            }
        }
    }
    // rates_k1_n1.print_rates();
    // std::cout << "*************************************************************************************************************" <<
    // std::endl;
    return rates_k1_n1;
}

void ElectronPhonon::compute_electron_phonon_rates_over_mesh() {
    auto indices_conduction_bands = m_indices_conduction_bands;
    auto min_idx_conduction_band  = *std::min_element(indices_conduction_bands.begin(), indices_conduction_bands.end());
    auto max_idx_conduction_band  = *std::max_element(indices_conduction_bands.begin(), indices_conduction_bands.end());
    std::cout << "Min index conduction band: " << min_idx_conduction_band << std::endl;
    std::ofstream file("rates.txt");

    std::random_device                     rd;
    std::mt19937                           gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    double                                 p_compute_rate = 1.0;

#pragma omp parallel for schedule(dynamic)
    for (std::size_t idx_k1 = 0; idx_k1 < m_list_vertices.size(); ++idx_k1) {
        double r = dis(gen);
        for (std::size_t idx_n1 = 0; idx_n1 < min_idx_conduction_band; ++idx_n1) {
            if (r > p_compute_rate) {
                std::array<double, 8> array = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                m_list_vertices[idx_k1].add_electron_phonon_rates(array);
                continue;
            }
            auto hole_rate = compute_hole_phonon_rate(idx_n1, idx_k1);
            auto array     = hole_rate.to_array();
            m_list_vertices[idx_k1].add_electron_phonon_rates(array);
        }

        for (std::size_t idx_n1 = min_idx_conduction_band; idx_n1 < max_idx_conduction_band; ++idx_n1) {
            if (r > p_compute_rate) {
                std::array<double, 8> array = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                m_list_vertices[idx_k1].add_electron_phonon_rates(array);
                continue;
            } else {
                auto rate  = compute_electron_phonon_rate(idx_n1, idx_k1);
                auto array = rate.to_array();
                m_list_vertices[idx_k1].add_electron_phonon_rates(array);
            }
        }
        if (omp_get_thread_num() == 0){
            std::cout << "\rComputing rates for k-point " << idx_k1 << std::flush;}
             
    }
    std::cout << std::endl;
    file.close();
    // Add rate to the mesh
}

void ElectronPhonon::compute_plot_electron_phonon_rates_vs_energy_over_mesh(int                nb_bands,
                                                                            double             max_energy,
                                                                            double             energy_step,
                                                                            const std::string& filename) {
    std::ofstream                      file(filename);
    std::vector<double>                energies;
    std::vector<double>                list_dos;
    std::vector<std::array<double, 8>> list_mean_rates;
    double                             total_dos = 0.0;
    for (double energy = 1.0; energy < max_energy; energy += energy_step) {
        std::cout << "\rEnergy: " << energy << " Max energy: " << max_energy << std::flush;
        double                dos = 0.0;
        std::array<double, 8> mean_rates;
        std::fill(mean_rates.begin(), mean_rates.end(), 0.0);
        for (auto&& tetra : m_list_tetrahedra) {
            for (int idx_band = 4; idx_band < nb_bands; ++idx_band) {
                double dos_tetra = tetra.compute_tetra_dos_energy_band(energy, idx_band);

                if (std::isnan(dos_tetra) || std::isinf(dos_tetra)) {
                    std::cout << "Energy: " << energy << " Band: " << idx_band << " DOS tetra: " << dos_tetra << std::endl;
                    throw std::runtime_error("DOS tetra is NaN or Inf");
                }
                dos += dos_tetra;
                std::array<double, 8> tetra_mean_rates = tetra.get_mean_electron_phonon_rates(idx_band);
                for (std::size_t idx_rate = 0; idx_rate < tetra_mean_rates.size(); ++idx_rate) {
                    mean_rates[idx_rate] += tetra_mean_rates[idx_rate] * dos_tetra;
                }
            }
        }
        energies.push_back(energy);
        list_mean_rates.push_back(mean_rates);
        list_dos.push_back(dos);
        total_dos += dos;
    }
    for (std::size_t idx_energy = 0; idx_energy < energies.size(); ++idx_energy) {
        file << energies[idx_energy] << " " << list_dos[idx_energy] << " ";
        for (auto&& rate : list_mean_rates[idx_energy]) {
            rate /= total_dos;
            file << rate << " ";
        }
        file << std::endl;
    }
    file.close();
}

void ElectronPhonon::plot_phonon_dispersion(const std::string& filename) const {
    std::ofstream file(filename);
    for (auto&& vtx : m_list_vertices) {
        auto k = vtx.get_position();
        file << k.x() << " " << k.y() << " " << k.z() << " ";
        for (auto&& mode_direction : m_phonon_dispersion) {
            auto q = k;
            q /= m_material.get_fourier_factor();
            auto e_ph = mode_direction.second.get_phonon_dispersion(q.norm()) * EmpiricalPseudopotential::Constants::h_bar_eV;
            file << e_ph << " ";
        }
        file << std::endl;
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
            if (carrierType == "electron") {
                if (mode == PhononMode::acoustic) {
                    m_acoustic_deformation_potential_e = deformationPotential;
                } else {
                    m_optical_deformation_potential_e = deformationPotential;
                }
            } else if (carrierType == "hole") {
                if (mode == PhononMode::acoustic) {
                    m_acoustic_deformation_potential_h = deformationPotential;
                } else {
                    m_optical_deformation_potential_h = deformationPotential;
                }
            } else {
                throw std::runtime_error("Unknown carrier type.");
            }
        }
    }
}

void ElectronPhonon::export_rate_values(const std::string& filename) const {
    std::ofstream file(filename);
    for (auto&& vertex : m_list_vertices) {
        std::vector<std::array<double, 8>> all_rates = vertex.get_electron_phonon_rates_all_bands();
        for (std::size_t idx_band = 0; idx_band < all_rates.size(); ++idx_band) {
            double energy = vertex.get_energy_at_band(idx_band);
            file << energy << " ";
            for (auto&& rate : all_rates[idx_band]) {
                file << rate << " ";
            }
            file << std::endl;
        }
    }
    file.close();
}

}  // namespace bz_mesh