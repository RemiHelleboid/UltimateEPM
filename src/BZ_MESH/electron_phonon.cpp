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
#include <filesystem>

#include "Constants.hpp"
#include "Vector3D.h"
#include "bz_states.hpp"
#include "gmsh.h"
#include "omp.h"

namespace bz_mesh {

inline double ElectronPhonon::bose_einstein_distribution(double energy, double temperature) {
    double x = energy / (EmpiricalPseudopotential::Constants::k_b_eV * temperature);
    return 1.0 / std::expm1(x);
}

double ElectronPhonon::electron_overlap_integral(const vector3& k1, const vector3& k2) {
    constexpr double R_Wigner_Seitz = 2.122e-10;
    const double     qRws           = (k1 - k2).norm() * R_Wigner_Seitz;
    double           integral       = 3.0 * (std::sin(qRws) - qRws * std::cos(qRws)) / (qRws * qRws * qRws);
    return integral;
}

double ElectronPhonon::hole_overlap_integral(int n1, const vector3& k1, int n2, const vector3& k2) {
    const double cos_angle_k1_k2   = compte_cos_angle(k1, k2);
    const double cos_angle_k1_k2_2 = cos_angle_k1_k2 * cos_angle_k1_k2;
    auto         A_B_params        = m_hole_overlap_int_params.get_params(n1, n2);
    double       integral          = (1.0 / 2.0) * std::sqrt(A_B_params[0] + A_B_params[1] * cos_angle_k1_k2_2);
    return integral;
}

RateValues ElectronPhonon::compute_electron_phonon_rate(int idx_n1, std::size_t idx_k1) {
    RateValues                rates_k1_n1;
    const std::vector<Tetra>& list_tetrahedra          = m_list_tetrahedra;
    const auto&               indices_conduction_bands = m_indices_conduction_bands;

    // Initial state energy (eV) and k (in SI units 1/m)
    const double energy_n1_k1 = m_list_vertices[idx_k1].get_energy_at_band(idx_n1);
    const auto   k1           = m_list_vertices[idx_k1].get_position();

    for (int idx_n2 : indices_conduction_bands) {
        for (const auto& tetra : list_tetrahedra) {
            // Final k chosen as tetra barycenter (SI)
            const vector3 k2 = tetra.compute_barycenter();

            // Overlap (dimensionless)
            // const double overlap2 = std::pow(electron_overlap_integral(k1, k2), 2);
            const double overlap  = electron_overlap_integral(k1, k2);
            const double overlap2 = overlap * overlap;

            // Phonon wavevector q = k2 - k1 (SI)
            vector3 q_ph = k2 - k1;

            // Fold q back to first BZ if needed (no Umklapp yet)
            if (!is_inside_mesh_geometry(q_ph)) {
                try {
                    // Fold q to the first BZ
                    q_ph = fold_ws_bcc(q_ph);
                } catch (const std::runtime_error& e) {
                    std::cerr << "Error folding q: for index " << idx_k1 << "\n";
                    continue;  // Skip this tetrahedron if folding fails
                }
            }
            const double q_ph_norm = q_ph.norm();
            // if (!is_inside_mesh_geometry(q_ph)) throw std::runtime_error("Q is not inside the BZ");
            assert(is_inside_mesh_geometry(q_ph) && "Q is not inside the BZ");

            // DBG
            // continue;
            // Loop phonon branches
            for (const auto& ph_mode : m_phonon_dispersion) {
                const PhononModeDirection mode_direction = ph_mode.first;
                const auto&               disp           = ph_mode.second;

                // ---- Phonon quantities
                // omega [1/s] from dispersion(|q|)
                const double omega = disp.get_phonon_dispersion_from_lookup(q_ph_norm);
                // if (omega <= 0.0) continue;  // skip unphysical/zero frequency points
                assert(omega >= 0.0 && "Phonon frequency is negative");

                // E_ph in eV for Bose factor
                const double Eph_eV = EmpiricalPseudopotential::Constants::h_bar_eV * omega;

                // Basic sanity
                // if (Eph_eV > 0.1)  // 100 meV cap as you had
                //     throw std::runtime_error("Energy phonon too high");
                assert(Eph_eV < 0.1 && "Phonon energy is too high");

                // ---- Loop emission/absorption
                for (double sign_phonon : {-1.0, 1.0}) {
                    // Bose factor in eV units
                    const double N0        = bose_einstein_distribution(Eph_eV, m_temperature);  // dimensionless
                    const double bose_part = (sign_phonon < 0.0) ? (N0 + 1.0) : N0;              // +1: emission, plain N0: absorption

                    // Final electronic energy (eV)
                    const double E_final_eV = energy_n1_k1 + sign_phonon * Eph_eV;
                    if (!tetra.is_energy_inside_band(E_final_eV, idx_n2)) {
                        continue;
                    }

                    const double dos_eV    = tetra.compute_tetra_dos_energy_band(E_final_eV, idx_n2);
                    const double dos_per_J = dos_eV / EmpiricalPseudopotential::Constants::q_e;

                    // Deformation potential in SI Joules
                    const DeformationPotential defpot = (mode_direction.first == PhononMode::acoustic) ? m_acoustic_deformation_potential_e
                                                                                                       : m_optical_deformation_potential_e;
                    const double               Delta_J =
                        defpot.get_deformation_potential(q_ph, E_final_eV) * EmpiricalPseudopotential::Constants::q_e;  // eV -> J

                    // P = (pi / (rho * omega)) * Delta^2 * |I|^2 * (bose) * DOS(E)
                    double rate_value = (EmpiricalPseudopotential::Constants::pi / (m_rho * omega)) * (Delta_J * Delta_J) * overlap2 *
                                        bose_part * dos_per_J;
                    rate_value /= m_reduce_bz_factor;  // Correct for BZ volume if mesh does not match theoretical BZ volume
                    rate_value *= m_spin_degeneracy;   // Account for spin degeneracy of conduction bands

                    // rates_k1_n1.add_rate(RateValue(phonon_mode, phonon_direction, phonon_event, rate_value));
                    rates_k1_n1.add_rate(RateValue(mode_direction.first,
                                                   mode_direction.second,
                                                   (sign_phonon < 0.0) ? PhononEvent::emission : PhononEvent::absorption,
                                                   rate_value));
                }
            }
        }
    }

    return rates_k1_n1;
}

RateValues ElectronPhonon::compute_hole_phonon_rate(int idx_n1, std::size_t idx_k1) {
    RateValues  rates_k1_n1;
    const auto& list_tetrahedra       = m_list_tetrahedra;
    const auto& indices_valence_bands = m_indices_valence_bands;

    // Initial state energy (eV) and k (SI: 1/m)
    const double  energy_n1_k1 = m_list_vertices[idx_k1].get_energy_at_band(idx_n1);
    const vector3 k1           = m_list_vertices[idx_k1].get_position();

    for (int idx_n2 : indices_valence_bands) {
        for (const auto& tetra : list_tetrahedra) {
            // Final k chosen as tetra barycenter (SI)
            const vector3 k2 = tetra.compute_barycenter();

            // Overlap (dimensionless)
            // Keep your hole overlap signature; if you also have hole_overlap_integral(k1,k2),
            // you can swap it in, but this call matches your old function.
            const double overlap  = hole_overlap_integral(idx_n1, k1, idx_n2, k2);
            const double overlap2 = overlap * overlap;

            // Phonon wavevector q = k2 - k1 (SI)
            vector3 q_ph = k2 - k1;

            // Fold q back to first BZ if needed (Normal processes only, same as electrons)
            if (!is_inside_mesh_geometry(q_ph)) {
                try {
                    q_ph = fold_ws_bcc(q_ph);
                } catch (const std::runtime_error&) {
                    // Skip this tetrahedron if folding fails
                    continue;
                }
            }
            const double q_ph_norm = q_ph.norm();
            assert(is_inside_mesh_geometry(q_ph) && "Q is not inside the BZ");

            // Loop phonon branches
            for (const auto& ph_mode : m_phonon_dispersion) {
                const PhononModeDirection mode_direction = ph_mode.first;
                const auto&               disp           = ph_mode.second;

                // ω(|q|) [1/s]
                const double omega = disp.get_phonon_dispersion(q_ph_norm);
                if (omega <= 0.0) continue;  // guard ω→0 to avoid division by zero

                // ħω in eV
                const double Eph_eV = EmpiricalPseudopotential::Constants::h_bar_eV * omega;

                // Bose factor (dimensionless)
                const double N0 = bose_einstein_distribution(Eph_eV, m_temperature);

                // Hole deformation potential (hole-specific)
                const DeformationPotential defpot =
                    (mode_direction.first == PhononMode::acoustic) ? m_acoustic_deformation_potential_h : m_optical_deformation_potential_h;

                // Emission / absorption branches — mapping matches electrons:
                // sign = -1 → emission (Ef = Ei - ħω, bose = N0+1)
                // sign = +1 → absorption (Ef = Ei + ħω, bose = N0)
                for (double sign_phonon : {-1.0, 1.0}) {
                    const double bose_part  = (sign_phonon < 0.0) ? (N0 + 1.0) : N0;
                    const double E_final_eV = energy_n1_k1 + sign_phonon * Eph_eV;
                    if (!tetra.is_energy_inside_band(E_final_eV, idx_n2)) continue;

                    // DOS(E) and unit conversion same pattern as electrons
                    const double dos_eV    = tetra.compute_tetra_dos_energy_band(E_final_eV, idx_n2);
                    const double dos_per_J = dos_eV / EmpiricalPseudopotential::Constants::q_e;

                    // Δ in Joules (eV → J) at (q, E_final)
                    const double Delta_J = defpot.get_deformation_potential(q_ph, E_final_eV) * EmpiricalPseudopotential::Constants::q_e;

                    // Prefactor preserved to match electron routine (π / (ρ ω)) × Δ² × overlap² × bose × DOS(J⁻¹)
                    double rate_value = (EmpiricalPseudopotential::Constants::pi / (m_rho * omega)) * (Delta_J * Delta_J) * overlap2 *
                                        bose_part * dos_per_J;
                    rate_value /= m_reduce_bz_factor;  // Correct for BZ volume if mesh does not match theoretical BZ volume
                    rate_value *= m_spin_degeneracy; // Account for spin degeneracy of valence bands

                    const PhononEvent event = (sign_phonon < 0.0) ? PhononEvent::emission : PhononEvent::absorption;

                    rates_k1_n1.add_rate(RateValue(mode_direction.first, mode_direction.second, event, rate_value));
                }
            }
        }
    }

    return rates_k1_n1;
}

void ElectronPhonon::compute_electron_phonon_rates_over_mesh() {
    auto indices_conduction_bands = m_indices_conduction_bands;
    auto min_idx_conduction_band  = *std::min_element(indices_conduction_bands.begin(), indices_conduction_bands.end());
    auto max_idx_conduction_band  = *std::max_element(indices_conduction_bands.begin(), indices_conduction_bands.end());
    std::cout << "Min index conduction band: " << min_idx_conduction_band << std::endl;
    std::cout << "Computing electron-phonon rates over mesh for " << m_list_vertices.size() << " k-points." << std::endl;
    
    
    // Counter for progress display
    std::cout << "Progress: 0%";
    std::atomic<std::size_t> counter{0};

#pragma omp parallel for schedule(static)
    for (std::size_t idx_k1 = 0; idx_k1 < m_list_vertices.size(); ++idx_k1) {
        auto done = ++counter;
        if ((done % 10) == 0 && omp_get_thread_num() == 0) {
            std::cout << "\rDone " << done << "/" << m_list_vertices.size() << " (" << std::fixed << std::setprecision(1)
                      << (100.0 * done / m_list_vertices.size()) << "%)" << std::flush;
        }
        for (std::size_t idx_n1 = 0; idx_n1 < min_idx_conduction_band; ++idx_n1) {
            auto hole_rate = compute_hole_phonon_rate(idx_n1, idx_k1);
            auto array_h   = hole_rate.to_array();
            m_list_vertices[idx_k1].add_electron_phonon_rates(array_h);
        }
        for (std::size_t idx_n1 = min_idx_conduction_band; idx_n1 <= max_idx_conduction_band; ++idx_n1) {
            auto rate  = compute_electron_phonon_rate(idx_n1, idx_k1);
            auto array = rate.to_array();
            m_list_vertices[idx_k1].add_electron_phonon_rates(array);
        }
    }
    std::cout << std::endl;
}

void ElectronPhonon::compute_plot_electron_phonon_rates_vs_energy_over_mesh(int                nb_bands,
                                                                            double             max_energy,
                                                                            double             energy_step,
                                                                            const std::string& filename) {
    if (energy_step <= 0.0) throw std::invalid_argument("energy_step must be > 0");
    if (max_energy < 0.0) throw std::invalid_argument("max_energy must be >= 0");

    // Require mesh present and discover total bands to clamp nb_bands
    if (m_list_vertices.empty()) throw std::runtime_error("No vertices in mesh.");
    const int total_bands = static_cast<int>(m_list_vertices.front().get_number_bands());
    if (total_bands <= 0) throw std::runtime_error("Mesh reports zero bands.");

    nb_bands = std::clamp(nb_bands, 0, total_bands);
    if (nb_bands == 0) throw std::runtime_error("nb_bands clamped to 0; nothing to process.");

    // Build output path with suffix `_8ch` before extension
    namespace fs = std::filesystem;

    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open " + filename + " for writing.");

    // Column labels (match your IDX_* / to_array() order exactly)
    static constexpr std::array<const char*, 8>
        kLabels{"ac_L_ab", "ac_T_ab", "op_L_ab", "op_T_ab", "ac_L_em", "ac_T_em", "op_L_em", "op_T_em"};

    // Header
    out << "# E[eV]  DOS(E)";
    for (auto* s : kLabels)
        out << ' ' << s;
    out << '\n';

    // Integer-stepped energy sweep to avoid FP drift
    const std::size_t n_steps = static_cast<std::size_t>(std::floor(max_energy / energy_step)) + 1;

    for (std::size_t istep = 0; istep < n_steps; ++istep) {
        const double E = std::min(max_energy, istep * energy_step);
        std::cout << "\rEnergy: " << E << " / " << max_energy << std::flush;

        double                dos_sum = 0.0;
        std::array<double, 8> num{};  // DOS-weighted numerators for each channel
        num.fill(0.0);

        for (const auto& tetra : m_list_tetrahedra) {
            for (int b = 0; b < nb_bands; ++b) {
                const double dos_t = tetra.compute_tetra_dos_energy_band(E, b);
                if (!std::isfinite(dos_t))
                    throw std::runtime_error("DOS is NaN/Inf at E=" + std::to_string(E) + " band=" + std::to_string(b));
                if (dos_t <= 0.0) continue;

                dos_sum += dos_t;

                const std::array<double, 8> rates = tetra.get_mean_electron_phonon_rates(b);
                for (int i = 0; i < 8; ++i)
                    num[i] += rates[i] * dos_t;
            }
        }

        std::array<double, 8> mean{};
        if (dos_sum > 0.0) {
            const double inv_dos = 1.0 / dos_sum;
            for (int i = 0; i < 8; ++i)
                mean[i] = num[i] * inv_dos;
        } else {
            mean.fill(0.0);
        }

        out << std::scientific << std::setprecision(10) << E << ' ' << dos_sum;
        for (double v : mean)
            out << ' ' << v;
        out << '\n';
    }

    std::cout << std::endl;
    out.close();
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

void ElectronPhonon::add_electron_phonon_rates_to_mesh(const std::string& initial_filename, const std::string& final_filename) {
    // If the file exists, remove it to avoid appending to an old file
    if (std::ifstream(final_filename)) {
        std::remove(final_filename.c_str());
    }

    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 0);
    gmsh::model::add("bz_mesh");
    gmsh::open(initial_filename);

    std::string model_file_name;
    gmsh::model::getCurrent(model_file_name);

    std::vector<std::size_t> node_tags;
    std::vector<double>      nodeCoords;
    std::vector<double>      nodeParams;
    gmsh::model::mesh::reclassifyNodes();
    gmsh::model::mesh::getNodes(node_tags, nodeCoords, nodeParams, -1, -1, false, false);

    auto indices_conduction_bands = m_indices_conduction_bands;
    auto min_idx_conduction_band  = *std::min_element(indices_conduction_bands.begin(), indices_conduction_bands.end());
    auto max_idx_conduction_band  = *std::max_element(indices_conduction_bands.begin(), indices_conduction_bands.end());
    std::cout << "Min index conduction band: " << min_idx_conduction_band << std::endl;
    int nb_bands = m_indices_conduction_bands.size() + m_indices_valence_bands.size();

    for (int idx_val_band = 0; idx_val_band < max_idx_conduction_band; ++idx_val_band) {
        std::vector<double> rates_ac_lo_em(m_list_vertices.size());
        std::vector<double> rates_ac_lo_ab(m_list_vertices.size());
        std::vector<double> rates_ac_tr_em(m_list_vertices.size());
        std::vector<double> rates_ac_tr_ab(m_list_vertices.size());
        std::vector<double> rates_opt_lo_em(m_list_vertices.size());
        std::vector<double> rates_opt_lo_ab(m_list_vertices.size());
        std::vector<double> rates_opt_tr_em(m_list_vertices.size());
        std::vector<double> rates_opt_tr_ab(m_list_vertices.size());

        for (std::size_t idx_k1 = 0; idx_k1 < m_list_vertices.size(); ++idx_k1) {
            auto rates               = m_list_vertices[idx_k1].get_electron_phonon_rates(idx_val_band);
            rates_ac_lo_em[idx_k1] = rates[0];
            rates_ac_lo_ab[idx_k1] = rates[1];
            rates_ac_tr_em[idx_k1] = rates[2];
            rates_ac_tr_ab[idx_k1] = rates[3];
            rates_opt_lo_em[idx_k1] = rates[4];
            rates_opt_lo_ab[idx_k1] = rates[5];
            rates_opt_tr_em[idx_k1] = rates[6];
            rates_opt_tr_ab[idx_k1] = rates[7];
        }
        std::string name_rate_ac_lo_em = "ac_lo_em_" + std::to_string(idx_val_band);
        std::string name_rate_ac_lo_ab = "ac_lo_ab_" + std::to_string(idx_val_band);
        std::string name_rate_ac_tr_em = "ac_tr_em_" + std::to_string(idx_val_band);
        std::string name_rate_ac_tr_ab = "ac_tr_ab_" + std::to_string(idx_val_band);
        std::string name_rate_opt_lo_em = "opt_lo_em_" + std::to_string(idx_val_band);
        std::string name_rate_opt_lo_ab = "opt_lo_ab_" + std::to_string(idx_val_band);
        std::string name_rate_opt_tr_em = "opt_tr_em_" + std::to_string(idx_val_band);
        std::string name_rate_opt_tr_ab = "opt_tr_ab_" + std::to_string(idx_val_band);

        int data_tag_ac_lo_em = gmsh::view::add(name_rate_ac_lo_em);
        int data_tag_ac_lo_ab = gmsh::view::add(name_rate_ac_lo_ab);
        int data_tag_ac_tr_em = gmsh::view::add(name_rate_ac_tr_em);
        int data_tag_ac_tr_ab = gmsh::view::add(name_rate_ac_tr_ab);
        int data_tag_opt_lo_em = gmsh::view::add(name_rate_opt_lo_em);
        int data_tag_opt_lo_ab = gmsh::view::add(name_rate_opt_lo_ab);
        int data_tag_opt_tr_em = gmsh::view::add(name_rate_opt_tr_em);
        int data_tag_opt_tr_ab = gmsh::view::add(name_rate_opt_tr_ab);

        gmsh::view::addHomogeneousModelData(data_tag_ac_lo_em, 0, model_file_name, "NodeData", node_tags, rates_ac_lo_em);
        gmsh::view::addHomogeneousModelData(data_tag_ac_lo_ab, 0, model_file_name, "NodeData", node_tags, rates_ac_lo_ab);
        gmsh::view::addHomogeneousModelData(data_tag_ac_tr_em, 0, model_file_name, "NodeData", node_tags, rates_ac_tr_em);
        gmsh::view::addHomogeneousModelData(data_tag_ac_tr_ab, 0, model_file_name, "NodeData", node_tags, rates_ac_tr_ab);
        gmsh::view::addHomogeneousModelData(data_tag_opt_lo_em, 0, model_file_name, "NodeData", node_tags, rates_opt_lo_em);
        gmsh::view::addHomogeneousModelData(data_tag_opt_lo_ab, 0, model_file_name, "NodeData", node_tags, rates_opt_lo_ab);
        gmsh::view::addHomogeneousModelData(data_tag_opt_tr_em, 0, model_file_name, "NodeData", node_tags, rates_opt_tr_em);
        gmsh::view::addHomogeneousModelData(data_tag_opt_tr_ab, 0, model_file_name, "NodeData", node_tags, rates_opt_tr_ab);

        gmsh::option::setNumber("PostProcessing.SaveMesh", 1);  // Save mesh only once
        gmsh::view::write(data_tag_ac_lo_em, final_filename, true);
        gmsh::option::setNumber("PostProcessing.SaveMesh", 0);
        gmsh::view::write(data_tag_ac_lo_ab, final_filename, true);
        gmsh::view::write(data_tag_ac_tr_em, final_filename, true);
        gmsh::view::write(data_tag_ac_tr_ab, final_filename, true);
        gmsh::view::write(data_tag_opt_lo_em, final_filename, true);
        gmsh::view::write(data_tag_opt_lo_ab, final_filename, true);
        gmsh::view::write(data_tag_opt_tr_em, final_filename, true);
        gmsh::view::write(data_tag_opt_tr_ab, final_filename, true);
    }
    gmsh::finalize();
    


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
            double q_max_norm = 1.5 / m_si2red;
            double n_points    = 200;
            std::cout << "Max q norm in reduced units: " << q_max_norm << std::endl;
            phononDispersion.fill_lookup_table(q_max_norm, n_points);
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

            PhononMode           mode;
            if (wave == "acoustic") {
                mode = PhononMode::acoustic;
            } else if (wave == "optic") {
                mode = PhononMode::optical;
            } else {
                throw std::runtime_error("Unknown wave type.");
            }
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
            std::cout << "Energy at band " << idx_band << ": " << energy << std::endl;
            file << idx_band << "," << energy << ",";
            for (std::size_t idx_rate = 0; idx_rate < all_rates[idx_band].size(); ++idx_rate) {
                double rate = all_rates[idx_band][idx_rate];
                file << rate << ((idx_rate < all_rates[idx_band].size() - 1) ? "," : "");
            }
            file << std::endl;
        }
    }
    file.close();
}

}  // namespace bz_mesh