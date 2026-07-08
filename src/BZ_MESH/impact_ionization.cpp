/**
 * @file bz_states.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-12-20
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "impact_ionization.hpp"

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>

#include "BandStructure.h"
#include "Hamiltonian.h"
#include "Options.h"
#include "bz_mesh.hpp"
#include "bz_states.hpp"
#include "epm_material.hpp"
#include "impact_ionization_matrix_element.hpp"
#include "omp.h"
#include "physical_constants.hpp"

namespace uepm::mesh_bz {

ImpactIonization::ImpactIonization(const uepm::pseudopotential::epm_material& material,
                                   const std::string&                         initial_mesh_path)
    : m_material(material),
      m_dielectric_mesh(material) {
    std::filesystem::path path(initial_mesh_path);
    if (!std::filesystem::exists(path)) {
        throw std::invalid_argument("Impact-ionization mesh file does not exist: " + initial_mesh_path);
    }
    std::ifstream file(initial_mesh_path);
    if (!file.is_open()) {
        throw std::invalid_argument("Could not open impact-ionization mesh file: " + initial_mesh_path);
    }
    m_initial_mesh_path = initial_mesh_path;
}

void ImpactIonization::read_dielectric_file(const std::string& filename) {
    m_dielectric_mesh.read_mesh_geometry_from_msh_file(filename, true);
    m_dielectric_mesh.build_search_tree();
    m_dielectric_mesh.read_dielectric_file(filename);
}

void ImpactIonization::interp_test_dielectric_function(std::string filename) {
    double        eps = 1e-6;
    double        x0  = eps;
    double        y0  = eps;
    double        z0  = eps;
    double        x1  = 3.0;
    double        y1  = 3.0;
    double        z1  = 3.0;
    int           nx  = 20;
    int           ny  = 20;
    int           nz  = 20;
    double        dx  = (x1 - x0) / nx;
    double        dy  = (y1 - y0) / ny;
    double        dz  = (z1 - z0) / nz;
    std::ofstream file(filename);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                double    x = x0 + i * dx;
                double    y = y0 + j * dy;
                double    z = z0 + k * dz;
                vector3   position(x, y, z);
                complex_d epsilon =
                    m_dielectric_mesh.interpolate_dielectric_function(m_dielectric_mesh.reduced_to_si_k(position),
                                                                      0.0102);
                file << x << ", " << y << ", " << z << ", " << epsilon.real() << ", " << epsilon.imag() << std::endl;
                std::cout << "Position: " << position << " epsilon: " << epsilon << std::endl;
            }
        }
    }
    file.close();
}

ImpactIonizationStateTable ImpactIonization::make_state_table(std::size_t bz_state_index) const {
    if (bz_state_index >= m_list_BZ_states.size()) {
        throw std::out_of_range("ImpactIonization::make_state_table: invalid BZ state index.");
    }
    return ImpactIonizationStateTable(*m_list_BZ_states[bz_state_index]);
}

void ImpactIonization::compute_eigenstates(int nb_threads) {
    int                      nb_bands_to_use = 16;
    uepm::mesh_bz::BZ_States my_bz_mesh(m_material);
    my_bz_mesh.set_nb_bands(nb_bands_to_use);
    uepm::pseudopotential::BandStructure band_structure{};
    int                                  nb_nearest_neighbors = 10;
    bool                                 nonlocal_epm         = false;
    bool                                 enable_soc           = false;
    band_structure.Initialize(m_material, nb_bands_to_use, {}, nb_nearest_neighbors, nonlocal_epm, enable_soc);
    auto basis = band_structure.get_basis_vectors();
    my_bz_mesh.set_basis_vectors(basis);

    my_bz_mesh.read_mesh_geometry_from_msh_file(m_initial_mesh_path);

    const vector3 b1 = {-1.0, 1.0, 1.0};
    const vector3 b2 = {1.0, -1.0, 1.0};
    const vector3 b3 = {1.0, 1.0, -1.0};

    // std::cout << "k: " << k << std::endl;
    double factor = 2.0 * M_PI / m_material.get_lattice_constant_meter();
    // test k + G
    std::vector<int> list_n_k = {0, 1, -1, 2, -2, 3, -3, 4, -4};
    // std::vector<int> list_n_k = {0, 1, -1, 2, -2};

    std::cout << "Max radius G0 BZ: " << m_max_radius_G0_BZ << std::endl;
    for (auto&& n_k_x : list_n_k) {
        for (auto&& n_k_y : list_n_k) {
            for (auto&& n_k_z : list_n_k) {
                vector3 G_BZ = n_k_x * b1 + n_k_y * b2 + n_k_z * b3;
                if (G_BZ.norm() > m_max_radius_G0_BZ) {
                    continue;
                }
                std::cout << "G_BZ: " << G_BZ << " --> " << G_BZ.norm() << std::endl;
                G_BZ = G_BZ * factor;

                auto ptr_BZ_states = std::make_unique<BZ_States>(m_material);
                ptr_BZ_states->set_nb_bands(nb_bands_to_use);
                ptr_BZ_states->set_basis_vectors(basis);
                ptr_BZ_states->read_mesh_geometry_from_msh_file(m_initial_mesh_path);
                ptr_BZ_states->shift_bz_center(G_BZ);
                ptr_BZ_states->compute_eigenstates(nb_threads);
                m_list_BZ_states.push_back(std::move(ptr_BZ_states));
            }
        }
    }
    std::cout << "Number of BZ states: " << m_list_BZ_states.size() << std::endl;
}

std::array<complex_d, 2> ImpactIonization::compute_direct_indirect_impact_ionization_matrix_element(
    int idx_n1,
    int idx_n1_prime,
    int idx_n2,
    int idx_n2_prime,
    int idx_k1,
    int idx_k1_prime,
    int idx_k2,
    int idx_k2_prime) const {
    if (m_list_BZ_states.empty()) {
        throw std::logic_error("ImpactIonization matrix element requires computed eigenstates.");
    }
    if (idx_n1 < 0 || idx_n1_prime < 0 || idx_n2 < 0 || idx_n2_prime < 0 || idx_k1 < 0 || idx_k1_prime < 0 ||
        idx_k2 < 0 || idx_k2_prime < 0) {
        throw std::invalid_argument("ImpactIonization matrix element indices must be non-negative.");
    }

    const auto state_table = make_state_table(0);
    const auto make_state = [&state_table](int idx_k, int idx_band) {
        const auto metadata = state_table.state(static_cast<std::size_t>(idx_k), static_cast<std::size_t>(idx_band));
        return ImpactIonizationPlaneWaveState{.k_SI         = metadata.k_SI,
                                              .energy_eV    = metadata.energy_eV,
                                              .coefficients = state_table.coefficients(metadata.k_index,
                                                                                      metadata.band_index)};
    };

    std::vector<vector3> basis_vectors_SI;
    const double         reciprocal_factor = m_material.get_fourier_factor();
    basis_vectors_SI.reserve(m_list_BZ_states[0]->get_basis_vectors().size());
    for (const auto& basis_vector : m_list_BZ_states[0]->get_basis_vectors()) {
        basis_vectors_SI.emplace_back(static_cast<double>(basis_vector.X) * reciprocal_factor,
                                      static_cast<double>(basis_vector.Y) * reciprocal_factor,
                                      static_cast<double>(basis_vector.Z) * reciprocal_factor);
    }

    const auto screened_coulomb = make_screened_coulomb_kernel();
    const auto interaction = [&screened_coulomb](const vector3& q_SI, double energy_transfer_eV) {
        return screened_coulomb.interaction_eV_m3(q_SI, energy_transfer_eV);
    };
    const ImpactIonizationMatrixElementConfig config{
        .momentum_tolerance_SI  = 1.0e6,
        .normalization_volume_m3 = m_material.get_atomic_volume(),
    };

    return compute_direct_exchange_impact_ionization_matrix_element(make_state(idx_k1, idx_n1),
                                                                    make_state(idx_k2, idx_n2),
                                                                    make_state(idx_k1_prime, idx_n1_prime),
                                                                    make_state(idx_k2_prime, idx_n2_prime),
                                                                    basis_vectors_SI,
                                                                    interaction,
                                                                    config);
}

double ImpactIonization::compute_impact_ionization_rate(int idx_n1, std::size_t idx_k1) {
    (void)idx_n1;
    (void)idx_k1;
    throw std::logic_error("ImpactIonization::compute_impact_ionization_rate is experimental and incomplete; no "
                           "validated rate is available.");

#if 0
    constexpr int                     nb_valence_bands    = 3;
    constexpr int                     nb_conduction_bands = 4;
    constexpr int                     min_conduction_band = 4;
    const std::size_t                 nb_vtx              = m_list_BZ_states[0]->get_list_vertices().size();
    const std::vector<Vector3D<int>>& basis_vector_G      = m_list_BZ_states[0]->get_basis_vectors();

    std::cout << "Start computing impact ionization rate for band " << idx_n1 << " and k-point " << idx_k1 << std::endl;
    auto start_precompute = std::chrono::high_resolution_clock::now();
    std::cout << "Nb of Bz states: " << m_list_BZ_states.size() << std::endl;
    std::cout << "Nb cols: " << m_list_BZ_states[0]->get_eigen_states()[0].cols() << std::endl;

    // Sum_2_prime[idx_band][idx_node]
    std::vector<std::vector<complex_d>> Sum_2_prime(nb_conduction_bands);
    for (int idx_n2_prime = 0; idx_n2_prime < nb_conduction_bands; ++idx_n2_prime) {
        int n2_prime = idx_n2_prime + min_conduction_band;
        Sum_2_prime[idx_n2_prime].resize(nb_vtx);
        for (std::size_t idx_node = 0; idx_node < nb_vtx; ++idx_node) {
            const Eigen::MatrixXcd& A_2_prime   = m_list_BZ_states[0]->get_eigen_states()[idx_node];
            Sum_2_prime[idx_n2_prime][idx_node] = A_2_prime.col(n2_prime).sum();
        }
    }
    std::cout << "Sum_2_prime done" << std::endl;
    // Sum_2[idx_band][idx_node]
    std::vector<std::vector<complex_d>> Sum_2(nb_valence_bands);
    for (int idx_n2 = 0; idx_n2 < nb_valence_bands; ++idx_n2) {
        Sum_2[idx_n2].resize(nb_vtx);
        for (std::size_t idx_node = 0; idx_node < nb_vtx; ++idx_node) {
            const Eigen::MatrixXcd& A_2 = m_list_BZ_states[0]->get_eigen_states()[idx_node];
            Sum_2[idx_n2][idx_node]     = A_2.col(idx_n2).sum();
        }
    }

    std::cout << "Sum_2 done" << std::endl;

    // Sum_2_prime[idx_band][idx_node] (n1, k1 are in an outter loop)
    std::vector<std::vector<complex_d>> Sum_1_prime_1(nb_conduction_bands);
    for (int idx_n1_prime = 0; idx_n1_prime < nb_conduction_bands; ++idx_n1_prime) {
        std::size_t n1_prime = idx_n1_prime + min_conduction_band;
        Sum_1_prime_1[idx_n1_prime].resize(nb_vtx);
        for (std::size_t idx_k1_prime = 0; idx_k1_prime < nb_vtx; ++idx_k1_prime) {
            const Eigen::MatrixXcd& A_1_prime = m_list_BZ_states[0]->get_eigen_states()[idx_k1_prime];
            const Eigen::MatrixXcd& A_1       = m_list_BZ_states[0]->get_eigen_states()[idx_k1];
            int                     nb_Gvect  = A_1_prime.rows();
            complex_d               sum       = 0.0;
            for (int idx_G1_prime = 0; idx_G1_prime < nb_Gvect; ++idx_G1_prime) {
                for (int idx_G1 = 0; idx_G1 < nb_Gvect; ++idx_G1) {
                    auto    G1_prime = basis_vector_G[idx_G1_prime];
                    auto    GG1      = basis_vector_G[idx_G1];
                    vector3 GA       = vector3(G1_prime.X + GG1.X, G1_prime.Y + GG1.Y, G1_prime.Z + GG1.Z);
                    vector3 GA_prime = vector3(G1_prime.X, G1_prime.Y, G1_prime.Z);
                    auto    k1       = m_list_BZ_states[0]->get_vertex_position(idx_k1);
                    auto    k1_prime = m_list_BZ_states[0]->get_vertex_position(idx_k1_prime);
                    // auto    q_a      = k1 - k1_prime + G1 + (-1 * G1_prime);
                    double energy_w =
                        m_list_BZ_states[0]->get_energies()[idx_n1] - m_list_BZ_states[0]->get_energies()[n1_prime];
                    // complex_d epsilon    = m_dielectric_mesh.interpolate_dielectric_function(q_a, energy_w);
                    // complex_d epsilon    = 1.0;
                    // complex_d factor_eps = uepm::constants::q_e * uepm::constants::q_e /
                    //                        (uepm::constants::eps_0 * epsilon * q_a.norm() * q_a.norm());
                    // sum += std::conj(A_1_prime(idx_G1_prime, n1_prime)) * A_1(idx_G1, idx_n1) * factor_eps;
                }
            }
            Sum_1_prime_1[idx_n1_prime][idx_k1_prime] = sum;
        }
    }

    auto end_precompute = std::chrono::high_resolution_clock::now();
    std::cout << "Precompute time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_precompute - start_precompute).count()
              << " ms" << std::endl;
    std::cout << "DONE PRECOMPUTE\n Start computing matrix element" << std::endl;

    auto start_compute = std::chrono::high_resolution_clock::now();

    const std::vector<Vertex>& list_vertices = m_list_BZ_states[0]->get_list_vertices();

    for (int idx_n2 = 0; idx_n2 < nb_valence_bands; ++idx_n2) {
        for (int idx_n1_prime = 0; idx_n1_prime < nb_conduction_bands; ++idx_n1_prime) {
            std::size_t n1_prime = idx_n1_prime + min_conduction_band;
            for (int idx_n2_prime = 0; idx_n2_prime < nb_conduction_bands; ++idx_n2_prime) {
                std::size_t n2_prime = idx_n2_prime + min_conduction_band;
                for (std::size_t idx_k1_prime = 0; idx_k1_prime < nb_vtx; ++idx_k1_prime) {
                    std::cout << "idx_n1: " << idx_n1 << " idx_n1_prime: " << idx_n1_prime << " idx_n2: " << idx_n2
                              << " idx_n2_prime: " << idx_n2_prime << " idx_k1: " << idx_k1
                              << " idx_k1_prime: " << idx_k1_prime << std::endl;
                    std::vector<double> FullMatrixElement(nb_vtx);
                    for (std::size_t idx_k2_prime = 0; idx_k2_prime < nb_vtx; ++idx_k2_prime) {
                        // std::cout << "\r" << idx_k1_prime << " / " << nb_vtx << " --> " << idx_k2_prime << " / " <<
                        // nb_vtx << std::flush;
                        vector3 k_2_momentum = list_vertices[idx_k1].get_position() -
                                               list_vertices[idx_k1_prime].get_position() -
                                               list_vertices[idx_k2_prime].get_position();
                        if (k_2_momentum.norm() > m_max_radius_G0_BZ || k_2_momentum.norm() < 1e-12) {
                            continue;
                        }
                        std::size_t idx_k2 = 18;
                        complex_d   Ma     = Sum_2_prime[idx_n2_prime][idx_k2_prime] *
                                       Sum_1_prime_1[idx_n1_prime][idx_k1_prime] * Sum_2[idx_n2][idx_k2];
                        complex_d Mb = Sum_2_prime[idx_n1_prime][idx_k1_prime] *
                                       Sum_1_prime_1[idx_n2_prime][idx_k2_prime] * Sum_2[idx_n2][idx_k2];
                        FullMatrixElement[idx_k2_prime] = std::abs(Ma) * std::abs(Ma) + std::abs(Mb) * std::abs(Mb) +
                                                          std::abs(Ma - Mb) * std::abs(Ma - Mb);
                    }
                }
            }
        }
    }
    std::cout << std::endl;
    auto end_compute = std::chrono::high_resolution_clock::now();
    std::cout << "Compute time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_compute - start_compute).count() << " ms"
              << std::endl;
    return 0.0;
#endif
}

}  // namespace uepm::mesh_bz
