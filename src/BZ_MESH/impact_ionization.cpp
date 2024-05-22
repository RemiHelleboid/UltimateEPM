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
#include "Constants.hpp"
#include "Hamiltonian.h"
#include "Material.h"
#include "Options.h"
#include "bz_mesh.hpp"
#include "bz_states.hpp"
#include "omp.h"

namespace bz_mesh {

ImpactIonization::ImpactIonization(const EmpiricalPseudopotential::Material& material, const std::string& initial_mesh_path) {
    std::filesystem::path path(initial_mesh_path);
    if (!std::filesystem::exists(path)) {
        std::cerr << "Error: file " << initial_mesh_path << " does not exist." << std::endl;
        exit(1);
    }
    std::ifstream file(initial_mesh_path);
    if (!file.is_open()) {
        std::cerr << "Error: could not open file " << initial_mesh_path << std::endl;
        exit(1);
    }
    m_material          = material;
    m_initial_mesh_path = initial_mesh_path;
}

void ImpactIonization::read_dielectric_file(const std::string& filename) {
    bool normalize_by_fourier_factor = false;
    m_dielectric_mesh.read_mesh_geometry_from_msh_file(filename, normalize_by_fourier_factor);
    m_dielectric_mesh.build_search_tree();
    m_dielectric_mesh.read_dielectric_file(filename);
}

void ImpactIonization::interp_test_dielectric_function(std::string filename) {
    double eps = 1e-6;
    double x0 = eps;
    double y0 = eps;
    double z0 = eps;
    double x1 = 3.0;
    double y1 = 3.0;
    double z1 = 3.0;
    int    nx = 20;
    int    ny = 20;
    int    nz = 20;
    double dx = (x1 - x0) / nx;
    double dy = (y1 - y0) / ny;
    double dz = (z1 - z0) / nz;
    std::ofstream file(filename);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                double x = x0 + i * dx;
                double y = y0 + j * dy;
                double z = z0 + k * dz;
                vector3     position(x, y, z);
                complex_d epsilon = m_dielectric_mesh.interpolate_dielectric_function(position, 0.0102);
                file << x << ", " << y << ", " << z << ", " << epsilon.real() << ", " << epsilon.imag() << std::endl;
                std::cout << "Position: " << position << " epsilon: " << epsilon << std::endl;
            }
        }
    }
    file.close();
}

void ImpactIonization::compute_eigenstates(int nb_threads) {
    int                nb_bands_to_use = 16;
    bz_mesh::BZ_States my_bz_mesh(m_material);
    my_bz_mesh.set_nb_bands(nb_bands_to_use);
    EmpiricalPseudopotential::BandStructure band_structure{};
    int                                     nb_nearest_neighbors = 10;
    bool                                    nonlocal_epm         = false;
    bool                                    enable_soc           = false;
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

double ImpactIonization::compute_impact_ionization_rate(int idx_n1, const Vector3D<double>& k1) {
    constexpr int nb_valence_bands = 3;
    constexpr int nb_conduction_bands = 4;

    std::size_t nb_vtx = m_list_BZ_states[0]->get_list_vertices().size();

    for (int n1_prime = 0; n1_prime < nb_valence_bands; ++n1_prime) {
        for (int n2 = 0; n2 < nb_conduction_bands; ++n2) {
            for (int n2_prime = 0; n2_prime < nb_conduction_bands; ++n2_prime) {
                for (int k1_prime = 0; k1_prime < nb_vtx; ++k1_prime) {
                    
                }
            }   
        }
    }
    




}

}  // namespace bz_mesh