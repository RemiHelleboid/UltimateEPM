/**
 * @file bz_dielectric.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-12-20
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <Eigen/Dense>

#include "Material.h"
#include "bz_mesh.hpp"

namespace uepm::mesh_bz {

class BZ_States : public MeshBZ {
 protected:

    std::vector<Vector3D<int>> m_basisVectors;

    std::vector<Eigen::VectorXd> m_eigenvalues_k;
    std::vector<Eigen::VectorXd> m_eigenvalues_k_plus_q;

    std::vector<Eigen::MatrixXcd> m_eigenvectors_k;
    std::vector<Eigen::MatrixXcd> m_eigenvectors_k_plus_q;

    Vector3D<double>    m_q_shift;
    std::vector<double> m_list_energies;

    /**
     * @brief Real part of the dielectric function.
     * m_dielectric_function_real[idx_energy] is the real part of the dielectric function at the energy m_energies[idx_energy].
     *
     */
    std::vector<double> m_dielectric_function_real;

    // m_vtx_dielectric_function_real[idx_vtx][idx_energy] is the real part of the dielectric function at the energy m_energies[idx_energy]
    // and at the vertex m_vertices[idx_vtx].
    std::vector<std::vector<double>> m_vtx_dielectric_function_real;

    /**
     * @brief The index of the first q-point this class is responsible for.
     * This is used to parallelize the computation of the dielectric function where
     * each instance of this class is responsible for a subset of the q-points of the mesh.
     *
     * For a calculation on a single CPU the offset is 0.
     *
     */
    std::size_t m_offset_q_index = 0;

    /**
     * @brief Number of q-points this class is responsible for.
     * This is used to parallelize the computation of the dielectric function where
     * each instance of this class is responsible for a subset of the q-points.
     *
     * For a calculation on a single CPU this parameter is equal to the size of the m_vertices.size().
     *
     */
    std::size_t m_nb_kpoints = 0;

 public:
    BZ_States(const uepm::pseudopotential::Material& material) : MeshBZ(material) {}
    BZ_States(const BZ_States& other) = default;

    void                              set_nb_bands(int nb_bands) { m_nb_bands = nb_bands; }
    void                              set_basis_vectors(const std::vector<Vector3D<int>>& basis_vectors) { m_basisVectors = basis_vectors; }
    const std::vector<Vector3D<int>>& get_basis_vectors() const { return m_basisVectors; }
    void                              compute_eigenstates(int nb_threads = 1);
    void                              compute_shifted_eigenstates(const Vector3D<double>& q_shift, int nb_threads = 1);

    const std::vector<double>& get_energies() const { return m_list_energies; }
    void                       set_energies(const std::vector<double>& energies) { m_list_energies = energies; }

    const std::vector<Eigen::MatrixXcd>& get_eigen_states() const { return m_eigenvectors_k; }

    const double compute_fermi_level(double doping_concentration, double temperature) const;

    void compute_dielectric_function(const std::vector<double>& energies, double eta_smearing, int nb_threads = 1);
    void export_dielectric_function(const std::string& prefix) const;

    void populate_vtx_dielectric_function(const std::vector<double>& energies, double eta_smearing);

    std::complex<double> get_dielectric_function(const vector3& q, double energy) const {
        // TODO: implement this function
        return 1.0;
    }

    void export_full_eigenstates() const;
};

}  // namespace uepm::mesh_bz