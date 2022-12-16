/**
 * @file DielectricFunction.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-11-24
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <atomic>
#include <vector>

#include "Material.h"
#include "SymmetryPoints.h"
#include "Vector3D.h"

namespace EmpiricalPseudopotential {

class DielectricFunction {
 protected:
    std::vector<Vector3D<int>>    m_basisVectors;
    std::vector<Vector3D<double>> m_kpoints;
    const Material&               m_material;
    const int                     m_nb_bands;
    bool                          m_nonlocal_epm = false;

    std::vector<Vector3D<double>> m_qpoints;
    std::vector<double>           m_energies;

    std::vector<Eigen::VectorXd>  m_eigenvalues_k;
    std::vector<Eigen::MatrixXcd> m_eigenvectors_k;

    std::string m_export_prefix = "dielectric_function";

    /**
     * @brief The index of the first k-point this class is responsible for.
     * This is used to parallelize the computation of the dielectric function where
     * each instance of this class is responsible for a subset of the k-points.
     *
     * For a calculation on a single CPU the offset is 0.
     *
     */
    std::size_t m_offset_k_index = 0;

    /**
     * @brief Number of k-points this class is responsible for.
     * This is used to parallelize the computation of the dielectric function where
     * each instance of this class is responsible for a subset of the k-points.
     *
     * For a calculation on a single CPU this parameter is equal to the size of the m_list_k_points.
     *
     */
    std::size_t m_nb_kpoints = 0;

    /**
     * @brief m_dielectric_function[idx_q][idx_energy] is the dielectric function
     * q = m_qpoints[idx_q]
     * energy = m_energies[idx_energy]
     *
     */
    std::vector<std::vector<double>> m_dielectric_function_real;

    /**
     * @brief m_dielectric_function_imag[idx_q][idx_energy] is the dielectric function
     * q = m_qpoints[idx_q]
     * energy = m_energies[idx_energy]
     *
     */
    std::vector<std::vector<double>> m_dielectric_function_imag;

 public:
    DielectricFunction(const Material& material, const std::vector<Vector3D<int>>& basisVectors, const int nb_bands);

    void set_non_local_epm(const bool new_value) {m_nonlocal_epm = new_value;}

    /**
     * @brief Randomly generate a list of k-points in the irreducible wedge of the first Brillouin zone.
     *
     * @param nb_points
     */
    void generate_k_points_random(std::size_t nb_points);

    /**
     * @brief Generate a list of k-points in the irreducible wedge of the first Brillouin zone, using the
     * Monkhorst-Pack algorithm.
     *
     * @param nb_points
     */
    void generate_k_points_grid(std::size_t Nx, std::size_t Ny, std::size_t Nz, double shift, bool irreducible_wedge);

    /** Get the list of k-points.
     * @return const std::vector<Vector3D<double>>&
     */
    const std::vector<Vector3D<double>>& get_kpoints() const { return m_kpoints; }

    /**
     * @brief Get the offset of the first k-point this class is responsible for.
     *
     * @param offset_k_index
     */
    std::size_t get_offset_k_index() const { return m_offset_k_index; }

    /**
     * @brief Get the nb kpoints the instance of this class is responsible for.
     *
     * @return std::size_t
     */
    std::size_t get_nb_kpoints() const { return m_nb_kpoints; }

    /**
     * @brief Set the offset k index.
     *
     * @param offset_k_index
     */
    void set_offset_k_index(std::size_t offset_k_index) { m_offset_k_index = offset_k_index; }

    /**
     * @brief Set the number of k-points this class is responsible for.
     *
     * @param nb_kpoints
     */
    void set_nb_kpoints(std::size_t nb_kpoints) { m_nb_kpoints = nb_kpoints; }

    /**
     * @brief Set the list of q-points for which the dielectric function will be computed.
     *
     * @param kpoints
     */
    void set_qpoints(const std::vector<Vector3D<double>>& qpoints) { m_qpoints = qpoints; }

    /**
     * @brief Set the list of energies for which the dielectric function will be computed.
     *
     * @param energies
     */
    void set_energies(const std::vector<double>& energies) { m_energies = energies; }

    void set_export_prefix(const std::string& prefix) { m_export_prefix = prefix; }

    /**
     * @brief Compute the dielectric function.
     *
     */
    void compute_dielectric_function(double eta_smearing = 1e-2);

    void clear_eigen_states() {
        m_eigenvalues_k.clear();
        m_eigenvectors_k.clear();
    }

    /**
     * @brief Get the dielectric function result.
     *
     * @return const std::vector<std::vector<double>>&
     */
    const std::vector<std::vector<double>>& get_dielectric_function() const { return m_dielectric_function_real; }

    const std::vector<double> get_flat_dielectric_function() const {
        std::vector<double> result;
        for (const auto& q : m_dielectric_function_real) {
            result.insert(result.end(), q.begin(), q.end());
        }
        return result;
    }

    /**
     * @brief Merge the results of multiple instances of this class.
     * This is used to parallelize the computation of the dielectric function where
     * each instance of this class is responsible for a subset of the k-points.
     *
     * @param dielectric_function_results
     * @param nb_kpoints_per_instance
     * @return std::vector<std::vector<double>>
     */
    static DielectricFunction merge_results(DielectricFunction                                  RootDielectricFunction,
                                            const std::vector<std::vector<std::vector<double>>> dielectric_function_results,
                                            std::vector<int>                                    nb_kpoints_per_instance);

    /**
     * @brief Apply Kramer's Kronig relations to the dielectric function to obtain the real part.
     *
     */
    void apply_kramers_kronig();

    /**
     * @brief Export the grid of k-points to a file.
     *
     * @param filename
     */
    void export_kpoints(const std::string& filename) const;

    /**
     * @brief Export the results of the computation of the dielectric function to a file.
     *
     * @param filename
     */
    void export_dielectric_function_at_q(const std::string& filename, std::size_t idx_q, bool name_auto) const;

    void export_dielectric_function(const std::string& filename, bool name_auto) const;
};

}  // namespace EmpiricalPseudopotential