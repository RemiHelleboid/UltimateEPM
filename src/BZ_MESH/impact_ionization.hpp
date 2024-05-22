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
#include <array>
#include <memory>

#include "Material.h"
#include "bz_mesh.hpp"
#include "bz_states.hpp"
#include "dielectric_mesh.hpp"

namespace bz_mesh {

typedef std::complex<double>       complex_d;
typedef std::unique_ptr<BZ_States> uptr_BZstates;

class ImpactIonization {
 private:
    /**
     * @brief List of Brillouin Zone with their precomputed electronic states.
     * Each of them correspond to a 1BZ shifted by a vector G of the reciprocal lattice.
     *
     */
    std::vector<uptr_BZstates> m_list_BZ_states;

    /**
     * @brief Path to the initial mesh file (centered at 0,0,0)
     *
     */
    std::string m_initial_mesh_path;

    /**
     * @brief Maximum radius of the Brillouin Zone in the reciprocal space.
     *
     */
    double m_max_radius_G0_BZ = 0.0;

    /**
     * @brief Material of the system.
     *
     */
    EmpiricalPseudopotential::Material m_material;

    /**
     * @brief Mesh of the dielectric function.
     *
     */
    DielectricMesh m_dielectric_mesh;

    /**
     * @brief List of the results of the impact ionization rate.
     *
     */
    std::vector<double> m_impact_ionization_results;

 public:
    ImpactIonization(const EmpiricalPseudopotential::Material& material, const std::string& initial_mesh_path);
    void read_dielectric_file(const std::string& filename);
    void interp_test_dielectric_function(std::string filename);

    double get_max_radius_G0_BZ() const { return m_max_radius_G0_BZ; }
    void   set_max_radius_G0_BZ(double max_radius_G0_BZ) { m_max_radius_G0_BZ = max_radius_G0_BZ; }
    void   compute_eigenstates(int nb_threads = 1);

    std::array<complex_d, 2> compute_direct_indirect_impact_ionization_matrix_element(int idx_n1,
                                                                                      int idx_n1_prime,
                                                                                      int idx_n2,
                                                                                      int idx_n2_prime,
                                                                                      int idx_k1,
                                                                                      int idx_k1_prime,
                                                                                      int idx_k2,
                                                                                      int idx_k2_prime) const;

    double compute_impact_ionization_rate(int idx_n1, const Vector3D<double>& k1);
};

}  // namespace bz_mesh