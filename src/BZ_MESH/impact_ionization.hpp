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
#include "bz_states.hpp"

namespace bz_mesh {

class ImpactIonization : private BZ_States {
 private:
 public:
    double compute_direct_impact_ionization_matrix_element(int                     idx_n1,
                                                           int                     idx_n1_prime,
                                                           int                     idx_n2,
                                                           int                     idx_n2_prime,
                                                           const Vector3D<double>& k1,
                                                           const Vector3D<double>& k1_prime,
                                                           const Vector3D<double>& k2,
                                                           const Vector3D<double>& k2_prime) const;

    void compute_impact_ionization_rate();
};


}  // namespace bz_mesh