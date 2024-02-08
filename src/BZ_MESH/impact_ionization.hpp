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

#include "Material.h"
#include "bz_mesh.hpp"
#include "bz_states.hpp"

namespace bz_mesh {

typedef std::complex<double> complex_d;

class ImpactIonization : private BZ_States {
 private:
    double m_impact_ionization_rate = 0.0;

 public:
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