/**
 * @file reciprocal_space.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-18
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <Eigen/Core>
#include <vector>

#include "vector_bz.hpp"

namespace uepm::mesh_bz {

class ReciprocalSpace {
 private:
    Eigen::Matrix3d      m_basis         = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d      m_inverse_basis = Eigen::Matrix3d::Identity();
    double               m_halfwidth     = 1.0;
    std::vector<vector3> m_reciprocal_shifts;

 public:
    void initialize_basis(const Eigen::Vector3d& b1_SI,
                          const Eigen::Vector3d& b2_SI,
                          const Eigen::Vector3d& b3_SI,
                          double                 halfwidth_reduced);

    void precompute_bcc_shifts(double si_to_reduced, int maximum_shell = 5);

    bool    contains_bcc(const vector3& k_SI, double si_to_reduced) const noexcept;
    vector3 retrieve_bcc_image(const vector3& k_SI, double si_to_reduced) const;
    vector3 ReciprocalSpace::fold_bcc_fast_SI(const vector3& k_SI, double si_to_reduced) const noexcept;
    vector3 fold_wigner_seitz(const vector3& k_SI) const noexcept;
    bool    inside_wigner_seitz(const vector3& k_SI, double si_to_reduced) const noexcept;
    bool    inside_irreducible_wedge(const vector3& k_SI, double si_to_reduced) const noexcept;
};

}  // namespace uepm::mesh_bz
