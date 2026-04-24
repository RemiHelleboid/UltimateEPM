/**
 * @file poisson_solver_3d.hpp
 * @author Rémi Helleboid (remi.helleboid@st.com)
 * @brief
 * @version 0.1
 * @date 2021-10-19
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "finite_element3d.hpp"
#include "materials.hpp"

namespace uepm::fem {

class poisson_solver_3d : public FiniteElementP1System3d {
 private:
    std::vector<std::shared_ptr<mesh::element>> m_list_bulk_elements;
    std::vector<double>                         m_pre_computed_second_member;
    physic::material::list_materials            m_list_materials;

 public:
    poisson_solver_3d(mesh::mesh* p_mesh, std::size_t system_size, const physic::material::list_materials& mat_list)
        : FiniteElementP1System3d(p_mesh, system_size),
          m_list_materials{mat_list} {}
    void            compute_stiffness_matrix() override;
    Eigen::Vector4d compute_charge_density_elementary_second_member(const std::shared_ptr<mesh::element>& sp_tetra);
    void            update_second_member();
};

}  // namespace uepm::fem
