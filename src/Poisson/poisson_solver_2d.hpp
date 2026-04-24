/**
 * @file poisson_solver_2d.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-06
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "finite_element2d.hpp"
#include "materials.hpp"

namespace uepm::fem {

class poisson_solver_2d : public FiniteElementP1System2d {
 private:
    std::vector<std::shared_ptr<mesh::element>> m_list_bulk_elements;
    std::vector<double>                         m_pre_computed_second_member;
    std::string                                 name_electron_charge_density_function = "eDensity";
    std::string                                 name_hole_charge_density_function     = "hDensity";
    physic::material::list_materials            m_list_materials;

 public:
    poisson_solver_2d(mesh::mesh* p_mesh, std::size_t system_size, const physic::material::list_materials& list_materials)
        : FiniteElementP1System2d(p_mesh, system_size),
          m_list_materials(list_materials){};

    void                   compute_stiffness_matrix() override;
    static Eigen::Vector3d compute_charge_density_elementary_second_member(const std::shared_ptr<mesh::element>& sp_triangle);
    void                   update_second_member();
};

}  // namespace uepm::fem
