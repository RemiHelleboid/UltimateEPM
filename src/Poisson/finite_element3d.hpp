/**
 * @file finite_element3d.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-06
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include <iostream>
#include <vector>

#include "finite_element.hpp"
#include "mesh.hpp"

namespace uepm::fem {

class FiniteElementP1System3d : public FiniteElementSystem {
 protected:
 public:
    static const Eigen::Matrix4d ElementaryMassMatrixRefElement;

    FiniteElementP1System3d() : FiniteElementSystem() {}
    explicit FiniteElementP1System3d(std::size_t system_size) : FiniteElementSystem(system_size) {}
    FiniteElementP1System3d(mesh::mesh* p_mesh, std::size_t system_size) : FiniteElementSystem(p_mesh, system_size) {}
    ~FiniteElementP1System3d() {}

    static Eigen::Matrix4d compute_elementary_stiffness_matrix(std::shared_ptr<mesh::element>);
    static Eigen::Matrix4d compute_elementary_mass_matrix(std::shared_ptr<mesh::element>);
    static Eigen::Vector4d compute_elementary_second_member(std::shared_ptr<mesh::element>, double constant_value);
    static Eigen::Vector4d compute_elementary_second_member(std::shared_ptr<mesh::element>,
                                                            std::function<double(double, double, double)> function);
    void                   compute_stiffness_matrix() override;
    void                   compute_mass_matrix() override;
    void                   compute_second_member(double constant_value) override;
    void                   compute_second_member(std::function<double(double, double, double)> function);
    void                   apply_dirichlet_condition(const std::string& region_name, const double boundary_value);
    void                   apply_dirichlet_condition_second_member(const std::string& region_name, const double boundary_value);
    void                   apply_neuman_condition(const std::string& region_name, const double boundary_value);
};

}  // namespace uepm::fem