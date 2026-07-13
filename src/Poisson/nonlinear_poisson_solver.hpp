#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cstddef>
#include <vector>

namespace uepm::fem {

struct nonlinear_poisson_options {
    std::size_t max_iterations           = 30;
    double      relative_residual_tolerance = 1.0e-8;
    double      potential_tolerance_V    = 1.0e-6;
    double      maximum_update_V         = 0.05;
    double      minimum_damping          = 1.0 / 1024.0;
    double      armijo_coefficient       = 1.0e-4;
};

struct nonlinear_poisson_dirichlet {
    Eigen::Index index = 0;
    double       value_V = 0.0;
};

struct nonlinear_poisson_result {
    Eigen::VectorXd potential_V;
    std::size_t     iterations             = 0;
    double          final_residual_norm     = 0.0;
    double          maximum_correction_V   = 0.0;
    double          last_damping            = 1.0;
    bool            converged               = false;
};

/**
 * Solve a mass-lumped nonlinear Poisson correction.
 *
 * The supplied linear problem uses the convention
 *
 *     stiffness * psi = fixed_source + lumped_volume * (p - n)
 *
 * and the mobile reference densities are corrected according to
 *
 *     n = n_reference * exp((psi - reference_potential) / thermal_voltage)
 *     p = p_reference * exp(-(psi - reference_potential) / thermal_voltage).
 *
 * All density and volume units are caller-defined, but their product must use
 * the same units as fixed_source. The stiffness matrix must use the matching
 * Poisson scaling.
 */
nonlinear_poisson_result solve_nonlinear_poisson_lumped(
    const Eigen::SparseMatrix<double>&             stiffness,
    const Eigen::VectorXd&                         fixed_source,
    const Eigen::VectorXd&                         lumped_volume,
    const Eigen::VectorXd&                         electron_density_reference,
    const Eigen::VectorXd&                         hole_density_reference,
    const Eigen::VectorXd&                         reference_potential_V,
    const Eigen::VectorXd&                         initial_potential_V,
    const Eigen::VectorXd&                         thermal_voltage_V,
    const std::vector<nonlinear_poisson_dirichlet>& dirichlet_conditions,
    const nonlinear_poisson_options&               options = {});

}  // namespace uepm::fem
