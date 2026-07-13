#include "nonlinear_poisson_solver.hpp"

#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace uepm::fem {
namespace {

void require_size(const Eigen::VectorXd& vector, Eigen::Index expected, const char* name) {
    if (vector.size() != expected) {
        throw std::invalid_argument(std::string("Nonlinear Poisson size mismatch for ") + name + ".");
    }
    if (!vector.allFinite()) {
        throw std::invalid_argument(std::string("Nonlinear Poisson received non-finite ") + name + ".");
    }
}

void validate_inputs(const Eigen::SparseMatrix<double>& stiffness,
                     const Eigen::VectorXd&             fixed_source,
                     const Eigen::VectorXd&             lumped_volume,
                     const Eigen::VectorXd&             electron_density,
                     const Eigen::VectorXd&             hole_density,
                     const Eigen::VectorXd&             reference_potential,
                     const Eigen::VectorXd&             initial_potential,
                     const Eigen::VectorXd&             thermal_voltage,
                     const std::vector<nonlinear_poisson_dirichlet>& dirichlet,
                     const nonlinear_poisson_options& options) {
    if (stiffness.rows() == 0 || stiffness.rows() != stiffness.cols()) {
        throw std::invalid_argument("Nonlinear Poisson stiffness matrix must be non-empty and square.");
    }
    const Eigen::Index size = stiffness.rows();
    require_size(fixed_source, size, "fixed source");
    require_size(lumped_volume, size, "lumped volume");
    require_size(electron_density, size, "electron density");
    require_size(hole_density, size, "hole density");
    require_size(reference_potential, size, "reference potential");
    require_size(initial_potential, size, "initial potential");
    require_size(thermal_voltage, size, "thermal voltage");
    if (!stiffness.coeffs().allFinite()) {
        throw std::invalid_argument("Nonlinear Poisson received a non-finite stiffness matrix.");
    }
    if ((lumped_volume.array() < 0.0).any() || (electron_density.array() < 0.0).any() ||
        (hole_density.array() < 0.0).any() || (thermal_voltage.array() <= 0.0).any()) {
        throw std::invalid_argument("Nonlinear Poisson requires non-negative volumes/densities and positive thermal voltage.");
    }
    if (options.max_iterations == 0 || !(options.relative_residual_tolerance > 0.0) ||
        !(options.potential_tolerance_V > 0.0) || !(options.maximum_update_V > 0.0) ||
        !(options.minimum_damping > 0.0 && options.minimum_damping <= 1.0) ||
        !(options.armijo_coefficient > 0.0 && options.armijo_coefficient < 1.0)) {
        throw std::invalid_argument("Invalid nonlinear Poisson options.");
    }
    for (const auto& condition : dirichlet) {
        if (condition.index < 0 || condition.index >= size || !std::isfinite(condition.value_V)) {
            throw std::invalid_argument("Invalid nonlinear Poisson Dirichlet condition.");
        }
    }
}

struct state {
    Eigen::VectorXd residual;
    Eigen::VectorXd electron_density;
    Eigen::VectorXd hole_density;
};

state evaluate(const Eigen::SparseMatrix<double>& stiffness,
               const Eigen::VectorXd& fixed_source,
               const Eigen::VectorXd& lumped_volume,
               const Eigen::VectorXd& electron_reference,
               const Eigen::VectorXd& hole_reference,
               const Eigen::VectorXd& reference_potential,
               const Eigen::VectorXd& thermal_voltage,
               const std::vector<nonlinear_poisson_dirichlet>& dirichlet,
               const Eigen::VectorXd& potential) {
    const Eigen::ArrayXd eta = (potential - reference_potential).array() / thermal_voltage.array();
    // Updates are limited before acceptance. Reaching this broad guard means
    // the nonlinear iterate is already far outside its useful physical range.
    if ((eta.abs() > 500.0).any()) {
        throw std::runtime_error("Nonlinear Poisson exponential argument exceeded its safety range.");
    }

    state value;
    value.electron_density = electron_reference.array() * eta.exp();
    value.hole_density     = hole_reference.array() * (-eta).exp();
    const Eigen::VectorXd mobile_source =
        (lumped_volume.array() * (value.hole_density - value.electron_density).array()).matrix();
    value.residual = stiffness * potential - fixed_source - mobile_source;
    for (const auto& condition : dirichlet) {
        value.residual(condition.index) = potential(condition.index) - condition.value_V;
    }
    if (!value.residual.allFinite() || !value.electron_density.allFinite() || !value.hole_density.allFinite()) {
        throw std::runtime_error("Nonlinear Poisson residual evaluation produced a non-finite value.");
    }
    return value;
}

Eigen::SparseMatrix<double> jacobian(const Eigen::SparseMatrix<double>& stiffness,
                                     const Eigen::VectorXd& lumped_volume,
                                     const Eigen::VectorXd& electron_density,
                                     const Eigen::VectorXd& hole_density,
                                     const Eigen::VectorXd& thermal_voltage,
                                     const std::vector<nonlinear_poisson_dirichlet>& dirichlet) {
    Eigen::SparseMatrix<double> result = stiffness;
    const Eigen::VectorXd diagonal = lumped_volume.array() *
                                     (electron_density + hole_density).array() / thermal_voltage.array();
    for (Eigen::Index i = 0; i < result.rows(); ++i) {
        result.coeffRef(i, i) += diagonal(i);
    }
    std::vector<bool> is_dirichlet(static_cast<std::size_t>(result.rows()), false);
    for (const auto& condition : dirichlet) {
        is_dirichlet[static_cast<std::size_t>(condition.index)] = true;
    }
    std::vector<Eigen::Triplet<double>> entries;
    entries.reserve(static_cast<std::size_t>(result.nonZeros()) + dirichlet.size());
    for (Eigen::Index outer = 0; outer < result.outerSize(); ++outer) {
        for (Eigen::SparseMatrix<double>::InnerIterator entry(result, outer); entry; ++entry) {
            if (!is_dirichlet[static_cast<std::size_t>(entry.row())]) {
                entries.emplace_back(entry.row(), entry.col(), entry.value());
            }
        }
    }
    for (const auto& condition : dirichlet) {
        entries.emplace_back(condition.index, condition.index, 1.0);
    }
    Eigen::SparseMatrix<double> constrained(result.rows(), result.cols());
    constrained.setFromTriplets(entries.begin(), entries.end());
    constrained.makeCompressed();
    return constrained;
}

}  // namespace

nonlinear_poisson_result solve_nonlinear_poisson_lumped(
    const Eigen::SparseMatrix<double>& stiffness,
    const Eigen::VectorXd& fixed_source,
    const Eigen::VectorXd& lumped_volume,
    const Eigen::VectorXd& electron_density_reference,
    const Eigen::VectorXd& hole_density_reference,
    const Eigen::VectorXd& reference_potential_V,
    const Eigen::VectorXd& initial_potential_V,
    const Eigen::VectorXd& thermal_voltage_V,
    const std::vector<nonlinear_poisson_dirichlet>& dirichlet_conditions,
    const nonlinear_poisson_options& options) {
    validate_inputs(stiffness, fixed_source, lumped_volume, electron_density_reference, hole_density_reference,
                    reference_potential_V, initial_potential_V, thermal_voltage_V, dirichlet_conditions, options);

    nonlinear_poisson_result result;
    result.potential_V = initial_potential_V;
    for (const auto& condition : dirichlet_conditions) {
        result.potential_V(condition.index) = condition.value_V;
    }

    state current = evaluate(stiffness, fixed_source, lumped_volume, electron_density_reference,
                             hole_density_reference, reference_potential_V, thermal_voltage_V,
                             dirichlet_conditions, result.potential_V);
    const double initial_norm = current.residual.norm();
    const double residual_scale = std::max(initial_norm, std::numeric_limits<double>::min());
    result.final_residual_norm = initial_norm;
    if (initial_norm == 0.0) {
        result.converged = true;
        return result;
    }

    for (std::size_t iteration = 0; iteration < options.max_iterations; ++iteration) {
        Eigen::SparseMatrix<double> matrix = jacobian(stiffness, lumped_volume, current.electron_density,
                                                      current.hole_density, thermal_voltage_V, dirichlet_conditions);
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
        solver.analyzePattern(matrix);
        solver.factorize(matrix);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Nonlinear Poisson Jacobian factorization failed.");
        }
        Eigen::VectorXd correction = solver.solve(-current.residual);
        if (solver.info() != Eigen::Success || !correction.allFinite()) {
            throw std::runtime_error("Nonlinear Poisson Newton solve failed.");
        }

        const double raw_maximum = correction.cwiseAbs().maxCoeff();
        double damping = raw_maximum > options.maximum_update_V
                             ? options.maximum_update_V / raw_maximum
                             : 1.0;
        const double current_merit = 0.5 * current.residual.squaredNorm();
        state trial;
        bool accepted = false;
        while (damping >= options.minimum_damping) {
            Eigen::VectorXd trial_potential = result.potential_V + damping * correction;
            trial = evaluate(stiffness, fixed_source, lumped_volume, electron_density_reference,
                             hole_density_reference, reference_potential_V, thermal_voltage_V,
                             dirichlet_conditions, trial_potential);
            const double trial_merit = 0.5 * trial.residual.squaredNorm();
            if (trial_merit <= (1.0 - options.armijo_coefficient * damping) * current_merit) {
                result.potential_V = std::move(trial_potential);
                accepted = true;
                break;
            }
            damping *= 0.5;
        }
        if (!accepted) {
            throw std::runtime_error("Nonlinear Poisson line search failed to reduce the residual.");
        }

        current = std::move(trial);
        result.iterations = iteration + 1;
        result.maximum_correction_V = damping * raw_maximum;
        result.last_damping = damping;
        result.final_residual_norm = current.residual.norm();
        const bool residual_converged = result.final_residual_norm / residual_scale <=
                                        options.relative_residual_tolerance;
        const bool update_converged = result.maximum_correction_V <= options.potential_tolerance_V;
        if (residual_converged && update_converged) {
            result.converged = true;
            break;
        }
    }
    return result;
}

}  // namespace uepm::fem
