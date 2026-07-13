#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest/doctest.h"
#include "nonlinear_poisson_solver.hpp"

#include <Eigen/Sparse>

namespace {

Eigen::SparseMatrix<double> two_node_stiffness() {
    Eigen::SparseMatrix<double> matrix(2, 2);
    matrix.insert(0, 0) = 1.0;
    matrix.insert(0, 1) = -1.0;
    matrix.insert(1, 0) = -1.0;
    matrix.insert(1, 1) = 1.0;
    matrix.makeCompressed();
    return matrix;
}

}  // namespace

TEST_CASE("Nonlinear Poisson preserves a neutral uniform equilibrium") {
    const auto stiffness = two_node_stiffness();
    const Eigen::Vector2d fixed_source = Eigen::Vector2d::Zero();
    const Eigen::Vector2d volume(0.5, 0.5);
    const Eigen::Vector2d electrons(2.0, 2.0);
    const Eigen::Vector2d holes(2.0, 2.0);
    const Eigen::Vector2d reference = Eigen::Vector2d::Constant(0.2);
    const Eigen::Vector2d thermal_voltage = Eigen::Vector2d::Constant(0.02585);
    const std::vector<uepm::fem::nonlinear_poisson_dirichlet> boundary{{0, 0.2}};

    const auto result = uepm::fem::solve_nonlinear_poisson_lumped(
        stiffness, fixed_source, volume, electrons, holes, reference, reference, thermal_voltage, boundary);

    CHECK(result.converged);
    CHECK(result.iterations == 0);
    CHECK(result.potential_V(0) == doctest::Approx(0.2));
    CHECK(result.potential_V(1) == doctest::Approx(0.2));
}

TEST_CASE("Nonlinear Poisson screens a perturbed mobile charge") {
    const auto stiffness = two_node_stiffness();
    const Eigen::Vector2d fixed_source(0.0, 0.75);
    const Eigen::Vector2d volume(0.5, 0.5);
    const Eigen::Vector2d electrons(0.0, 1.0);
    const Eigen::Vector2d holes = Eigen::Vector2d::Zero();
    const Eigen::Vector2d reference = Eigen::Vector2d::Zero();
    const Eigen::Vector2d initial = Eigen::Vector2d::Zero();
    const Eigen::Vector2d thermal_voltage = Eigen::Vector2d::Constant(0.02585);
    const std::vector<uepm::fem::nonlinear_poisson_dirichlet> boundary{{0, 0.0}};

    uepm::fem::nonlinear_poisson_options options;
    options.potential_tolerance_V = 1.0e-10;
    options.relative_residual_tolerance = 1.0e-10;
    const auto result = uepm::fem::solve_nonlinear_poisson_lumped(
        stiffness, fixed_source, volume, electrons, holes, reference, initial, thermal_voltage, boundary, options);

    CHECK(result.converged);
    CHECK(result.iterations > 0);
    CHECK(result.final_residual_norm < 1.0e-9);
    CHECK(result.potential_V(0) == doctest::Approx(0.0));
    CHECK(result.potential_V(1) > 0.0);
    CHECK(result.potential_V(1) < 0.1);
}

TEST_CASE("Nonlinear Poisson rejects negative carrier density") {
    const auto stiffness = two_node_stiffness();
    const Eigen::Vector2d zero = Eigen::Vector2d::Zero();
    const Eigen::Vector2d volume = Eigen::Vector2d::Ones();
    const Eigen::Vector2d electrons(0.0, -1.0);
    const Eigen::Vector2d thermal_voltage = Eigen::Vector2d::Constant(0.02585);

    CHECK_THROWS_AS(uepm::fem::solve_nonlinear_poisson_lumped(
                        stiffness, zero, volume, electrons, zero, zero, zero, thermal_voltage, {}),
                    std::invalid_argument);
}
