#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest/doctest.h"
#include "impact_ionization_matrix_element.hpp"

namespace {

uepm::mesh_bz::ImpactIonizationPlaneWaveState make_state(double kx, double energy_eV, int active_basis_index) {
    uepm::mesh_bz::ImpactIonizationPlaneWaveState state;
    state.k_SI                             = uepm::mesh_bz::vector3(kx, 0.0, 0.0);
    state.energy_eV                        = energy_eV;
    state.coefficients                     = Eigen::VectorXcd::Zero(2);
    state.coefficients[active_basis_index] = 1.0;
    return state;
}

}  // namespace

TEST_CASE("screened two-body matrix element enforces momentum conservation") {
    const std::vector<uepm::mesh_bz::vector3> basis = {
        uepm::mesh_bz::vector3(0.0, 0.0, 0.0),
        uepm::mesh_bz::vector3(10.0, 0.0, 0.0),
    };
    const auto interaction = [](const uepm::mesh_bz::vector3& q, double energy_eV) {
        CHECK(q.x() == doctest::Approx(1.0));
        CHECK(q.y() == doctest::Approx(0.0));
        CHECK(q.z() == doctest::Approx(0.0));
        CHECK(energy_eV == doctest::Approx(2.0));
        return uepm::mesh_bz::complex_d{4.0, 2.0};
    };
    const uepm::mesh_bz::ImpactIonizationMatrixElementConfig config{
        .momentum_tolerance_SI   = 1.0e-6,
        .normalization_volume_m3 = 2.0,
    };

    const auto initial_a = make_state(1.0, 5.0, 0);
    const auto final_a   = make_state(0.0, 3.0, 0);
    const auto initial_b = make_state(-1.0, -1.0, 0);
    const auto final_b   = make_state(0.0, 1.0, 0);

    const auto matrix_element = uepm::mesh_bz::compute_screened_two_body_matrix_element(initial_a,
                                                                                        initial_b,
                                                                                        final_a,
                                                                                        final_b,
                                                                                        basis,
                                                                                        interaction,
                                                                                        config);
    CHECK(matrix_element.real() == doctest::Approx(2.0));
    CHECK(matrix_element.imag() == doctest::Approx(1.0));

    const auto momentum_forbidden = uepm::mesh_bz::compute_screened_two_body_matrix_element(initial_a,
                                                                                            initial_b,
                                                                                            final_a,
                                                                                            make_state(0.5, 1.0, 0),
                                                                                            basis,
                                                                                            interaction,
                                                                                            config);
    CHECK(momentum_forbidden.real() == doctest::Approx(0.0));
    CHECK(momentum_forbidden.imag() == doctest::Approx(0.0));
}

TEST_CASE("direct exchange matrix element swaps final electrons") {
    const std::vector<uepm::mesh_bz::vector3> basis = {
        uepm::mesh_bz::vector3(0.0, 0.0, 0.0),
        uepm::mesh_bz::vector3(10.0, 0.0, 0.0),
    };
    const auto interaction = [](const uepm::mesh_bz::vector3& q, double) {
        return uepm::mesh_bz::complex_d{q.x(), 0.0};
    };
    const uepm::mesh_bz::ImpactIonizationMatrixElementConfig config{
        .momentum_tolerance_SI   = 1.0e-6,
        .normalization_volume_m3 = 1.0,
    };

    const auto initial_hot     = make_state(3.0, 5.0, 0);
    const auto initial_valence = make_state(-1.0, -1.0, 0);
    const auto final_1         = make_state(0.0, 3.0, 0);
    const auto final_2         = make_state(2.0, 1.0, 0);

    const auto [direct, exchange] =
        uepm::mesh_bz::compute_direct_exchange_impact_ionization_matrix_element(initial_hot,
                                                                                initial_valence,
                                                                                final_1,
                                                                                final_2,
                                                                                basis,
                                                                                interaction,
                                                                                config);

    CHECK(direct.real() == doctest::Approx(3.0));
    CHECK(exchange.real() == doctest::Approx(1.0));
    CHECK(uepm::mesh_bz::antisymmetrized_impact_ionization_strength(direct, exchange) == doctest::Approx(4.0));
}

TEST_CASE("screened two-body matrix element validates inputs") {
    const std::vector<uepm::mesh_bz::vector3> basis = {uepm::mesh_bz::vector3(0.0, 0.0, 0.0)};
    const auto interaction = [](const uepm::mesh_bz::vector3&, double) { return uepm::mesh_bz::complex_d{1.0, 0.0}; };
    auto       state       = make_state(0.0, 0.0, 0);
    state.coefficients     = Eigen::VectorXcd::Ones(2);

    CHECK_THROWS_AS(
        uepm::mesh_bz::compute_screened_two_body_matrix_element(state, state, state, state, basis, interaction),
        std::invalid_argument);

    const uepm::mesh_bz::ImpactIonizationMatrixElementConfig bad_config{
        .momentum_tolerance_SI   = 0.0,
        .normalization_volume_m3 = 1.0,
    };
    auto valid_state         = make_state(0.0, 0.0, 0);
    valid_state.coefficients = Eigen::VectorXcd::Ones(1);
    CHECK_THROWS_AS(uepm::mesh_bz::compute_screened_two_body_matrix_element(valid_state,
                                                                            valid_state,
                                                                            valid_state,
                                                                            valid_state,
                                                                            basis,
                                                                            interaction,
                                                                            bad_config),
                    std::invalid_argument);
}
