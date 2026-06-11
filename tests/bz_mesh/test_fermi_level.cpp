#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest/doctest.h"
#include "fermi_level.hpp"

using uepm::mesh_bz::fermi::charge_neutrality_residual;

TEST_CASE("charge neutrality uses physical dopant signs") {
    CHECK(charge_neutrality_residual(10.0, 2.0, 8.0, 0.0) == doctest::Approx(0.0));
    CHECK(charge_neutrality_residual(2.0, 10.0, 0.0, 8.0) == doctest::Approx(0.0));
    CHECK(charge_neutrality_residual(10.0, 2.0, 0.0, 0.0) > 0.0);
    CHECK(charge_neutrality_residual(2.0, 10.0, 0.0, 0.0) < 0.0);
}
