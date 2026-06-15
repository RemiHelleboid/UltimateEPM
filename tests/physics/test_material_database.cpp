#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>

#include <string>

#include "materials.hpp"

TEST_CASE("material database loads typed common values") {
    const std::string filename = PROJECT_SRC_DIR + std::string("/examples/materials/materials.yaml");

    uepm::physics::material_database materials(filename);
    const auto&                      silicon = materials.require(uepm::physics::material_id::silicon);

    CHECK(silicon.name == "Silicon");
    CHECK(silicon.symbol == "Si");
    CHECK(silicon.lattice_constant_m == doctest::Approx(5.43e-10));
    CHECK(silicon.mass_density_kg_m3 == doctest::Approx(2329.0));
    CHECK(silicon.static_relative_permittivity == doctest::Approx(11.7));
    CHECK(&silicon == &materials.require("Silicon"));
    CHECK(&silicon == &materials.require("Si"));
}

TEST_CASE("missing material lookup is explicit") {
    uepm::physics::material_database materials;
    materials.add(uepm::physics::silicon_material_info());

    CHECK_THROWS_WITH_AS(materials.require("Ge"),
                         "epm_material 'Ge' is not available in material database.",
                         std::runtime_error);
}
