#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>

#include "materials.hpp"

TEST_CASE("material database loads typed common values") {
    const uepm::physics::material_repository repository;
    const auto                               materials = repository.load_all_materials();
    const auto&                              silicon   = materials.require(uepm::physics::material_id::silicon);

    CHECK(silicon.name == "Silicon");
    CHECK(silicon.symbol == "Si");
    CHECK(silicon.lattice_constant_m == doctest::Approx(5.431e-10));
    CHECK(silicon.mass_density_kg_m3 == doctest::Approx(2329.0));
    CHECK(silicon.static_relative_permittivity == doctest::Approx(11.7));
    CHECK(&silicon == &materials.require("Silicon"));
    CHECK(&silicon == &materials.require("Si"));
}

TEST_CASE("material repository discovers named parameter sets") {
    const uepm::physics::material_repository repository;

    CHECK(repository.material_file("Si").filename() == "material.yaml");
    CHECK(repository.has_parameter_set("Si", "epm", "local-cohen"));
    CHECK(repository.has_parameter_set("Si", "electron_phonon", "kamakura"));
    CHECK(repository.has_parameter_set("Si", "admc", "arora-canali"));
    CHECK(repository.has_parameter_set("Si", "pbmc", "default"));
    CHECK_FALSE(repository.has_parameter_set("Si", "epm", "unknown"));
}

TEST_CASE("missing material lookup is explicit") {
    uepm::physics::material_database materials;
    materials.add(uepm::physics::silicon_material_info());

    CHECK_THROWS_WITH_AS(materials.require("Ge"),
                         "Material 'Ge' is not available in material database.",
                         std::runtime_error);
}
