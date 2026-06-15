#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>

#include "epm_material.hpp"
#include "materials.hpp"

TEST_CASE("EPP parameters reuse the common material identity and lattice constant") {
    const uepm::physics::material_repository repository;
    uepm::pseudopotential::Materials         epm_materials;
    epm_materials.load_material(repository, "Si", "potz-vogl");

    const auto& silicon = epm_materials.materials.at("Si");
    CHECK(silicon.get_id() == uepm::physics::material_id::silicon);
    CHECK(silicon.get_material_info().symbol == "Si");
    CHECK(silicon.get_material_info().name == "Silicon");
    CHECK(silicon.get_lattice_constant_meter() == doctest::Approx(repository.load_material("Si").lattice_constant_m));
    CHECK(silicon.get_lattice_constant_meter() == doctest::Approx(5.431e-10));
}
