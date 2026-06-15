#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>

#include <string>

#include "epm_material.hpp"
#include "materials.hpp"

TEST_CASE("EPP parameters reuse the common material identity and lattice constant") {
    const std::string filename = PROJECT_SRC_DIR + std::string("/parameter_files/materials.yaml");

    uepm::physics::material_database common_materials(filename);
    uepm::pseudopotential::Materials epm_materials;
    epm_materials.load_material_parameters(filename, common_materials);

    const auto& silicon = epm_materials.materials.at("Si");
    CHECK(silicon.get_id() == uepm::physics::material_id::silicon);
    CHECK(silicon.get_material_info().symbol == "Si");
    CHECK(silicon.get_material_info().name == "Silicon");
    CHECK(silicon.get_lattice_constant_meter() == doctest::Approx(common_materials.require("Si").lattice_constant_m));
    CHECK(silicon.get_lattice_constant_meter() == doctest::Approx(5.43e-10));
}
