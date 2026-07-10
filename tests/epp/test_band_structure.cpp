#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "BandStructure.h"
#include "Vector3D.h"
#include "epm_material.hpp"
#include "materials.hpp"

namespace {

uepm::pseudopotential::epm_material load_si_material(const std::string& parameter_set) {
    const uepm::physics::material_repository repository;
    uepm::pseudopotential::Materials         materials;
    materials.load_material(repository, "Si", parameter_set);
    return materials.materials.at("Si");
}

std::vector<Vector3D<double>> tiny_k_path() {
    return {
        Vector3D<double>(0.0, 0.0, 0.0),
        Vector3D<double>(0.1, 0.0, 0.0),
    };
}

}  // namespace

TEST_CASE("serial band structure computes one gradient row per k-point") {
    const auto material = load_si_material("local-cohen");
    const auto kpoints  = tiny_k_path();

    uepm::pseudopotential::BandStructure band_structure;
    band_structure.Initialize(material, 2, kpoints, 2, false, false);
    band_structure.Compute(true);

    REQUIRE(band_structure.get_band_energies().size() == kpoints.size());
    REQUIRE(band_structure.get_band_energy_gradients().size() == kpoints.size());
    CHECK(band_structure.get_band_energy_gradients().at(0).size() == 2);
    CHECK(band_structure.get_band_energy_gradients().at(1).size() == 2);
}

TEST_CASE("parallel band structure stores gradients by k-point and honors thread-count overload") {
    const auto material = load_si_material("local-cohen");
    const auto kpoints  = tiny_k_path();

    uepm::pseudopotential::BandStructure band_structure;
    band_structure.Initialize(material, 2, kpoints, 2, false, false);
    band_structure.Compute_parallel(true, 2);

    REQUIRE(band_structure.get_band_energies().size() == kpoints.size());
    REQUIRE(band_structure.get_band_energy_gradients().size() == kpoints.size());
    CHECK(band_structure.get_band_energy_gradients().at(0).size() == 2);
    CHECK(band_structure.get_band_energy_gradients().at(1).size() == 2);

    band_structure.Compute_parallel(2);
    CHECK(band_structure.get_band_energies().size() == kpoints.size());
    CHECK(band_structure.get_band_energy_gradients().empty());
}

TEST_CASE("band structure rejects unsupported gradient modes and invalid access") {
    const auto material = load_si_material("potz-vogl");
    const auto kpoints  = tiny_k_path();

    uepm::pseudopotential::BandStructure band_structure;
    CHECK_THROWS_AS(band_structure.get_band(0), std::runtime_error);

    band_structure.Initialize(material, 2, kpoints, 2, true, false);
    CHECK_THROWS_AS(band_structure.Compute(true), std::runtime_error);
    CHECK_THROWS_AS(band_structure.Compute_parallel(true, 1), std::runtime_error);

    band_structure.Initialize(material, 2, kpoints, 2, false, false);
    band_structure.Compute_parallel(1);
    CHECK_THROWS_AS(band_structure.get_band(2), std::out_of_range);
}
