#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>

#include "materials.hpp"
#include "pbmc_material_model.hpp"
#include "pbmc_scattering_model.hpp"

TEST_CASE("PBMC silicon model loads all parameter families from the material repository") {
    const uepm::physics::material_repository repository;
    auto                                     common_material = repository.load_material("Si");
    common_material.mass_density_kg_m3                       = 2400.0;
    common_material.static_relative_permittivity             = 12.25;

    const auto material = uepm::PBMC::load_pbmc_material_model(repository, common_material);

    CHECK(material.m_dielectric.epsilon_r == doctest::Approx(12.25));
    CHECK(material.m_electron_acoustic.mass_density_kg_per_m3 == doctest::Approx(2400.0));
    CHECK(material.m_electron_acoustic.sound_velocity_m_per_s == doctest::Approx(6606.666666666667));

    REQUIRE(material.m_electron_valleys.size() == 6);
    CHECK(material.m_electron_valleys.front().name() == "Delta_x_plus");
    CHECK(material.m_electron_valleys.front().transverse_effective_mass() ==
          doctest::Approx(0.1905 * uepm::constants::m_e));

    REQUIRE(material.m_hole_bands.size() == 2);
    CHECK(material.m_hole_bands.back().name() == "light_hole");

    REQUIRE(material.m_electron_intervalley_transitions.size() == 5);
    CHECK(material.m_electron_intervalley_transitions[2].m_name == "g3_LO");
    CHECK(material.m_electron_intervalley_transitions[2].m_deformation_potential_0 == doctest::Approx(3.4e10));

    REQUIRE(material.m_hole_optical_transitions.size() == 4);
    CHECK(material.m_hole_optical_transitions[1].name == "hh_to_lh");
    CHECK(material.m_hole_optical_transitions[1].deformation_potential_eV_per_m == doctest::Approx(8.0e10));

    CHECK(material.m_impurity_mobility.m_electron.m_mu0_cm2_per_V_s == doctest::Approx(1417.0));
    CHECK(material.m_impact_ionization.m_hole.m_prefactor_s_1 == doctest::Approx(1.4e12));
}

TEST_CASE("repository-backed PBMC parameters preserve silicon scattering rates") {
    const auto material = uepm::PBMC::load_pbmc_material_model();

    const double acoustic_rate = uepm::PBMC::acoustic_scattering_rate(material.m_electron_valleys.front(),
                                                                      material.m_electron_acoustic,
                                                                      0.1,
                                                                      300.0);
    CHECK(acoustic_rate == doctest::Approx(3.694230372508554e12).epsilon(1.0e-12));

    const double intervalley_rate =
        uepm::PBMC::intervalley_scattering_rate(material.m_electron_valleys.front(),
                                                material.m_electron_intervalley_transitions[2],
                                                material.m_electron_acoustic.mass_density_kg_per_m3,
                                                0.1,
                                                true,
                                                300.0);
    CHECK(intervalley_rate == doctest::Approx(1.500525823504873e11).epsilon(1.0e-12));
}
