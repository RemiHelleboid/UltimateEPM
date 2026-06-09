#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <algorithm>
#include <string_view>

#include "amc_material_model.hpp"
#include "amc_scattering_model.hpp"
#include "amc_transport_kernel.hpp"

TEST_CASE("silicon material model preserves scattering parameters") {
    const auto material = uepm::amc::make_silicon_amc_material_model();

    CHECK(material.m_name == "Silicon");
    CHECK(material.m_symbol == "Si");
    CHECK(material.m_dielectric.epsilon_r == doctest::Approx(11.7));
    CHECK(material.m_electron_acoustic.mass_density_kg_per_m3 == doctest::Approx(2.329e3));
    CHECK(material.m_electron_acoustic.sound_velocity_m_per_s == doctest::Approx(6.6e3));
    CHECK(material.m_electron_acoustic.deformation_potential_eV == doctest::Approx(9.0));
    CHECK(material.m_hole_acoustic.sound_velocity_m_per_s == doctest::Approx(6.6e3));
    CHECK(material.m_hole_acoustic.deformation_potential_eV == doctest::Approx(5.5));
    CHECK(material.m_hole_acoustic.overlap_factor == doctest::Approx(0.5));

    REQUIRE(material.m_electron_intervalley_transitions.size() == 5);
    CHECK(material.m_electron_intervalley_transitions[0].m_name == "g1_TA");
    CHECK(material.m_electron_intervalley_transitions[4].m_name == "f2_LA");

    REQUIRE(material.m_hole_optical_transitions.size() == 4);
    CHECK(material.m_hole_optical_transitions[0].name == "hh_to_hh");
    CHECK(material.m_hole_optical_transitions[3].name == "lh_to_lh");
}

TEST_CASE("parameterized rates reproduce silicon reference values") {
    const auto material = uepm::amc::make_silicon_amc_material_model();

    const double electron_acoustic = uepm::amc::acoustic_scattering_rate(
        material.m_electron_valleys.front(), material.m_electron_acoustic, 0.1, 300.0);
    CHECK(electron_acoustic == doctest::Approx(6.988811279204761e12).epsilon(1.0e-12));

    const double hole_acoustic = uepm::amc::acoustic_scattering_rate(
        material.m_hole_bands.front(), material.m_hole_acoustic, 0.1, 300.0);
    CHECK(hole_acoustic == doctest::Approx(5.152184403984006e12).epsilon(1.0e-12));

    const auto& g3_transition = material.m_electron_intervalley_transitions[2];
    const double g3_absorption = uepm::amc::intervalley_scattering_rate(
        material.m_electron_valleys.front(),
        g3_transition,
        material.m_electron_acoustic.mass_density_kg_per_m3,
        0.1,
        true,
        300.0);
    CHECK(g3_absorption == doctest::Approx(1.500525823504873e11).epsilon(1.0e-12));

    const double screened_impurity = uepm::amc::screened_coulomb_impurity_momentum_relaxation_rate(
        material.m_electron_valleys.front(), material.m_dielectric.epsilon_r, 0.1, 1.0e17, 1.0e17, 300.0);
    CHECK(screened_impurity == doctest::Approx(8.669406979196892e11).epsilon(1.0e-12));
}

TEST_CASE("transition names travel with channels and are recorded") {
    uepm::amc::amc_transport_config config;
    config.m_carrier_type             = uepm::amc::particle_type::electron;
    config.m_enable_impact_ionization = false;

    uepm::amc::amc_transport_kernel transport(config, uepm::amc::make_silicon_amc_material_model(), 1234);
    transport.initialize();

    uepm::amc::particle_amc particle(0, uepm::amc::particle_type::electron);
    transport.initialize_particle_state(particle);
    particle.state().kinetic_energy = 0.2;

    const auto channels = transport.build_scattering_channels(particle);
    const auto transition = std::find_if(channels.begin(), channels.end(), [](const auto& channel) {
        return channel.transition_name == std::string_view{"g3_LO"} &&
               channel.process == uepm::amc::intervalley_process::absorption;
    });

    REQUIRE(transition != channels.end());
    transport.apply_scattering_channel(particle, *transition);

    const auto& recorded_transitions = particle.history().transition_events();
    REQUIRE(recorded_transitions.contains("g3_LO"));
    CHECK(recorded_transitions.at("g3_LO") == 1);

    config.m_carrier_type = uepm::amc::particle_type::hole;
    uepm::amc::amc_transport_kernel hole_transport(
        config, uepm::amc::make_silicon_amc_material_model(), 5678);
    hole_transport.initialize();

    uepm::amc::particle_amc hole(1, uepm::amc::particle_type::hole);
    hole_transport.initialize_particle_state(hole);
    hole.state().valley_index    = 0;
    hole.state().kinetic_energy = 0.2;

    const auto hole_channels = hole_transport.build_scattering_channels(hole);
    const auto hole_transition =
        std::find_if(hole_channels.begin(), hole_channels.end(), [](const auto& channel) {
            return channel.transition_name == std::string_view{"hh_to_lh"} &&
                   channel.process == uepm::amc::intervalley_process::absorption;
        });

    REQUIRE(hole_transition != hole_channels.end());
    hole_transport.apply_scattering_channel(hole, *hole_transition);
    REQUIRE(hole.history().transition_events().contains("hh_to_lh"));
    CHECK(hole.history().transition_events().at("hh_to_lh") == 1);
}
