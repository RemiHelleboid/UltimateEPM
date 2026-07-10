#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <string_view>

#include "pbmc_material_model.hpp"
#include "pbmc_scattering_model.hpp"
#include "pbmc_transport_kernel.hpp"

TEST_CASE("silicon material model preserves scattering parameters") {
    auto common_material                         = uepm::physics::silicon_material_info();
    common_material.mass_density_kg_m3           = 2400.0;
    common_material.static_relative_permittivity = 12.25;
    const auto material                          = uepm::PBMC::make_silicon_pbmc_material_model(common_material);

    CHECK(material.m_id == uepm::physics::material_id::silicon);
    CHECK(material.m_dielectric.epsilon_r == doctest::Approx(12.25));
    CHECK(material.m_electron_acoustic.mass_density_kg_per_m3 == doctest::Approx(2400.0));
    CHECK(material.m_hole_acoustic.mass_density_kg_per_m3 == doctest::Approx(2400.0));
    CHECK(material.m_electron_acoustic.sound_velocity_m_per_s == doctest::Approx(6606.666666666667));
    CHECK(material.m_electron_acoustic.deformation_potential_eV == doctest::Approx(6.55));
    CHECK(material.m_hole_acoustic.sound_velocity_m_per_s == doctest::Approx(6606.0));
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
    const auto material = uepm::PBMC::make_silicon_pbmc_material_model();

    const double electron_acoustic = uepm::PBMC::acoustic_scattering_rate(material.m_electron_valleys.front(),
                                                                          material.m_electron_acoustic,
                                                                          0.1,
                                                                          300.0);
    CHECK(electron_acoustic == doctest::Approx(3.694230372508554e12).epsilon(1.0e-12));

    const double hole_acoustic =
        uepm::PBMC::acoustic_scattering_rate(material.m_hole_bands.front(), material.m_hole_acoustic, 0.1, 300.0);
    CHECK(hole_acoustic == doctest::Approx(5.142829554521578e12).epsilon(1.0e-12));

    const auto&  g3_transition = material.m_electron_intervalley_transitions[2];
    const double g3_absorption =
        uepm::PBMC::intervalley_scattering_rate(material.m_electron_valleys.front(),
                                                g3_transition,
                                                material.m_electron_acoustic.mass_density_kg_per_m3,
                                                0.1,
                                                true,
                                                300.0);
    CHECK(g3_absorption == doctest::Approx(1.500525823504873e11).epsilon(1.0e-12));

    const double screened_impurity =
        uepm::PBMC::screened_coulomb_impurity_momentum_relaxation_rate(material.m_electron_valleys.front(),
                                                                       material.m_dielectric.epsilon_r,
                                                                       0.1,
                                                                       1.0e17,
                                                                       1.0e17,
                                                                       300.0);
    CHECK(screened_impurity == doctest::Approx(8.669406979196892e11).epsilon(1.0e-12));
}

TEST_CASE("finite-temperature impurity screening matches Dawson references") {
    CHECK(uepm::PBMC::impurity_screening_function(0.0) == doctest::Approx(1.0));
    CHECK(uepm::PBMC::impurity_screening_function(0.1) == doctest::Approx(0.9933599239785286).epsilon(2.0e-7));
    CHECK(uepm::PBMC::impurity_screening_function(0.5) == doctest::Approx(0.8488727670040446).epsilon(2.0e-7));
    CHECK(uepm::PBMC::impurity_screening_function(1.0) == doctest::Approx(0.5380795069127684).epsilon(2.0e-7));
    CHECK(uepm::PBMC::impurity_screening_function(4.0) == doctest::Approx(0.03233700030900128).epsilon(2.0e-7));
    CHECK(uepm::PBMC::impurity_screening_function(10.0) == doctest::Approx(0.005025384718759853).epsilon(2.0e-6));
}

TEST_CASE("full impurity screening reproduces independent numerical references") {
    const auto  material = uepm::PBMC::make_silicon_pbmc_material_model();
    const auto& valley   = material.m_electron_valleys.front();

    struct reference_case {
        double energy_eV;
        double density_cm_3;
        double expected_rate_s_1;
    };
    const reference_case cases[] = {
        {0.01, 1.0e15, 3.813939568092906e11},
        {0.10, 1.0e17, 8.726862850796857e11},
        {0.50, 1.0e19, 4.981710987949476e12},
        {1.00, 1.0e17, 4.379353408058030e10},
    };

    for (const auto& test_case : cases) {
        const double rate = uepm::PBMC::screened_coulomb_impurity_momentum_relaxation_rate(
            valley,
            material.m_dielectric.epsilon_r,
            test_case.energy_eV,
            test_case.density_cm_3,
            test_case.density_cm_3,
            300.0,
            uepm::PBMC::impurity_screening_model::finite_temperature_full);
        CHECK(rate == doctest::Approx(test_case.expected_rate_s_1).epsilon(2.0e-6));
    }
}

TEST_CASE("full impurity screening is finite across the validation grid") {
    const auto   material         = uepm::PBMC::make_silicon_pbmc_material_model();
    const auto&  valley           = material.m_electron_valleys.front();
    const double energies_eV[]    = {0.01, 0.05, 0.10, 0.50, 1.00};
    const double densities_cm_3[] = {1.0e15, 1.0e17, 1.0e19};

    for (double energy_eV : energies_eV) {
        for (double density_cm_3 : densities_cm_3) {
            const double analytic_rate = uepm::PBMC::screened_coulomb_impurity_momentum_relaxation_rate(
                valley,
                material.m_dielectric.epsilon_r,
                energy_eV,
                density_cm_3,
                density_cm_3,
                300.0,
                uepm::PBMC::impurity_screening_model::debye_analytic);
            const double full_rate = uepm::PBMC::screened_coulomb_impurity_momentum_relaxation_rate(
                valley,
                material.m_dielectric.epsilon_r,
                energy_eV,
                density_cm_3,
                density_cm_3,
                300.0,
                uepm::PBMC::impurity_screening_model::finite_temperature_full);

            CHECK(std::isfinite(full_rate));
            CHECK(full_rate >= analytic_rate);
        }
    }
}

TEST_CASE("transition names travel with channels and are recorded") {
    uepm::PBMC::pbmc_transport_config config;
    config.m_carrier_type             = uepm::PBMC::particle_type::electron;
    config.m_enable_impact_ionization = false;

    uepm::PBMC::pbmc_transport_kernel transport(config, uepm::PBMC::make_silicon_pbmc_material_model(), 1234);
    transport.initialize();

    uepm::PBMC::pbmc_particle particle(0, uepm::PBMC::particle_type::electron);
    transport.initialize_particle_state(particle);
    particle.state().kinetic_energy = 0.2;

    const auto channels   = transport.build_scattering_channels(particle);
    const auto transition = std::find_if(channels.begin(), channels.end(), [](const auto& channel) {
        return channel.transition_name == std::string_view{"g3_LO"} &&
               channel.process == uepm::PBMC::intervalley_process::absorption;
    });

    REQUIRE(transition != channels.end());
    transport.apply_scattering_channel(particle, *transition);

    const auto& recorded_transitions = particle.history().transition_events();
    REQUIRE(recorded_transitions.contains("g3_LO"));
    CHECK(recorded_transitions.at("g3_LO") == 1);

    config.m_carrier_type = uepm::PBMC::particle_type::hole;
    uepm::PBMC::pbmc_transport_kernel hole_transport(config, uepm::PBMC::make_silicon_pbmc_material_model(), 5678);
    hole_transport.initialize();

    uepm::PBMC::pbmc_particle hole(1, uepm::PBMC::particle_type::hole);
    hole_transport.initialize_particle_state(hole);
    hole.state().valley_index   = 0;
    hole.state().kinetic_energy = 0.2;

    const auto hole_channels   = hole_transport.build_scattering_channels(hole);
    const auto hole_transition = std::find_if(hole_channels.begin(), hole_channels.end(), [](const auto& channel) {
        return channel.transition_name == std::string_view{"hh_to_lh"} &&
               channel.process == uepm::PBMC::intervalley_process::absorption;
    });

    REQUIRE(hole_transition != hole_channels.end());
    hole_transport.apply_scattering_channel(hole, *hole_transition);
    REQUIRE(hole.history().transition_events().contains("hh_to_lh"));
    CHECK(hole.history().transition_events().at("hh_to_lh") == 1);
}
