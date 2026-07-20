#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>

#include "admc_mobility.hpp"
#include "admc_transport.hpp"
#include "bulk_admc_simulation.hpp"
#include "physical_constants.hpp"

namespace {

constexpr double temperature_K = 300.0;

}  // namespace

TEST_CASE("ADMC silicon low-field mobility reproduces the 300 K Arora limits") {
    const uepm::ADMC::silicon_arora_canali_mobility model;
    const auto&                                     parameters = model.parameters();

    const double electron_mobility =
        model.low_field_mobility_m2_per_V_s(uepm::ADMC::carrier_type::electron, temperature_K, 0.0);
    const double hole_mobility =
        model.low_field_mobility_m2_per_V_s(uepm::ADMC::carrier_type::hole, temperature_K, 0.0);

    CHECK(electron_mobility == doctest::Approx(0.1340));
    CHECK(hole_mobility == doctest::Approx(0.04613));
    CHECK(parameters.reference_temperature_K == doctest::Approx(300.0));
    CHECK(parameters.electron.low_field.reference_doping_300_cm_3 == doctest::Approx(1.25e17));
    CHECK(parameters.hole.high_field.saturation_velocity_300_cm_per_s == doctest::Approx(8.37e6));
}

TEST_CASE("ADMC Canali mobility limits drift velocity at high field") {
    const uepm::ADMC::silicon_arora_canali_mobility model;
    constexpr double                                field_V_per_m = 1.0e9;

    const double mobility =
        model.mobility_m2_per_V_s(uepm::ADMC::carrier_type::electron, temperature_K, 0.0, field_V_per_m);
    const double drift_velocity_m_per_s = mobility * field_V_per_m;

    CHECK(drift_velocity_m_per_s == doctest::Approx(1.07e5).epsilon(2.0e-3));
}

TEST_CASE("ADMC Einstein diffusion uses SI units") {
    constexpr double mobility_m2_per_V_s = 0.1;
    const double     diffusion           = uepm::ADMC::einstein_diffusion_m2_per_s(mobility_m2_per_V_s, temperature_K);
    const double     expected = mobility_m2_per_V_s * uepm::constants::k_B * temperature_K / uepm::constants::q_e;

    CHECK(diffusion == doctest::Approx(expected));
}

TEST_CASE("ADMC deterministic drift has opposite electron and hole directions") {
    uepm::ADMC::admc_transport_kernel        kernel;
    const uepm::ADMC::admc_local_environment environment{
        .electric_field_V_per_m    = {1.0e5, 0.0, 0.0},
        .doping_concentration_cm_3 = 0.0,
        .lattice_temperature_K     = temperature_K,
    };
    uepm::ADMC::admc_particle electron(0, uepm::ADMC::carrier_type::electron);
    uepm::ADMC::admc_particle hole(1, uepm::ADMC::carrier_type::hole);

    kernel.step(electron, environment, 1.0e-12, {});
    kernel.step(hole, environment, 1.0e-12, {});

    CHECK(electron.state().position_m.x() < 0.0);
    CHECK(hole.state().position_m.x() > 0.0);
    CHECK(electron.state().position_m.y() == doctest::Approx(0.0));
    CHECK(hole.state().position_m.y() == doctest::Approx(0.0));
}

TEST_CASE("ADMC diffusion displacement uses the supplied normal draw") {
    uepm::ADMC::admc_transport_kernel        kernel;
    const uepm::ADMC::admc_local_environment environment{
        .electric_field_V_per_m    = {},
        .doping_concentration_cm_3 = 0.0,
        .lattice_temperature_K     = temperature_K,
    };
    uepm::ADMC::admc_particle particle(0, uepm::ADMC::carrier_type::electron);
    constexpr double          time_step_s = 2.0e-12;

    kernel.step(particle, environment, time_step_s, {1.0, -2.0, 0.5});

    const double sigma = std::sqrt(2.0 * particle.state().diffusion_m2_per_s * time_step_s);
    CHECK(particle.state().position_m.x() == doctest::Approx(sigma));
    CHECK(particle.state().position_m.y() == doctest::Approx(-2.0 * sigma));
    CHECK(particle.state().position_m.z() == doctest::Approx(0.5 * sigma));
}

TEST_CASE("bulk ADMC runs to the exact requested final time and is reproducible") {
    uepm::ADMC::bulk_admc_simulation_config config;
    config.environment.electric_field_V_per_m = {2.0e5, 0.0, 0.0};
    config.number_electrons                   = 2;
    config.number_holes                       = 1;
    config.time_step_s                        = 3.0e-13;
    config.final_time_s                       = 1.0e-12;
    config.random_seed                        = 42;

    uepm::ADMC::bulk_admc_simulation first(config);
    uepm::ADMC::bulk_admc_simulation second(config);
    first.run();
    second.run();

    REQUIRE(first.particles().size() == 3);
    REQUIRE(second.particles().size() == first.particles().size());
    CHECK(first.current_time_s() == doctest::Approx(config.final_time_s));
    const double first_run_x = first.particles().front().state().position_m.x();
    for (std::size_t i = 0; i < first.particles().size(); ++i) {
        CHECK(first.particles()[i].state().position_m.x() ==
              doctest::Approx(second.particles()[i].state().position_m.x()));
        CHECK(first.particles()[i].state().time_s == doctest::Approx(config.final_time_s));
    }

    first.run();
    CHECK(first.particles().front().state().position_m.x() == doctest::Approx(first_run_x));
}

TEST_CASE("bulk ADMC measured directional diffusion recovers the Einstein coefficient") {
    uepm::ADMC::bulk_admc_simulation_config config;
    config.environment.electric_field_V_per_m = {2.0e5, 0.0, 0.0};
    config.number_electrons                   = 20000;
    config.number_holes                       = 0;
    config.time_step_s                        = 1.0e-11;
    config.final_time_s                       = 1.0e-11;
    config.random_seed                        = 1234;

    uepm::ADMC::bulk_admc_simulation simulation(config);
    simulation.run();

    const double expected = simulation.particles().front().state().diffusion_m2_per_s;
    const auto&  measured = simulation.diffusion_observables().coefficient_m2_per_s;
    CHECK(measured[0] == doctest::Approx(expected).epsilon(0.03));
    CHECK(measured[1] == doctest::Approx(expected).epsilon(0.03));
    CHECK(measured[2] == doctest::Approx(expected).epsilon(0.03));
}

TEST_CASE("bulk ADMC records and exports PBMC-style per-particle histories on request") {
    uepm::ADMC::bulk_admc_simulation_config config;
    config.number_electrons                   = 2;
    config.number_holes                       = 0;
    config.time_step_s                        = 3.0e-13;
    config.final_time_s                       = 1.0e-12;
    config.record_history                     = true;
    config.initial_position_m                 = {2.0e-6, -3.0e-6, 4.0e-6};
    config.environment.electric_field_V_per_m = {2.0e5, 0.0, 0.0};

    uepm::ADMC::bulk_admc_simulation simulation(config);
    simulation.run();

    REQUIRE(simulation.particles().size() == 2);
    for (const auto& particle : simulation.particles()) {
        REQUIRE(particle.history().recorded_number_of_steps() == 5);
        CHECK(particle.history().snapshots().front().time_s == doctest::Approx(0.0));
        CHECK(particle.history().snapshots().front().electric_field_V_per_m.x() ==
              doctest::Approx(config.environment.electric_field_V_per_m.x()));
        CHECK(particle.history().snapshots().front().mobility_m2_per_V_s > 0.0);
        CHECK(particle.history().snapshots().back().time_s == doctest::Approx(config.final_time_s));
    }

    const auto directory = std::filesystem::temp_directory_path() / "ultimate_epm_admc_bulk_history";
    std::filesystem::create_directories(directory);
    const auto prefix = directory / "simulation_results_test";
    simulation.export_particles_history_to_csv(prefix.string());

    for (std::size_t index = 0; index < simulation.particles().size(); ++index) {
        const auto    filename = prefix.string() + "_particle_" + std::to_string(index) + ".csv";
        std::ifstream file(filename);
        REQUIRE(file.is_open());
        std::string header;
        std::getline(file, header);
        CHECK(header == "time,position_x,position_y,position_z,local_k_x,local_k_y,local_k_z,"
                        "velocity_x,velocity_y,velocity_z,lattice_temperature_K,kinetic_energy,gamma,valley_index");
        std::size_t row_count = 0;
        std::string row;
        while (std::getline(file, row)) {
            CHECK(std::count(row.begin(), row.end(), ',') == 13);
            if (row_count == 0) {
                CHECK(row.starts_with("0,2,-3,4,"));
            }
            ++row_count;
        }
        CHECK(row_count == 5);
        std::filesystem::remove(filename);
    }
    std::filesystem::remove(directory);
}

TEST_CASE("bulk ADMC does not retain particle histories by default") {
    uepm::ADMC::bulk_admc_simulation_config config;
    config.number_electrons = 2;
    config.final_time_s     = config.time_step_s;

    uepm::ADMC::bulk_admc_simulation simulation(config);
    simulation.run();

    CHECK(simulation.particles().front().history().snapshots().empty());
}

TEST_CASE("ADMC rejects invalid physical inputs") {
    uepm::ADMC::admc_transport_kernel kernel;
    uepm::ADMC::admc_particle         particle(0, uepm::ADMC::carrier_type::electron);

    CHECK_THROWS_AS(kernel.step(particle, {}, 0.0, {}), std::invalid_argument);
    CHECK_THROWS_AS(uepm::ADMC::einstein_diffusion_m2_per_s(-1.0, temperature_K), std::invalid_argument);
}
