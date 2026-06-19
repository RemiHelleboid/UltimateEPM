#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <cmath>
#include <limits>
#include <random>

#include "keldysh_impactio.hpp"
#include "particle.hpp"
#include "single_part_fbmc.hpp"

TEST_CASE("Keldysh impact ionization is zero below threshold and follows the excess-energy law") {
    const uepm::fbmc::KeldyshImpactIonization model(2.0e11, 2.0, 1.5);

    CHECK(model.compute_rate(1.49) == doctest::Approx(0.0));
    CHECK(model.compute_rate(1.5) == doctest::Approx(0.0));
    CHECK(model.compute_rate(2.0) == doctest::Approx(5.0e10));
}

TEST_CASE("Keldysh impact ionization rejects invalid inputs") {
    const uepm::fbmc::KeldyshImpactIonization model;

    CHECK_THROWS_AS(model.compute_rate(std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
    CHECK_THROWS_AS(uepm::fbmc::KeldyshImpactIonization(-1.0, 2.0, 1.0).compute_rate(2.0), std::invalid_argument);
    CHECK_THROWS_AS(uepm::fbmc::KeldyshImpactIonization(1.0, -2.0, 1.0).compute_rate(2.0), std::invalid_argument);
    CHECK_THROWS_AS(uepm::fbmc::KeldyshImpactIonization(1.0, 2.0, -1.0).compute_rate(2.0), std::invalid_argument);
}

TEST_CASE("FBMC particle free-flight draw is reproducible and stores gamma") {
    uepm::fbmc::particle particle(7, uepm::fbmc::particle_type::electron, nullptr);
    particle.set_random_generator(std::mt19937(1234));

    std::mt19937                           expected_rng(1234);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    const double                           expected_u  = distribution(expected_rng);
    const double                           gamma_s_1   = 4.0e12;
    const double                           expected_dt = -std::log(expected_u) / gamma_s_1;

    particle.draw_free_flight_time(gamma_s_1);

    CHECK(particle.state().m_gamma == doctest::Approx(gamma_s_1));
    CHECK(particle.state().m_free_flight_time == doctest::Approx(expected_dt));
}

TEST_CASE("FBMC particle history-derived statistics are safe for empty histories") {
    const uepm::fbmc::particle particle(0, uepm::fbmc::particle_type::electron, nullptr);

    CHECK(particle.compute_mean_energy() == doctest::Approx(0.0));
    CHECK(particle.compute_mean_energy(1.0e-12) == doctest::Approx(0.0));
    CHECK(particle.extract_impact_ionization_coeff() == doctest::Approx(0.0));
    CHECK(particle.extract_global_average_velocity() == doctest::Approx(0.0));
    CHECK(particle.extract_global_average_velocity(1.0e-12) == doctest::Approx(0.0));
}

TEST_CASE("FBMC particle rejects invalid free-flight gamma") {
    uepm::fbmc::particle particle(0, uepm::fbmc::particle_type::electron, nullptr);

    CHECK_THROWS_AS(particle.draw_free_flight_time(0.0), std::runtime_error);
    CHECK_THROWS_AS(particle.draw_free_flight_time(-1.0), std::runtime_error);
    CHECK_THROWS_AS(particle.draw_free_flight_time(std::numeric_limits<double>::infinity()), std::runtime_error);
}

TEST_CASE("FBMC electron and hole acceleration have opposite signs") {
    uepm::fbmc::particle electron(0, uepm::fbmc::particle_type::electron, nullptr);
    uepm::fbmc::particle hole(1, uepm::fbmc::particle_type::hole, nullptr);
    electron.state().m_free_flight_time = 1.0e-15;
    hole.state().m_free_flight_time     = 1.0e-15;

    const uepm::fbmc::vector3 field(1.0e5, 0.0, 0.0);
    electron.update_k_vector(field);
    hole.update_k_vector(field);

    CHECK(electron.get_signed_charge() == doctest::Approx(-1.0));
    CHECK(hole.get_signed_charge() == doctest::Approx(1.0));
    CHECK(electron.state().m_k_vector.x() < 0.0);
    CHECK(hole.state().m_k_vector.x() > 0.0);
    CHECK(hole.state().m_k_vector.x() == doctest::Approx(-electron.state().m_k_vector.x()));
}

TEST_CASE("FBMC particle statistics can discard an initial warmup window") {
    uepm::fbmc::particle particle(0, uepm::fbmc::particle_type::electron, nullptr);

    particle.state().m_time     = 0.0;
    particle.state().m_position = uepm::fbmc::vector3(0.0, 0.0, 0.0);
    particle.state().m_energy   = 0.0;
    particle.update_history();

    particle.state().m_time     = 1.0;
    particle.state().m_position = uepm::fbmc::vector3(10.0, 0.0, 0.0);
    particle.state().m_energy   = 10.0;
    particle.update_history();

    particle.state().m_time     = 3.0;
    particle.state().m_position = uepm::fbmc::vector3(30.0, 0.0, 0.0);
    particle.state().m_energy   = 30.0;
    particle.update_history();

    CHECK(particle.compute_mean_energy() == doctest::Approx(70.0 / 3.0));
    CHECK(particle.compute_mean_energy(1.0) == doctest::Approx(30.0));
    CHECK(particle.compute_mean_energy(0.5) == doctest::Approx(26.0));

    CHECK(particle.extract_global_average_velocity() == doctest::Approx(10.0));
    CHECK(particle.extract_global_average_velocity(1.0) == doctest::Approx(10.0));
    CHECK(particle.extract_global_average_velocity(0.5) == doctest::Approx(10.0));
}

TEST_CASE("FBMC bulk impact-ionization statistics use carrier time and cm^-1 output") {
    uepm::fbmc::impact_ionization_coefficient_statistics stats;
    stats.m_events                         = 20;
    stats.m_carrier_time_s                 = 2.0e-9;
    stats.m_drift_velocity_time_integral_m = 4.0e-4;

    CHECK(stats.event_rate_per_carrier_s_1() == doctest::Approx(1.0e10));
    CHECK(stats.average_drift_velocity_m_per_s() == doctest::Approx(2.0e5));
    CHECK(stats.ionization_coefficient_cm_1() == doctest::Approx(500.0));
}

TEST_CASE("FBMC bulk impact-ionization statistics are safe without samples") {
    const uepm::fbmc::impact_ionization_coefficient_statistics stats;

    CHECK(stats.event_rate_per_carrier_s_1() == doctest::Approx(0.0));
    CHECK(stats.average_drift_velocity_m_per_s() == doctest::Approx(0.0));
    CHECK(stats.ionization_coefficient_cm_1() == doctest::Approx(0.0));
}
