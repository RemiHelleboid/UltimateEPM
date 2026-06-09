#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "amc_avalanche_detector.hpp"

TEST_CASE("avalanche is detected immediately at the voltage-drop threshold") {
    uepm::amc::voltage_drop_avalanche_detector detector(1.0);

    detector.update(true, 1.0e-12, 0.99);
    CHECK_FALSE(detector.state().m_detected);

    detector.update(true, 2.0e-12, 1.0);
    REQUIRE(detector.state().m_detected);
    REQUIRE(detector.state().m_time_s.has_value());
    REQUIRE(detector.state().m_voltage_drop_V.has_value());
    CHECK(*detector.state().m_time_s == doctest::Approx(2.0e-12));
    CHECK(*detector.state().m_voltage_drop_V == doctest::Approx(1.0));
}

TEST_CASE("avalanche detection uses the absolute voltage drop and latches the first event") {
    uepm::amc::voltage_drop_avalanche_detector detector(0.5);

    detector.update(true, 3.0e-12, -0.6);
    detector.update(true, 4.0e-12, 2.0);

    REQUIRE(detector.state().m_detected);
    CHECK(*detector.state().m_time_s == doctest::Approx(3.0e-12));
    CHECK(*detector.state().m_voltage_drop_V == doctest::Approx(-0.6));
}

TEST_CASE("avalanche detection is inactive when the quench circuit is disabled") {
    uepm::amc::voltage_drop_avalanche_detector detector(0.5);

    detector.update(false, 1.0e-12, 1.0);

    CHECK_FALSE(detector.state().m_detected);
    CHECK_FALSE(detector.state().m_time_s.has_value());
}
