#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "pbmc_quench_detector.hpp"

TEST_CASE("successful quench requires an uninterrupted quiet window after avalanche") {
    uepm::PBMC::successful_quench_detector detector(3.0e-12);

    detector.update(false, 1.0e-12, false, false);
    detector.update(true, 2.0e-12, false, false);
    detector.update(true, 4.0e-12, false, false);
    CHECK_FALSE(detector.state().m_detected);

    detector.update(true, 5.0e-12, false, false);
    REQUIRE(detector.state().m_detected);
    REQUIRE(detector.state().m_time_s.has_value());
    CHECK(*detector.state().m_time_s == doctest::Approx(5.0e-12));
}

TEST_CASE("high-field particles reset the successful-quench quiet window") {
    uepm::PBMC::successful_quench_detector detector(2.0e-12);

    detector.update(true, 1.0e-12, false, false);
    detector.update(true, 2.0e-12, true, false);
    detector.update(true, 3.0e-12, false, false);
    detector.update(true, 4.0e-12, false, false);
    CHECK_FALSE(detector.state().m_detected);

    detector.update(true, 5.0e-12, false, false);
    CHECK(detector.state().m_detected);
}

TEST_CASE("impact ionization resets the successful-quench quiet window") {
    uepm::PBMC::successful_quench_detector detector(2.0e-12);

    detector.update(true, 1.0e-12, false, false);
    detector.update(true, 2.0e-12, false, true);
    detector.update(true, 3.0e-12, false, false);
    detector.update(true, 5.0e-12, false, false);

    REQUIRE(detector.state().m_detected);
    CHECK(*detector.state().m_quiet_period_start_s == doctest::Approx(3.0e-12));
    CHECK(*detector.state().m_time_s == doctest::Approx(5.0e-12));
}
