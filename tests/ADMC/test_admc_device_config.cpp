#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include "admc_device_config.hpp"

namespace {

std::filesystem::path write_config(const std::string& contents) {
    const auto path = std::filesystem::temp_directory_path() / "ultimate_epm_admc_device_config.yaml";
    std::ofstream(path) << contents;
    return path;
}

}  // namespace

TEST_CASE("ADMC config parses boundary reflection model") {
    const auto config_file = write_config(R"(
input:
  device_mesh: device.msh
transport:
  boundary_reflection: diffuse
contacts:
  voltages_V:
    anode: 0.0
    cathode: 1.0
  collecting:
    anode: true
    cathode: true
  ramo_electrode: anode
)");

    const auto config = uepm::ADMC::load_device_admc_config(config_file);
    CHECK(config.device_options.m_boundary_reflection_model == uepm::mesh::boundary_reflection_model::diffuse);
}

TEST_CASE("ADMC config parses scheduled contact voltage events") {
    const auto config_file = write_config(R"(
input:
  device_mesh: device.msh
contacts:
  voltages_V:
    anode: 0.0
    cathode: 0.0
    gate: 0.0
  collecting:
    anode: true
    cathode: true
  ramo_electrode: anode
contact_voltage_schedule:
  enabled: true
  events:
    - time_s: 5.0e-12
      voltages_V:
        cathode: 0.5
    - time_s: 2.0e-12
      voltages_V:
        gate: 1.0
)");

    const auto config = uepm::ADMC::load_device_admc_config(config_file);
    const auto& schedule = config.self_consistent_options_2d.m_common.m_contact_voltage_schedule;

    REQUIRE(schedule.size() == 2);
    CHECK(schedule[0].m_time_s == doctest::Approx(2.0e-12));
    CHECK(schedule[0].m_contact_voltages_V.at("gate") == doctest::Approx(1.0));
    CHECK(schedule[1].m_time_s == doctest::Approx(5.0e-12));
    CHECK(schedule[1].m_contact_voltages_V.at("cathode") == doctest::Approx(0.5));
}
