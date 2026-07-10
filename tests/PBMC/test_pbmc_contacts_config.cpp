#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include "pbmc_device_config.hpp"
#include "pbmc_self_consistent_device_simulation_base.hpp"

namespace {

std::filesystem::path write_contact_config(const std::string& contents) {
    const auto path = std::filesystem::temp_directory_path() / "ultimate_epm_pbmc_contacts.yaml";
    std::ofstream(path) << contents;
    return path;
}

}  // namespace

TEST_CASE("PBMC config accepts named NMOS contacts and one Ramo electrode") {
    const auto path = write_contact_config(R"(
input:
  device_mesh: nmos.msh
geometry_2d:
  effective_depth_um: 0.25
contacts:
  voltages_V:
    source: 0.0
    drain: 0.1
    body: 0.0
    gate: 1.2
  collecting:
    source: true
    drain: true
  ramo_electrode: drain
quench_circuit:
  enabled: false
)");

    const auto  config   = uepm::PBMC::load_device_pbmc_config(path, {"contacts.voltages_V.drain=0.25"});
    const auto& contacts = config.self_consistent_options_2d.m_common.m_contact_voltages_V;

    CHECK(contacts.size() == 4);
    CHECK(contacts.at("drain") == doctest::Approx(0.25));
    CHECK(contacts.at("gate") == doctest::Approx(1.2));
    CHECK(config.self_consistent_options_2d.m_common.m_ramo_electrode == "drain");
    CHECK(config.collecting_contacts.size() == 2);
    CHECK(config.collecting_contacts[0] == "source");
    CHECK(config.collecting_contacts[1] == "drain");
}

TEST_CASE("PBMC config accepts named PN contacts") {
    const auto path = write_contact_config(R"(
input:
  device_mesh: pn.msh
contacts:
  voltages_V:
    anode: 2.0
    cathode: 0.0
  collecting:
    anode: true
    cathode: true
  ramo_electrode: anode
quench_circuit:
  enabled: false
)");

    const auto  config = uepm::PBMC::load_device_pbmc_config(path);
    const auto& common = config.self_consistent_options_2d.m_common;

    CHECK(common.m_contact_voltages_V.size() == 2);
    CHECK(common.m_contact_voltages_V.at("anode") == doctest::Approx(2.0));
    CHECK(common.m_ramo_electrode == "anode");
    CHECK(config.collecting_contacts.size() == 2);
    CHECK(config.collecting_contacts[0] == "anode");
    CHECK(config.collecting_contacts[1] == "cathode");
}

TEST_CASE("PBMC config parses scheduled contact voltage events") {
    const auto path = write_contact_config(R"(
input:
  device_mesh: nmos.msh
contacts:
  voltages_V:
    source: 0.0
    drain: 0.0
    gate: 0.0
  collecting:
    source: true
    drain: true
  ramo_electrode: drain
quench_circuit:
  enabled: false
contact_voltage_schedule:
  enabled: true
  events:
    - time_s: 5.0e-12
      voltages_V:
        drain: 0.5
    - time_s: 2.0e-12
      voltages_V:
        gate: 1.0
)");

    const auto  config   = uepm::PBMC::load_device_pbmc_config(path);
    const auto& schedule = config.self_consistent_options_2d.m_common.m_contact_voltage_schedule;

    REQUIRE(schedule.size() == 2);
    CHECK(schedule[0].m_time_s == doctest::Approx(2.0e-12));
    CHECK(schedule[0].m_contact_voltages_V.at("gate") == doctest::Approx(1.0));
    CHECK(schedule[1].m_time_s == doctest::Approx(5.0e-12));
    CHECK(schedule[1].m_contact_voltages_V.at("drain") == doctest::Approx(0.5));
}

TEST_CASE("PBMC config rejects scheduled voltage for unknown contact") {
    const auto path = write_contact_config(R"(
input:
  device_mesh: nmos.msh
contacts:
  voltages_V:
    source: 0.0
    drain: 0.0
  collecting:
    source: true
    drain: true
  ramo_electrode: drain
quench_circuit:
  enabled: false
contact_voltage_schedule:
  enabled: true
  events:
    - time_s: 2.0e-12
      voltages_V:
        gate: 1.0
)");

    CHECK_THROWS_WITH_AS(uepm::PBMC::load_device_pbmc_config(path),
                         "Scheduled contact-voltage event references unknown contact 'gate'.",
                         std::invalid_argument);
}

TEST_CASE("PBMC config parses boundary reflection model") {
    const auto path = write_contact_config(R"(
input:
  device_mesh: pn.msh
transport:
  boundary_reflection: specular
contacts:
  voltages_V:
    anode: 2.0
    cathode: 0.0
  collecting:
    anode: true
    cathode: true
  ramo_electrode: anode
quench_circuit:
  enabled: false
)");

    const auto config = uepm::PBMC::load_device_pbmc_config(path);
    CHECK(config.device_options.m_boundary_reflection_model == uepm::mesh::boundary_reflection_model::specular);
}

TEST_CASE("PBMC config rejects an unknown Ramo electrode") {
    const auto path = write_contact_config(R"(
input:
  device_mesh: device.msh
contacts:
  voltages_V:
    source: 0.0
    drain: 0.1
  collecting:
    source: true
    drain: true
  ramo_electrode: gate
quench_circuit:
  enabled: false
)");

    CHECK_THROWS_WITH_AS(uepm::PBMC::load_device_pbmc_config(path),
                         "Ramo electrode 'gate' is not present in the contact voltage map.",
                         std::invalid_argument);
}

TEST_CASE("PBMC silicon intrinsic concentration follows temperature") {
    CHECK(uepm::PBMC::silicon_intrinsic_concentration_cm_3(300.0) == doctest::Approx(1.0e10));
    CHECK(uepm::PBMC::silicon_intrinsic_concentration_cm_3(350.0) >
          uepm::PBMC::silicon_intrinsic_concentration_cm_3(300.0));
    CHECK(uepm::PBMC::silicon_intrinsic_concentration_cm_3(250.0) <
          uepm::PBMC::silicon_intrinsic_concentration_cm_3(300.0));
    CHECK_THROWS_AS(uepm::PBMC::silicon_intrinsic_concentration_cm_3(0.0), std::invalid_argument);
}
