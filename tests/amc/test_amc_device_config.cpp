#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "pbmc_device_config.hpp"

namespace {

std::filesystem::path write_config(const std::string& contents) {
    const auto    path = std::filesystem::temp_directory_path() / "ultimate_epm_pbmc_device_config.yaml";
    std::ofstream stream(path);
    stream << contents;
    return path;
}

}  // namespace

TEST_CASE("PBMC device config loads YAML and applies CLI overrides last") {
    const auto config_file = write_config(R"(
input:
  device_mesh: mesh/device.msh
  material_root: materials
run:
  threads: 3
simulation:
  final_time_s: 2.0e-12
contacts:
  voltages_V:
    anode: 0.0
    cathode: 25.0
scheduled_injection:
  enabled: true
  time_s: 1.0e-12
  position:
    x_um: 1.0
  type: hole
)");

    const auto config = uepm::PBMC::load_device_pbmc_config(
        config_file,
        {"run.threads=8", "contacts.voltages_V.cathode=30.0", "scheduled_injection.weight=4.5"});

    CHECK(config.device_options.m_nb_threads == 8);
    CHECK(config.device_options.m_t_max == doctest::Approx(2.0e-12));
    CHECK(config.self_consistent_options_2d.m_common.m_contact_voltages_V.at("cathode") == doctest::Approx(30.0));
    CHECK(config.device_options.m_scheduled_particle_injection.m_particle_type == uepm::PBMC::particle_type::hole);
    CHECK(config.device_options.m_scheduled_particle_injection.m_weight == doctest::Approx(4.5));
    CHECK(std::filesystem::path(config.mesh_file) == config_file.parent_path() / "mesh/device.msh");
}

TEST_CASE("PBMC device config supports arbitrary named contacts and a Ramo electrode") {
    const auto config_file = write_config(R"(
input:
  device_mesh: nmos.msh
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

    const auto  config = uepm::PBMC::load_device_pbmc_config(config_file, {"contacts.voltages_V.drain=0.25"});
    const auto& common = config.self_consistent_options_2d.m_common;

    CHECK(common.m_contact_voltages_V.size() == 4);
    CHECK(common.m_contact_voltages_V.at("drain") == doctest::Approx(0.25));
    CHECK(common.m_contact_voltages_V.at("gate") == doctest::Approx(1.2));
    CHECK(common.m_ramo_electrode == "drain");
}

TEST_CASE("PBMC device config rejects a Ramo electrode that is not a configured contact") {
    const auto config_file = write_config(R"(
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

    CHECK_THROWS_WITH_AS(uepm::PBMC::load_device_pbmc_config(config_file),
                         "Ramo electrode 'gate' is not present in the contact voltage map.",
                         std::invalid_argument);
}

TEST_CASE("PBMC device config rejects unknown YAML and override keys") {
    const auto config_file = write_config(R"(
input:
  device_mesh: device.msh
simulation:
  typo_time: 1.0
)");

    CHECK_THROWS_WITH_AS(uepm::PBMC::load_device_pbmc_config(config_file),
                         "Unknown configuration key 'simulation.typo_time'.",
                         std::invalid_argument);

    const auto valid_config = write_config(R"(
input:
  device_mesh: device.msh
)");
    CHECK_THROWS_WITH_AS(uepm::PBMC::load_device_pbmc_config(valid_config, {"run.typo=2"}),
                         "Unknown configuration override 'run.typo'.",
                         std::invalid_argument);
}

TEST_CASE("generated PBMC device config contains the complete schema and is loadable") {
    const auto config_file = std::filesystem::temp_directory_path() / "ultimate_epm_basic_pbmc_device_config.yaml";
    uepm::PBMC::write_basic_device_pbmc_config(config_file);

    std::ifstream stream(config_file);
    REQUIRE(stream.is_open());
    const std::string contents{std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};

    CHECK(contents.find("transport:") != std::string::npos);
    CHECK(contents.find("particles:") != std::string::npos);
    CHECK(contents.find("scheduled_injection:") != std::string::npos);
    CHECK(contents.find("geometry_2d:") != std::string::npos);
    CHECK(contents.find("output:") != std::string::npos);
    CHECK(contents.find("quench_circuit:") != std::string::npos);
    CHECK(contents.find("avalanche_detection:") != std::string::npos);
    CHECK(contents.find("quench_detection:") != std::string::npos);

    const auto config = uepm::PBMC::load_device_pbmc_config(config_file);
    CHECK(config.device_options.m_time_step == doctest::Approx(1.0e-15));
    CHECK(config.device_options.m_nb_threads == 1);
}
