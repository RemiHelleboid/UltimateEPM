#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "amc_device_config.hpp"

namespace {

std::filesystem::path write_config(const std::string& contents) {
    const auto path =
        std::filesystem::temp_directory_path() / "ultimate_epm_amc_device_config.yaml";
    std::ofstream stream(path);
    stream << contents;
    return path;
}

}  // namespace

TEST_CASE("AMC device config loads YAML and applies CLI overrides last") {
    const auto config_file = write_config(R"(
input:
  device_mesh: mesh/device.msh
  material_file: materials.yaml
run:
  threads: 3
simulation:
  final_time_s: 2.0e-12
contacts:
  cathode_voltage_V: 25.0
scheduled_injection:
  enabled: true
  time_s: 1.0e-12
  position:
    x_um: 1.0
  type: hole
)");

    const auto config = uepm::amc::load_device_amc_config(
        config_file,
        {"run.threads=8", "contacts.cathode_voltage_V=30.0", "scheduled_injection.weight=4.5"});

    CHECK(config.device_options.m_nb_threads == 8);
    CHECK(config.device_options.m_t_max == doctest::Approx(2.0e-12));
    CHECK(config.self_consistent_options_2d.m_common.m_cathode_voltage == doctest::Approx(30.0));
    CHECK(config.device_options.m_scheduled_particle_injection.m_particle_type ==
          uepm::amc::particle_type::hole);
    CHECK(config.device_options.m_scheduled_particle_injection.m_weight == doctest::Approx(4.5));
    CHECK(std::filesystem::path(config.mesh_file) == config_file.parent_path() / "mesh/device.msh");
}

TEST_CASE("AMC device config rejects unknown YAML and override keys") {
    const auto config_file = write_config(R"(
input:
  device_mesh: device.msh
simulation:
  typo_time: 1.0
)");

    CHECK_THROWS_WITH_AS(uepm::amc::load_device_amc_config(config_file),
                         "Unknown configuration key 'simulation.typo_time'.",
                         std::invalid_argument);

    const auto valid_config = write_config(R"(
input:
  device_mesh: device.msh
)");
    CHECK_THROWS_WITH_AS(uepm::amc::load_device_amc_config(valid_config, {"run.typo=2"}),
                         "Unknown configuration override 'run.typo'.",
                         std::invalid_argument);
}

TEST_CASE("generated AMC device config contains the complete schema and is loadable") {
    const auto config_file =
        std::filesystem::temp_directory_path() / "ultimate_epm_basic_amc_device_config.yaml";
    uepm::amc::write_basic_device_amc_config(config_file);

    std::ifstream stream(config_file);
    REQUIRE(stream.is_open());
    const std::string contents{std::istreambuf_iterator<char>(stream),
                               std::istreambuf_iterator<char>()};

    CHECK(contents.find("transport:") != std::string::npos);
    CHECK(contents.find("particles:") != std::string::npos);
    CHECK(contents.find("scheduled_injection:") != std::string::npos);
    CHECK(contents.find("geometry_2d:") != std::string::npos);
    CHECK(contents.find("output:") != std::string::npos);
    CHECK(contents.find("quench_circuit:") != std::string::npos);
    CHECK(contents.find("avalanche_detection:") != std::string::npos);
    CHECK(contents.find("quench_detection:") != std::string::npos);

    const auto config = uepm::amc::load_device_amc_config(config_file);
    CHECK(config.device_options.m_time_step == doctest::Approx(1.0e-15));
    CHECK(config.device_options.m_nb_threads == 1);
}
