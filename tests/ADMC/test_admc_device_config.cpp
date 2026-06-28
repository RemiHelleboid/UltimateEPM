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
