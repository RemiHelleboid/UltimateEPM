#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include "amc_run_manifest.hpp"

TEST_CASE("simulation manifest groups sections and preserves precise values") {
    uepm::amc::simulation_manifest manifest;
    manifest.add("run", "status", "completed");
    manifest.add("simulation", "time_step_s", 1.0e-15);
    manifest.add("run", "elapsed_seconds", 1.25);
    manifest.add("simulation", "enabled", true);
    manifest.add("input", "path", "directory with spaces/input.msh");

    const auto filename = std::filesystem::temp_directory_path() / "ultimate_epm_amc_manifest_test.txt";
    manifest.write(filename);

    std::ifstream stream(filename);
    REQUIRE(stream.is_open());
    const std::string contents{std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};

    CHECK(contents.find("[run]\nstatus = completed\nelapsed_seconds = 1.25") != std::string::npos);
    CHECK(contents.find("[simulation]\ntime_step_s = 1.0000000000000001e-15\nenabled = true") != std::string::npos);
    CHECK(contents.find("[input]\npath = directory with spaces/input.msh") != std::string::npos);
}

TEST_CASE("command line reconstruction quotes arguments containing spaces") {
    const char* arguments[] = {"device_amc.epm", "--outdir", "result directory"};
    CHECK(uepm::amc::command_line_from_arguments(3, arguments) == "device_amc.epm --outdir \"result directory\"");
}
