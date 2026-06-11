#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include "amc_device_history.hpp"

TEST_CASE("empty AMC device history exports a header-only CSV") {
    uepm::amc::history_device_amc history;
    const auto filename =
        std::filesystem::temp_directory_path() / "ultimate_epm_empty_amc_device_history.csv";

    CHECK_NOTHROW(history.export_to_csv(filename.string()));

    std::ifstream stream(filename);
    REQUIRE(stream.is_open());
    const std::string contents{std::istreambuf_iterator<char>(stream),
                               std::istreambuf_iterator<char>()};

    CHECK(contents ==
          "time,nb_electrons,nb_holes,nb_impact_ionization,ramo_current_electron,ramo_current_hole,ramo_current,"
          "max_electric_field,anode_voltage_V,cathode_voltage_V,quench_bias_voltage_V,quench_device_current_A,"
          "quench_resistor_current_A,quench_voltage_drop_V\n");
}
