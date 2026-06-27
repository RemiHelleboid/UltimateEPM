/**
 * @file admc_device_runner.hpp
 * @brief Runner support for self-consistent ADMC device simulations.
 */

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "self_consistent_device_admc_simulation_2d.hpp"

namespace uepm::ADMC {

struct self_consistent_device_admc_run_config {
    std::string mesh_file;
    std::string material_root;
    std::string material_symbol = "Si";
    std::string output_dir;
    std::string simulation_name = "self_consistent_ADMC";
    std::string command_line;
    std::vector<std::string> collecting_contacts;

    mesh::vector3 starting_position_um{0.0, 0.0, 0.0};
    std::size_t   number_electrons_start = 0;
    std::size_t   number_holes_start     = 0;
    std::uint64_t random_seed            = 5489u;

    options_device_ADMC                    device_options{};
    options_self_consistent_device_ADMC_2d self_consistent_options_2d{};
};

void run_self_consistent_device_admc_simulation(const self_consistent_device_admc_run_config& config);

std::string admc_command_line_from_arguments(int argc, const char* const* argv);

}  // namespace uepm::ADMC
