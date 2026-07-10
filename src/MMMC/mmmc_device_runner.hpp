/**
 * @file mmmc_device_runner.hpp
 * @brief Runner support for self-consistent MMMC device simulations.
 */

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "self_consistent_device_mmmc_simulation_2d.hpp"

namespace uepm::MMMC {

struct self_consistent_device_mmmc_run_config {
    std::string              mesh_file;
    std::string              material_root;
    std::string              material_symbol = "Si";
    std::string              output_dir;
    std::string              simulation_name = "self_consistent_MMMC";
    std::string              command_line;
    std::vector<std::string> collecting_contacts;

    mesh::vector3 starting_position{0.0, 0.0, 0.0};
    std::size_t   number_electrons_start = 1;
    std::size_t   number_holes_start     = 0;
    int           seed_random_generator  = 0;

    options_device_MMMC                    device_options{};
    options_self_consistent_device_MMMC_2d self_consistent_options_2d{};
};

void run_self_consistent_device_mmmc_simulation(const self_consistent_device_mmmc_run_config& config);

std::string mmmc_command_line_from_arguments(int argc, const char* const* argv);

}  // namespace uepm::MMMC
