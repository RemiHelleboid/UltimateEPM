/**
 * @file amc_device_runner.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-09
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <cstddef>
#include <string>

#include "amc_self_consistent_device_simulation_2d.hpp"
#include "amc_self_consistent_device_simulation_3d.hpp"
#include "device_amc_simulation.hpp"
#include "vector.hpp"

namespace uepm::amc {

struct self_consistent_device_amc_run_config {
    std::string mesh_file;
    std::string material_file;
    std::string material_symbol = "Si";
    std::string output_dir;
    std::string simulation_name = "self_consistent_amc";
    std::string command_line;

    mesh::vector3 starting_position{0.0, 0.0, 0.0};
    std::size_t   number_electrons_start = 1;
    std::size_t   number_holes_start     = 0;
    int           seed_random_generator  = 0;

    options_device_amc                   device_options{};
    options_self_consistent_device_amc_2d self_consistent_options_2d{};
    options_self_consistent_device_amc_3d self_consistent_options_3d{};
};

void run_self_consistent_device_amc_simulation(const self_consistent_device_amc_run_config& config);

}  // namespace uepm::amc
