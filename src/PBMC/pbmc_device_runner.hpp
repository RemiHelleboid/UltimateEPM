/**
 * @file pbmc_device_runner.hpp
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
#include <vector>

#include "pbmc_self_consistent_device_simulation_2d.hpp"
#include "pbmc_self_consistent_device_simulation_3d.hpp"
#include "device_pbmc_simulation.hpp"
#include "vector.hpp"

namespace uepm::PBMC {

struct self_consistent_device_pbmc_run_config {
    std::string mesh_file;
    std::string material_root;
    std::string material_symbol = "Si";
    std::string output_dir;
    std::string simulation_name = "self_consistent_PBMC";
    std::string command_line;
    std::vector<std::string> collecting_contacts;

    mesh::vector3 starting_position{0.0, 0.0, 0.0};
    std::size_t   number_electrons_start = 1;
    std::size_t   number_holes_start     = 0;
    int           seed_random_generator  = 0;

    options_device_PBMC                    device_options{};
    options_self_consistent_device_pbmc_2d self_consistent_options_2d{};
    options_self_consistent_device_pbmc_3d self_consistent_options_3d{};
};

void run_self_consistent_device_pbmc_simulation(const self_consistent_device_pbmc_run_config& config);

}  // namespace uepm::PBMC
