/**
 * @file mmmc_device_config.hpp
 * @brief YAML configuration support for self-consistent MMMC device simulations.
 */

#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "mmmc_device_runner.hpp"

namespace uepm::MMMC {

self_consistent_device_mmmc_run_config load_device_mmmc_config(const std::filesystem::path&    config_file,
                                                               const std::vector<std::string>& overrides = {});

void write_basic_device_mmmc_config(const std::filesystem::path& config_file);

}  // namespace uepm::MMMC
