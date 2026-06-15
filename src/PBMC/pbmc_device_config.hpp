/**
 * @file pbmc_device_config.hpp
 * @brief YAML configuration support for the self-consistent PBMC device runner.
 */

#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "pbmc_device_runner.hpp"

namespace uepm::PBMC {

self_consistent_device_pbmc_run_config load_device_pbmc_config(const std::filesystem::path&    config_file,
                                                             const std::vector<std::string>& overrides = {});

void write_basic_device_pbmc_config(const std::filesystem::path& config_file);

}  // namespace uepm::PBMC
