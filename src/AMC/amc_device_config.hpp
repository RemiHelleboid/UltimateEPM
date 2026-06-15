/**
 * @file amc_device_config.hpp
 * @brief YAML configuration support for the self-consistent AMC device runner.
 */

#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "amc_device_runner.hpp"

namespace uepm::amc {

self_consistent_device_amc_run_config load_device_amc_config(const std::filesystem::path&    config_file,
                                                             const std::vector<std::string>& overrides = {});

void write_basic_device_amc_config(const std::filesystem::path& config_file);

}  // namespace uepm::amc
