/**
 * @file admc_device_config.hpp
 * @brief YAML configuration support for self-consistent ADMC device simulation.
 */

#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "admc_device_runner.hpp"

namespace uepm::ADMC {

self_consistent_device_admc_run_config load_device_admc_config(const std::filesystem::path&    config_file,
                                                               const std::vector<std::string>& overrides = {});

void write_basic_device_admc_config(const std::filesystem::path& config_file);

}  // namespace uepm::ADMC
