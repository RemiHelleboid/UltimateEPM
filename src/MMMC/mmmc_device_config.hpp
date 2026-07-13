/**
 * @file mmmc_device_config.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-07-10
 * 
 * @copyright Copyright (c) 2026
 * 
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
