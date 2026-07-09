/**
 * @file mmmc_device_config.cpp
 * @brief YAML configuration support for self-consistent MMMC device simulations.
 */

#include "mmmc_device_config.hpp"

#include <fmt/core.h>
#include <yaml-cpp/yaml.h>

#include <filesystem>
#include <fstream>
#include <stdexcept>

#include "pbmc_device_config.hpp"

namespace uepm::MMMC {
namespace {

double yaml_value_or(const YAML::Node& node, const char* key, double fallback) {
    return node && node[key] ? node[key].as<double>() : fallback;
}

void apply_scalar_override(YAML::Node& root, const std::string& assignment) {
    const auto equals = assignment.find('=');
    if (equals == std::string::npos) {
        return;
    }
    const std::string path  = assignment.substr(0, equals);
    const std::string value = assignment.substr(equals + 1);
    if (!path.starts_with("mmmc.pbmc_bbox_um.")) {
        return;
    }
    const std::string key = path.substr(std::string("mmmc.pbmc_bbox_um.").size());
    root["mmmc"]["pbmc_bbox_um"][key] = value;
}

bbox_transport_policy parse_policy(const YAML::Node& config) {
    const auto bbox_node = config["mmmc"] ? config["mmmc"]["pbmc_bbox_um"] : YAML::Node{};
    if (!bbox_node) {
        throw std::invalid_argument("MMMC configuration requires mmmc.pbmc_bbox_um.");
    }
    bbox_transport_policy policy;
    policy.m_pbmc_region_um = mesh::bbox{yaml_value_or(bbox_node, "x_min", 0.0),
                                         yaml_value_or(bbox_node, "x_max", 0.0),
                                         yaml_value_or(bbox_node, "y_min", 0.0),
                                         yaml_value_or(bbox_node, "y_max", 0.0),
                                         yaml_value_or(bbox_node, "z_min", 0.0),
                                         yaml_value_or(bbox_node, "z_max", 0.0)};
    policy.validate();
    return policy;
}

}  // namespace

self_consistent_device_mmmc_run_config load_device_mmmc_config(const std::filesystem::path&    config_file,
                                                               const std::vector<std::string>& overrides) {
    auto pbmc_config = PBMC::load_device_pbmc_config(config_file, overrides);

    YAML::Node raw_config;
    try {
        raw_config = YAML::LoadFile(config_file.string());
    } catch (const std::exception& error) {
        throw std::runtime_error(
            fmt::format("Could not read MMMC configuration '{}': {}", config_file.string(), error.what()));
    }
    for (const auto& override : overrides) {
        apply_scalar_override(raw_config, override);
    }

    self_consistent_device_mmmc_run_config result;
    result.mesh_file                = pbmc_config.mesh_file;
    result.material_root            = pbmc_config.material_root;
    result.material_symbol          = pbmc_config.material_symbol;
    result.output_dir               = pbmc_config.output_dir;
    result.simulation_name          = pbmc_config.simulation_name == "self_consistent_PBMC"
                                          ? "self_consistent_MMMC"
                                          : pbmc_config.simulation_name;
    result.command_line             = pbmc_config.command_line;
    result.collecting_contacts      = pbmc_config.collecting_contacts;
    result.starting_position        = pbmc_config.starting_position;
    result.number_electrons_start   = pbmc_config.number_electrons_start;
    result.number_holes_start       = pbmc_config.number_holes_start;
    result.seed_random_generator    = pbmc_config.seed_random_generator;
    result.device_options.m_pbmc    = pbmc_config.device_options;
    result.device_options.synchronize_from_pbmc();
    result.self_consistent_options_2d.m_common = pbmc_config.self_consistent_options_2d.m_common;
    result.self_consistent_options_2d.m_effective_depth_um =
        pbmc_config.self_consistent_options_2d.m_effective_depth_um;
    result.self_consistent_options_2d.m_policy = parse_policy(raw_config);
    result.device_options.validate();
    result.self_consistent_options_2d.validate();
    return result;
}

void write_basic_device_mmmc_config(const std::filesystem::path& config_file) {
    PBMC::write_basic_device_pbmc_config(config_file);
}

}  // namespace uepm::MMMC
