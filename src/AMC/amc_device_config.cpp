/**
 * @file amc_device_config.cpp
 * @brief YAML configuration support for the self-consistent AMC device runner.
 */

#include "amc_device_config.hpp"

#include <fmt/format.h>
#include <yaml-cpp/yaml.h>

#include <fstream>
#include <stdexcept>
#include <string_view>

#include "amc_device_setup.hpp"

namespace uepm::amc {
namespace {

YAML::Node make_default_config() {
    YAML::Node config;

    config["input"]["device_mesh"] = "";
    config["input"]["material_file"] =
        (std::filesystem::path(PROJECT_SRC_DIR) / "examples/materials/materials.yaml").string();
    config["input"]["material"] = "Si";

    config["run"]["name"]             = "self_consistent_amc";
    config["run"]["output_directory"] = "";
    config["run"]["threads"]          = 1;
    config["run"]["seed"]             = 0;

    config["simulation"]["final_time_s"]           = 1.0e-12;
    config["simulation"]["time_step_s"]            = 1.0e-15;
    config["simulation"]["temperature_K"]          = 300.0;
    config["simulation"]["max_particles"]          = 1000000000;
    config["simulation"]["poisson_frequency"]      = 10;
    config["simulation"]["stop_when_no_electrons"] = true;

    config["transport"]["max_energy_eV"]       = 10.0;
    config["transport"]["gamma_safety_factor"] = 1.2;
    config["transport"]["gamma_samples"]       = 1000;
    config["transport"]["impact_ionization"]   = true;
    config["transport"]["particle_creation"]   = true;
    config["transport"]["impurity_scattering"] = false;
    config["transport"]["impurity_model"]      = "mobility";
    config["transport"]["impurity_screening"]  = "debye";

    config["contacts"]["anode_voltage_V"]   = 0.0;
    config["contacts"]["cathode_voltage_V"] = 0.0;

    config["particles"]["initial_electrons"]        = 1;
    config["particles"]["initial_holes"]            = 0;
    config["particles"]["initial_position"]["x_um"] = 0.0;
    config["particles"]["initial_position"]["y_um"] = 0.0;
    config["particles"]["initial_position"]["z_um"] = 0.0;
    config["particles"]["initialize_from_doping"]   = true;
    config["particles"]["initial_weight"]           = 2.0;
    config["particles"]["contact_injection_weight"] = 2.0;

    config["geometry_2d"]["effective_depth_um"]   = 1.0;
    config["geometry_2d"]["particle_z_period_um"] = 1.0;

    config["scheduled_injection"]["enabled"]          = false;
    config["scheduled_injection"]["time_s"]           = 0.0;
    config["scheduled_injection"]["position"]["x_um"] = 0.0;
    config["scheduled_injection"]["position"]["y_um"] = 0.0;
    config["scheduled_injection"]["position"]["z_um"] = 0.0;
    config["scheduled_injection"]["type"]             = "electron";
    config["scheduled_injection"]["weight"]           = 1.0;

    config["output"]["keep_particle_history"] = false;
    config["output"]["export_time_steps"]     = false;
    config["output"]["export_frequency"]      = 100;

    config["quench_circuit"]["enabled"]                   = true;
    config["quench_circuit"]["resistance_ohm"]            = 1.0;
    config["quench_circuit"]["capacitance_F"]             = 1.0;
    config["quench_circuit"]["biased_contact"]            = "cathode";
    config["quench_circuit"]["ramo_current_sign"]         = -1.0;
    config["quench_circuit"]["background_ramo_current_A"] = 0.0;

    config["avalanche_detection"]["voltage_drop_V"]   = 1.0;
    config["quench_detection"]["high_field_V_per_cm"] = 1.0e5;
    config["quench_detection"]["quiet_time_s"]        = 1.0e-11;

    return config;
}

void merge_config(YAML::Node target, const YAML::Node& source, const std::string& path = "") {
    if (!source.IsMap()) {
        throw std::invalid_argument(
            fmt::format("Configuration section '{}' must be a mapping.", path.empty() ? "<root>" : path));
    }

    for (const auto& entry : source) {
        const std::string key       = entry.first.as<std::string>();
        const std::string full_path = path.empty() ? key : path + "." + key;
        if (!target[key]) {
            throw std::invalid_argument(fmt::format("Unknown configuration key '{}'.", full_path));
        }

        const YAML::Node value = entry.second;
        if (target[key].IsMap()) {
            if (!value.IsMap()) {
                throw std::invalid_argument(fmt::format("Configuration section '{}' must be a mapping.", full_path));
            }
            merge_config(target[key], value, full_path);
        } else {
            if (!value.IsScalar()) {
                throw std::invalid_argument(fmt::format("Configuration value '{}' must be a scalar.", full_path));
            }
            target[key] = value;
        }
    }
}

void apply_override(YAML::Node config, const std::string& override_text) {
    const std::size_t equals = override_text.find('=');
    if (equals == std::string::npos || equals == 0 || equals + 1 >= override_text.size()) {
        throw std::invalid_argument(fmt::format("Invalid override '{}'. Expected path.to.value=value.", override_text));
    }

    const std::string path       = override_text.substr(0, equals);
    const std::string value_text = override_text.substr(equals + 1);
    YAML::Node        node       = config;
    std::size_t       begin      = 0;

    while (true) {
        const std::size_t dot = path.find('.', begin);
        const std::string key = path.substr(begin, dot == std::string::npos ? dot : dot - begin);
        if (key.empty() || !node[key]) {
            throw std::invalid_argument(fmt::format("Unknown configuration override '{}'.", path));
        }
        if (dot == std::string::npos) {
            if (!node[key].IsScalar()) {
                throw std::invalid_argument(fmt::format("Configuration override '{}' does not name a value.", path));
            }
            YAML::Node value = YAML::Load(value_text);
            if (!value.IsScalar()) {
                throw std::invalid_argument(fmt::format("Override value for '{}' must be a scalar.", path));
            }
            node[key] = value;
            return;
        }
        if (!node[key].IsMap()) {
            throw std::invalid_argument(fmt::format("Configuration override '{}' traverses a scalar value.", path));
        }
        node.reset(node[key]);
        begin = dot + 1;
    }
}

template <typename T>
T value_at(const YAML::Node& config, std::string_view section, std::string_view key) {
    try {
        return config[std::string(section)][std::string(key)].as<T>();
    } catch (const YAML::Exception& error) {
        throw std::invalid_argument(fmt::format("Invalid value for '{}.{}': {}", section, key, error.what()));
    }
}

template <typename T>
T nested_value_at(const YAML::Node& config,
                  std::string_view  section,
                  std::string_view  subsection,
                  std::string_view  key) {
    try {
        return config[std::string(section)][std::string(subsection)][std::string(key)].as<T>();
    } catch (const YAML::Exception& error) {
        throw std::invalid_argument(
            fmt::format("Invalid value for '{}.{}.{}': {}", section, subsection, key, error.what()));
    }
}

std::string resolve_input_path(const std::filesystem::path& config_file, const std::string& value) {
    if (value.empty()) {
        return value;
    }
    const std::filesystem::path path(value);
    return path.is_absolute() ? path.lexically_normal().string()
                              : (config_file.parent_path() / path).lexically_normal().string();
}

quench_biased_contact parse_biased_contact(const std::string& value) {
    if (value == "anode") {
        return quench_biased_contact::anode;
    }
    if (value == "cathode") {
        return quench_biased_contact::cathode;
    }
    throw std::invalid_argument("quench_circuit.biased_contact must be either anode or cathode.");
}

}  // namespace

self_consistent_device_amc_run_config load_device_amc_config(const std::filesystem::path&    config_file,
                                                             const std::vector<std::string>& overrides) {
    YAML::Node config = make_default_config();

    try {
        merge_config(config, YAML::LoadFile(config_file.string()));
        for (const auto& override_text : overrides) {
            apply_override(config, override_text);
        }
    } catch (const YAML::Exception& error) {
        throw std::invalid_argument(
            fmt::format("Could not read AMC configuration '{}': {}", config_file.string(), error.what()));
    }

    self_consistent_device_amc_run_config result;
    result.mesh_file       = resolve_input_path(config_file, value_at<std::string>(config, "input", "device_mesh"));
    result.material_file   = resolve_input_path(config_file, value_at<std::string>(config, "input", "material_file"));
    result.material_symbol = value_at<std::string>(config, "input", "material");
    result.output_dir      = value_at<std::string>(config, "run", "output_directory");
    result.simulation_name = value_at<std::string>(config, "run", "name");
    result.seed_random_generator = value_at<int>(config, "run", "seed");

    if (result.mesh_file.empty()) {
        throw std::invalid_argument("input.device_mesh is required.");
    }

    auto& device                                  = result.device_options;
    device.m_t_max                                = value_at<double>(config, "simulation", "final_time_s");
    device.m_time_step                            = value_at<double>(config, "simulation", "time_step_s");
    device.m_lattice_temperature                  = value_at<double>(config, "simulation", "temperature_K");
    device.m_max_number_particle                  = value_at<std::size_t>(config, "simulation", "max_particles");
    device.m_stop_simu_when_no_electron_remaining = value_at<bool>(config, "simulation", "stop_when_no_electrons");
    device.m_nb_threads                           = value_at<int>(config, "run", "threads");
    device.m_max_energy_eV                        = value_at<double>(config, "transport", "max_energy_eV");
    device.m_self_scattering_safety_factor        = value_at<double>(config, "transport", "gamma_safety_factor");
    device.m_gamma_max_energy_samples             = value_at<std::size_t>(config, "transport", "gamma_samples");
    device.m_activate_impact_ionization           = value_at<bool>(config, "transport", "impact_ionization");
    device.m_particle_creation_activated          = value_at<bool>(config, "transport", "particle_creation");
    device.m_enable_impurity_scattering           = value_at<bool>(config, "transport", "impurity_scattering");
    device.m_impurity_scattering_model =
        parse_impurity_model(value_at<std::string>(config, "transport", "impurity_model"));
    device.m_impurity_screening_model =
        parse_impurity_screening_model(value_at<std::string>(config, "transport", "impurity_screening"));
    device.m_keep_particles_history      = value_at<bool>(config, "output", "keep_particle_history");
    device.m_export_time_step            = value_at<bool>(config, "output", "export_time_steps");
    device.m_frequency_export_trajectory = value_at<int>(config, "output", "export_frequency");

    device.m_enable_scheduled_particle_injection = value_at<bool>(config, "scheduled_injection", "enabled");
    auto& injection                              = device.m_scheduled_particle_injection;
    injection.m_time_s                           = value_at<double>(config, "scheduled_injection", "time_s");
    injection.m_position_um = mesh::vector3{nested_value_at<double>(config, "scheduled_injection", "position", "x_um"),
                                            nested_value_at<double>(config, "scheduled_injection", "position", "y_um"),
                                            nested_value_at<double>(config, "scheduled_injection", "position", "z_um")};
    injection.m_particle_type = parse_particle_type(value_at<std::string>(config, "scheduled_injection", "type"));
    injection.m_weight        = value_at<double>(config, "scheduled_injection", "weight");

    result.starting_position = mesh::vector3{nested_value_at<double>(config, "particles", "initial_position", "x_um"),
                                             nested_value_at<double>(config, "particles", "initial_position", "y_um"),
                                             nested_value_at<double>(config, "particles", "initial_position", "z_um")};
    result.number_electrons_start = value_at<std::size_t>(config, "particles", "initial_electrons");
    result.number_holes_start     = value_at<std::size_t>(config, "particles", "initial_holes");

    options_self_consistent_device_amc_common common;
    common.m_poisson_frequency                 = value_at<std::size_t>(config, "simulation", "poisson_frequency");
    common.m_anode_voltage                     = value_at<double>(config, "contacts", "anode_voltage_V");
    common.m_cathode_voltage                   = value_at<double>(config, "contacts", "cathode_voltage_V");
    common.m_initialize_particles_from_doping  = value_at<bool>(config, "particles", "initialize_from_doping");
    common.m_initial_particle_weight           = value_at<double>(config, "particles", "initial_weight");
    common.m_contact_injection_particle_weight = value_at<double>(config, "particles", "contact_injection_weight");
    common.m_passive_quench_circuit.m_enabled  = value_at<bool>(config, "quench_circuit", "enabled");
    common.m_passive_quench_circuit.m_resistance_ohm = value_at<double>(config, "quench_circuit", "resistance_ohm");
    common.m_passive_quench_circuit.m_capacitance_F  = value_at<double>(config, "quench_circuit", "capacitance_F");
    common.m_quench_biased_contact =
        parse_biased_contact(value_at<std::string>(config, "quench_circuit", "biased_contact"));
    common.m_ramo_current_to_quench_current_sign = value_at<double>(config, "quench_circuit", "ramo_current_sign");
    common.m_background_ramo_current_A = value_at<double>(config, "quench_circuit", "background_ramo_current_A");
    common.m_avalanche_voltage_drop_threshold_V   = value_at<double>(config, "avalanche_detection", "voltage_drop_V");
    common.m_quench_high_field_threshold_V_per_cm = value_at<double>(config, "quench_detection", "high_field_V_per_cm");
    common.m_quench_quiet_time_s                  = value_at<double>(config, "quench_detection", "quiet_time_s");

    const auto biased_voltage                        = common.m_quench_biased_contact == quench_biased_contact::anode
                                                           ? common.m_anode_voltage
                                                           : common.m_cathode_voltage;
    common.m_passive_quench_circuit.m_bias_voltage_V = biased_voltage;
    common.m_passive_quench_circuit.m_initial_device_voltage_V = biased_voltage;

    result.self_consistent_options_2d.m_common = common;
    result.self_consistent_options_2d.m_effective_depth_um =
        value_at<double>(config, "geometry_2d", "effective_depth_um");
    result.self_consistent_options_2d.m_particle_z_period_um =
        value_at<double>(config, "geometry_2d", "particle_z_period_um");
    result.self_consistent_options_3d.m_common = common;

    device.validate();
    common.validate();
    return result;
}

void write_basic_device_amc_config(const std::filesystem::path& config_file) {
    if (config_file.has_parent_path()) {
        std::filesystem::create_directories(config_file.parent_path());
    }

    YAML::Node config              = make_default_config();
    config["input"]["device_mesh"] = "device.msh";

    std::ofstream stream(config_file);
    if (!stream.is_open()) {
        throw std::runtime_error(fmt::format("Could not write '{}'.", config_file.string()));
    }
    stream << "# Complete AMC device simulation configuration\n"
           << "# Override settings with --set path.to.setting=value or --set path.to.setting value.\n"
           << "# Relative input paths are resolved from this file's directory.\n"
           << config << '\n';
}

}  // namespace uepm::amc
