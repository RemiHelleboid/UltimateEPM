/**
 * @file admc_device_config.cpp
 * @brief YAML configuration support for self-consistent ADMC device simulation.
 */

#include "admc_device_config.hpp"

#include <fmt/format.h>
#include <yaml-cpp/yaml.h>

#include <fstream>
#include <stdexcept>
#include <string_view>

namespace uepm::ADMC {
namespace {

YAML::Node make_default_config() {
    YAML::Node config;
    config["input"]["device_mesh"]                           = "";
    config["input"]["material_root"]                         = "";
    config["input"]["material"]                              = "Si";
    config["run"]["name"]                                    = "self_consistent_ADMC";
    config["run"]["output_directory"]                        = "";
    config["run"]["seed"]                                    = 0;
    config["run"]["threads"]                                 = 1;
    config["simulation"]["final_time_s"]                     = 100.0e-12;
    config["simulation"]["time_step_s"]                      = 1.0e-15;
    config["simulation"]["temperature_K"]                    = 300.0;
    config["simulation"]["max_particles"]                    = 100000000;
    config["simulation"]["poisson_frequency"]                = 10;
    config["simulation"]["stop_when_no_electrons"]           = false;
    config["transport"]["impact_ionization"]                 = false;
    config["transport"]["particle_creation"]                 = false;
    config["transport"]["impurity_scattering"]               = false;
    config["transport"]["impurity_model"]                    = "mobility";
    config["transport"]["impurity_screening"]                = "debye";
    config["transport"]["max_energy_eV"]                     = 1.0;
    config["transport"]["gamma_safety_factor"]               = 1.2;
    config["transport"]["gamma_samples"]                     = 1000;
    config["transport"]["boundary_reflection"]               = "reverse";
    config["contacts"]["voltages_V"]["anode"]                = 0.0;
    config["contacts"]["voltages_V"]["cathode"]              = 0.0;
    config["contacts"]["collecting"]["anode"]                = true;
    config["contacts"]["collecting"]["cathode"]              = true;
    config["contacts"]["ramo_electrode"]                     = "anode";
    config["contacts"]["apply_built_in_potential"]           = true;
    config["contacts"]["built_in_voltage_scale"]             = 1.0;
    config["poisson_mixing"]["enabled"]                      = false;
    config["poisson_mixing"]["old_solution_fraction"]        = 0.0;
    config["particles"]["initial_electrons"]                 = 0;
    config["particles"]["initial_holes"]                     = 0;
    config["particles"]["initial_position"]["x_um"]          = 0.0;
    config["particles"]["initial_position"]["y_um"]          = 0.0;
    config["particles"]["initial_position"]["z_um"]          = 0.0;
    config["particles"]["initialize_from_doping"]            = true;
    config["particles"]["initial_weight"]                    = 1.0;
    config["particles"]["initial_state_file"]                = "";
    config["particles"]["contact_injection_weight"]          = 1.0;
    config["geometry_2d"]["effective_depth_um"]              = 1.0;
    config["output"]["export_time_steps"]                    = true;
    config["output"]["export_frequency"]                     = 1000;
    config["output"]["keep_particle_history"]                = false;
    config["output"]["mesh_particle_local_averages"]         = true;
    config["scheduled_injection"]["enabled"]                 = false;
    config["scheduled_injection"]["time_s"]                  = 0.0;
    config["scheduled_injection"]["position"]["x_um"]        = 0.0;
    config["scheduled_injection"]["position"]["y_um"]        = 0.0;
    config["scheduled_injection"]["position"]["z_um"]        = 0.0;
    config["scheduled_injection"]["type"]                    = "electron";
    config["scheduled_injection"]["weight"]                  = 1.0;
    config["current_probe"]["enabled"]                       = false;
    config["current_probe"]["x_min_um"]                      = 0.0;
    config["current_probe"]["x_max_um"]                      = 0.0;
    config["current_probe"]["y_min_um"]                      = 0.0;
    config["current_probe"]["y_max_um"]                      = 0.0;
    config["current_probe"]["z_min_um"]                      = 0.0;
    config["current_probe"]["z_max_um"]                      = 0.0;
    config["quench_circuit"]["enabled"]                      = false;
    config["quench_circuit"]["resistance_ohm"]               = 1.0;
    config["quench_circuit"]["capacitance_F"]                = 1.0;
    config["quench_circuit"]["biased_contact"]               = "cathode";
    config["quench_circuit"]["ramo_current_sign"]            = -1.0;
    config["quench_circuit"]["background_ramo_current_A"]    = 0.0;
    config["quench_circuit"]["auto_background_ramo_current"] = false;
    config["avalanche_detection"]["voltage_drop_V"]          = 1.0;
    config["quench_detection"]["high_field_V_per_cm"]        = 1.0e5;
    config["quench_detection"]["quiet_time_s"]               = 1.0e-11;
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
        if (path == "contacts.voltages_V" || path == "contacts.collecting") {
            if (!entry.second.IsScalar()) {
                throw std::invalid_argument(fmt::format("Configuration value '{}' must be a scalar.", full_path));
            }
            target[key] = entry.second;
            continue;
        }
        if (!target[key]) {
            throw std::invalid_argument(fmt::format("Unknown configuration key '{}'.", full_path));
        }
        if (target[key].IsMap()) {
            if (!entry.second.IsMap()) {
                throw std::invalid_argument(fmt::format("Configuration section '{}' must be a mapping.", full_path));
            }
            merge_config(target[key], entry.second, full_path);
        } else {
            if (!entry.second.IsScalar()) {
                throw std::invalid_argument(fmt::format("Configuration value '{}' must be a scalar.", full_path));
            }
            target[key] = entry.second;
        }
    }
}

void apply_override(YAML::Node config, const std::string& override_text) {
    const std::size_t equals = override_text.find('=');
    if (equals == std::string::npos || equals == 0 || equals + 1 >= override_text.size()) {
        throw std::invalid_argument(fmt::format("Invalid override '{}'. Expected path.to.value=value.", override_text));
    }
    const std::string          path       = override_text.substr(0, equals);
    const std::string          value_text = override_text.substr(equals + 1);
    constexpr std::string_view contact_voltage_prefix = "contacts.voltages_V.";
    if (path.starts_with(contact_voltage_prefix)) {
        const std::string contact_name = path.substr(contact_voltage_prefix.size());
        if (contact_name.empty()) {
            throw std::invalid_argument("Contact voltage override must include a contact name.");
        }
        YAML::Node value = YAML::Load(value_text);
        if (!value.IsScalar()) {
            throw std::invalid_argument(fmt::format("Override value for '{}' must be a scalar.", path));
        }
        config["contacts"]["voltages_V"][contact_name] = value;
        return;
    }
    constexpr std::string_view collecting_contact_prefix = "contacts.collecting.";
    if (path.starts_with(collecting_contact_prefix)) {
        const std::string contact_name = path.substr(collecting_contact_prefix.size());
        if (contact_name.empty()) {
            throw std::invalid_argument("Collecting-contact override must include a contact name.");
        }
        YAML::Node value = YAML::Load(value_text);
        if (!value.IsScalar()) {
            throw std::invalid_argument(fmt::format("Override value for '{}' must be a scalar.", path));
        }
        config["contacts"]["collecting"][contact_name] = value;
        return;
    }

    YAML::Node  node  = config;
    std::size_t begin = 0;
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
    return config[std::string(section)][std::string(key)].as<T>();
}

template <typename T>
T nested_value_at(const YAML::Node& config,
                  std::string_view  section,
                  std::string_view  subsection,
                  std::string_view  key) {
    return config[std::string(section)][std::string(subsection)][std::string(key)].as<T>();
}

std::string resolve_input_path(const std::filesystem::path& config_file, const std::string& value) {
    if (value.empty()) {
        return value;
    }
    const std::filesystem::path path(value);
    return path.is_absolute() ? path.lexically_normal().string()
                              : (config_file.parent_path() / path).lexically_normal().string();
}

}  // namespace

self_consistent_device_admc_run_config load_device_admc_config(const std::filesystem::path&    config_file,
                                                               const std::vector<std::string>& overrides) {
    YAML::Node config = make_default_config();
    YAML::Node explicit_contact_voltages;
    YAML::Node explicit_collecting_contacts;
    try {
        const YAML::Node user_config = YAML::LoadFile(config_file.string());
        if (user_config["contacts"] && user_config["contacts"]["voltages_V"]) {
            explicit_contact_voltages = YAML::Clone(user_config["contacts"]["voltages_V"]);
        }
        if (user_config["contacts"] && user_config["contacts"]["collecting"]) {
            explicit_collecting_contacts = YAML::Clone(user_config["contacts"]["collecting"]);
        }
        merge_config(config, user_config);
        for (const auto& override_text : overrides) {
            if (override_text.starts_with("contacts.voltages_V.")) {
                if (!explicit_contact_voltages) {
                    explicit_contact_voltages = YAML::Clone(config["contacts"]["voltages_V"]);
                }
                const std::size_t equals = override_text.find('=');
                const std::string contact_name =
                    override_text.substr(std::string_view("contacts.voltages_V.").size(),
                                         equals - std::string_view("contacts.voltages_V.").size());
                explicit_contact_voltages[contact_name] = YAML::Load(override_text.substr(equals + 1));
            }
            if (override_text.starts_with("contacts.collecting.")) {
                if (!explicit_collecting_contacts) {
                    explicit_collecting_contacts = YAML::Clone(config["contacts"]["collecting"]);
                }
                const std::size_t equals = override_text.find('=');
                const std::string contact_name =
                    override_text.substr(std::string_view("contacts.collecting.").size(),
                                         equals - std::string_view("contacts.collecting.").size());
                explicit_collecting_contacts[contact_name] = YAML::Load(override_text.substr(equals + 1));
            }
            apply_override(config, override_text);
        }
    } catch (const YAML::Exception& error) {
        throw std::invalid_argument(
            fmt::format("Could not read ADMC configuration '{}': {}", config_file.string(), error.what()));
    }

    self_consistent_device_admc_run_config result;
    result.mesh_file       = resolve_input_path(config_file, value_at<std::string>(config, "input", "device_mesh"));
    result.material_root   = resolve_input_path(config_file, value_at<std::string>(config, "input", "material_root"));
    result.material_symbol = value_at<std::string>(config, "input", "material");
    result.output_dir      = value_at<std::string>(config, "run", "output_directory");
    result.simulation_name = value_at<std::string>(config, "run", "name");
    result.random_seed     = value_at<std::uint64_t>(config, "run", "seed");
    if (result.mesh_file.empty()) {
        throw std::invalid_argument("input.device_mesh is required.");
    }

    const YAML::Node contacts =
        explicit_contact_voltages ? explicit_contact_voltages : config["contacts"]["voltages_V"];
    for (const auto& entry : contacts) {
        result.self_consistent_options_2d.m_common.m_contact_voltages_V.emplace(entry.first.as<std::string>(),
                                                                                entry.second.as<double>());
    }
    const YAML::Node collecting =
        explicit_collecting_contacts ? explicit_collecting_contacts : config["contacts"]["collecting"];
    for (const auto& entry : collecting) {
        if (entry.second.as<bool>()) {
            result.collecting_contacts.push_back(entry.first.as<std::string>());
        }
    }
    if (result.collecting_contacts.empty()) {
        throw std::invalid_argument("At least one contacts.collecting entry must be enabled.");
    }

    auto& device                                 = result.device_options;
    device.m_final_time_s                        = value_at<double>(config, "simulation", "final_time_s");
    device.m_time_step_s                         = value_at<double>(config, "simulation", "time_step_s");
    device.m_lattice_temperature_K               = value_at<double>(config, "simulation", "temperature_K");
    device.m_max_number_particles                = value_at<std::size_t>(config, "simulation", "max_particles");
    device.m_stop_when_no_electrons              = value_at<bool>(config, "simulation", "stop_when_no_electrons");
    device.m_export_time_step                    = value_at<bool>(config, "output", "export_time_steps");
    device.m_frequency_export                    = value_at<int>(config, "output", "export_frequency");
    device.m_export_mesh_particle_local_averages = value_at<bool>(config, "output", "mesh_particle_local_averages");
    device.m_boundary_reflection_model =
        mesh::parse_boundary_reflection_model(value_at<std::string>(config, "transport", "boundary_reflection"));
    device.m_enable_scheduled_particle_injection = value_at<bool>(config, "scheduled_injection", "enabled");
    device.m_scheduled_injection_time_s          = value_at<double>(config, "scheduled_injection", "time_s");
    device.m_scheduled_injection_position_um =
        mesh::vector3{nested_value_at<double>(config, "scheduled_injection", "position", "x_um"),
                      nested_value_at<double>(config, "scheduled_injection", "position", "y_um"),
                      nested_value_at<double>(config, "scheduled_injection", "position", "z_um")};
    const auto injection_type = value_at<std::string>(config, "scheduled_injection", "type");
    if (injection_type == "electron" || injection_type == "e") {
        device.m_scheduled_injection_type = carrier_type::electron;
    } else if (injection_type == "hole" || injection_type == "h") {
        device.m_scheduled_injection_type = carrier_type::hole;
    } else {
        throw std::invalid_argument("scheduled_injection.type must be electron/e or hole/h.");
    }
    device.m_scheduled_injection_weight = value_at<double>(config, "scheduled_injection", "weight");
    device.m_current_probe.m_enabled     = value_at<bool>(config, "current_probe", "enabled");
    device.m_current_probe.m_box_um =
        mesh::bbox{value_at<double>(config, "current_probe", "x_min_um"),
                   value_at<double>(config, "current_probe", "x_max_um"),
                   value_at<double>(config, "current_probe", "y_min_um"),
                   value_at<double>(config, "current_probe", "y_max_um"),
                   value_at<double>(config, "current_probe", "z_min_um"),
                   value_at<double>(config, "current_probe", "z_max_um")};

    result.number_electrons_start = value_at<std::size_t>(config, "particles", "initial_electrons");
    result.number_holes_start     = value_at<std::size_t>(config, "particles", "initial_holes");
    result.starting_position_um =
        mesh::vector3{nested_value_at<double>(config, "particles", "initial_position", "x_um"),
                      nested_value_at<double>(config, "particles", "initial_position", "y_um"),
                      nested_value_at<double>(config, "particles", "initial_position", "z_um")};

    auto& common                               = result.self_consistent_options_2d.m_common;
    common.m_poisson_frequency                 = value_at<std::size_t>(config, "simulation", "poisson_frequency");
    common.m_ramo_electrode                    = value_at<std::string>(config, "contacts", "ramo_electrode");
    common.m_enable_built_in_potential         = value_at<bool>(config, "contacts", "apply_built_in_potential");
    common.m_built_in_contact_voltage_scale    = value_at<double>(config, "contacts", "built_in_voltage_scale");
    common.m_enable_poisson_mixing             = value_at<bool>(config, "poisson_mixing", "enabled");
    common.m_poisson_mixing_old_solution_fraction =
        value_at<double>(config, "poisson_mixing", "old_solution_fraction");
    common.m_initialize_particles_from_doping  = value_at<bool>(config, "particles", "initialize_from_doping");
    common.m_initial_particle_weight           = value_at<double>(config, "particles", "initial_weight");
    common.m_initial_particle_state_file =
        resolve_input_path(config_file, value_at<std::string>(config, "particles", "initial_state_file"));
    if (!common.m_initial_particle_state_file.empty()) {
        throw std::invalid_argument(
            "ADMC accepts particles.initial_state_file for PBMC YAML compatibility, but loading an initial "
            "particle-state file is not implemented for ADMC yet.");
    }
    common.m_contact_injection_particle_weight = value_at<double>(config, "particles", "contact_injection_weight");

    result.self_consistent_options_2d.m_effective_depth_um =
        value_at<double>(config, "geometry_2d", "effective_depth_um");

    result.device_options.validate();
    result.self_consistent_options_2d.validate();
    return result;
}

void write_basic_device_admc_config(const std::filesystem::path& config_file) {
    std::ofstream stream(config_file);
    if (!stream.is_open()) {
        throw std::runtime_error("Could not open ADMC configuration file for writing: " + config_file.string());
    }
    stream << make_default_config();
}

}  // namespace uepm::ADMC
