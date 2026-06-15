#include "materials.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <utility>

#include "yaml-cpp/yaml.h"

namespace uepm::physics {

namespace {

struct material_name_entry {
    material_id      id;
    std::string_view name;
    std::string_view symbol;
};

constexpr material_name_entry material_names[] = {
    {material_id::custom, "Custom", "Custom"},
    {material_id::silicon, "Silicon", "Si"},
    {material_id::germanium, "Germanium", "Ge"},
    {material_id::tin, "Tin", "Sn"},
    {material_id::gallium_phosphide, "Gallium Phosphide", "GaP"},
    {material_id::gallium_arsenide, "Gallium Arsenide", "GaAs"},
    {material_id::aluminium_antimonide, "Aluminium Antimonide", "AlSb"},
    {material_id::aluminium_arsenide, "Aluminium Arsenide", "AlAs"},
    {material_id::indium_phosphide, "Indium Phosphide", "InP"},
    {material_id::indium_arsenide, "Indium Arsenide", "InAs"},
    {material_id::indium_antimonide, "Indium Antimonide", "InSb"},
    {material_id::zinc_selenide, "Zinc Selenide", "ZnSe"},
    {material_id::cadmium_telluride, "Cadmium Telluride", "CdTe"},
    {material_id::silicon_dioxide, "Silicon Dioxide", "SiO2"},
    {material_id::gas, "Gas", "Gas"},
};

void validate_non_negative(const material_info& material) {
    if (material.lattice_constant_m < 0.0 || material.mass_density_kg_m3 < 0.0 ||
        material.static_relative_permittivity < 0.0) {
        throw std::runtime_error("Material '" + material.symbol + "' contains a negative physical property.");
    }
}

void validate_path_component(std::string_view value, std::string_view label) {
    if (value.empty() || value == "." || value == ".." || value.find_first_of("/\\") != std::string_view::npos) {
        throw std::invalid_argument("Invalid " + std::string(label) + " '" + std::string(value) + "'.");
    }
}

material_info parse_material_node(const YAML::Node& node) {
    if (!node.IsMap() || !node["schema_version"] || !node["id"] || !node["name"] || !node["symbol"] ||
        !node["lattice_constant_m"] || !node["mass_density_kg_m3"] || !node["static_relative_permittivity"]) {
        throw std::runtime_error("Material files must use the complete schema_version 1 material schema.");
    }
    if (node["schema_version"].as<int>() != 1) {
        throw std::runtime_error("Unsupported material schema version.");
    }

    material_info material;
    material.name                         = node["name"].as<std::string>();
    material.symbol                       = node["symbol"].as<std::string>();
    material.id                           = material_id_from_name_or_symbol(material.symbol);
    material.lattice_constant_m           = node["lattice_constant_m"].as<double>();
    material.mass_density_kg_m3           = node["mass_density_kg_m3"].as<double>();
    material.static_relative_permittivity = node["static_relative_permittivity"].as<double>();
    validate_non_negative(material);
    return material;
}

}  // namespace

std::string_view to_string(material_id id) {
    const auto it = std::find_if(std::begin(material_names), std::end(material_names), [id](const auto& entry) {
        return entry.id == id;
    });
    if (it == std::end(material_names)) {
        throw std::invalid_argument("Unknown material id.");
    }
    return it->symbol;
}

material_id material_id_from_name_or_symbol(std::string_view name_or_symbol) {
    const auto it =
        std::find_if(std::begin(material_names), std::end(material_names), [name_or_symbol](const auto& entry) {
            return entry.name == name_or_symbol || entry.symbol == name_or_symbol;
        });
    if (it == std::end(material_names)) {
        throw std::invalid_argument("Material '" + std::string(name_or_symbol) + "' has no registered material_id.");
    }
    return it->id;
}

void material_database::load_from_file(const std::string& filename) {
    if (!std::filesystem::exists(filename)) {
        throw std::runtime_error("Material file '" + filename + "' does not exist.");
    }

    const YAML::Node config = YAML::LoadFile(filename);
    m_materials.clear();
    add(parse_material_node(config));
}

void material_database::add(material_info material) {
    if (material.name.empty() || material.symbol.empty()) {
        throw std::invalid_argument("Material name and symbol must not be empty.");
    }
    validate_non_negative(material);
    if (contains(material.id)) {
        throw std::invalid_argument("Material '" + material.symbol + "' already exists in material database.");
    }
    m_materials.push_back(std::move(material));
}

bool material_database::contains(material_id id) const {
    return std::any_of(m_materials.begin(), m_materials.end(), [id](const auto& material) {
        return material.id == id;
    });
}

bool material_database::contains(const std::string& name_or_symbol) const {
    return std::any_of(m_materials.begin(), m_materials.end(), [&name_or_symbol](const auto& material) {
        return material.name == name_or_symbol || material.symbol == name_or_symbol;
    });
}

const material_info& material_database::require(material_id id) const {
    const auto it =
        std::find_if(m_materials.begin(), m_materials.end(), [id](const auto& material) { return material.id == id; });
    if (it == m_materials.end()) {
        throw std::runtime_error("Material '" + std::string(to_string(id)) +
                                 "' is not available in material database.");
    }
    return *it;
}

const material_info& material_database::require(const std::string& name_or_symbol) const {
    const auto it = std::find_if(m_materials.begin(), m_materials.end(), [&name_or_symbol](const auto& material) {
        return material.name == name_or_symbol || material.symbol == name_or_symbol;
    });
    if (it == m_materials.end()) {
        throw std::runtime_error("Material '" + name_or_symbol + "' is not available in material database.");
    }
    return *it;
}

void material_database::print_materials() const {
    for (const auto& material : m_materials) {
        std::cout << "Material " << material.name << " (" << material.symbol << ")\n"
                  << "  lattice_constant_m = " << material.lattice_constant_m << '\n'
                  << "  mass_density_kg_m3 = " << material.mass_density_kg_m3 << '\n'
                  << "  static_relative_permittivity = " << material.static_relative_permittivity << '\n';
    }
}

material_repository::material_repository() : material_repository(default_root()) {}

material_repository::material_repository(std::filesystem::path root)
    : m_root(std::filesystem::absolute(std::move(root)).lexically_normal()) {
    if (!std::filesystem::is_directory(m_root)) {
        throw std::runtime_error("Material repository '" + m_root.string() + "' does not exist.");
    }
}

std::filesystem::path material_repository::default_root() {
#ifdef UEPM_MATERIAL_DATA_DIR
    return UEPM_MATERIAL_DATA_DIR;
#else
    return std::filesystem::path("data") / "materials";
#endif
}

std::filesystem::path material_repository::material_file(const std::string& material_symbol) const {
    validate_path_component(material_symbol, "material symbol");
    const auto path = m_root / material_symbol / "material.yaml";
    if (!std::filesystem::is_regular_file(path)) {
        throw std::runtime_error("Material '" + material_symbol + "' is not available in repository '" +
                                 m_root.string() + "'.");
    }
    return path;
}

std::filesystem::path material_repository::parameter_file(const std::string& material_symbol,
                                                          const std::string& module,
                                                          const std::string& parameter_set) const {
    validate_path_component(material_symbol, "material symbol");
    validate_path_component(module, "material module");
    validate_path_component(parameter_set, "parameter set");
    const auto path = m_root / material_symbol / module / (parameter_set + ".yaml");
    if (!std::filesystem::is_regular_file(path)) {
        throw std::runtime_error("Parameter set '" + module + "/" + parameter_set +
                                 "' is not available for material '" + material_symbol + "'.");
    }
    return path;
}

bool material_repository::has_parameter_set(const std::string& material_symbol,
                                            const std::string& module,
                                            const std::string& parameter_set) const {
    try {
        (void)parameter_file(material_symbol, module, parameter_set);
        return true;
    } catch (const std::runtime_error&) {
        return false;
    }
}

material_info material_repository::load_material(const std::string& material_symbol) const {
    material_database database(material_file(material_symbol).string());
    return database.require(material_symbol);
}

material_database material_repository::load_all_materials() const {
    material_database database;
    for (const auto& symbol : material_symbols()) {
        database.add(load_material(symbol));
    }
    return database;
}

std::vector<std::string> material_repository::material_symbols() const {
    std::vector<std::string> symbols;
    for (const auto& entry : std::filesystem::directory_iterator(m_root)) {
        if (entry.is_directory() && std::filesystem::is_regular_file(entry.path() / "material.yaml")) {
            symbols.push_back(entry.path().filename().string());
        }
    }
    std::sort(symbols.begin(), symbols.end());
    return symbols;
}

std::vector<std::string> material_repository::parameter_sets(const std::string& material_symbol,
                                                             const std::string& module) const {
    validate_path_component(material_symbol, "material symbol");
    validate_path_component(module, "material module");
    const auto directory = m_root / material_symbol / module;

    std::vector<std::string> profiles;
    if (!std::filesystem::is_directory(directory)) {
        return profiles;
    }
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (entry.is_regular_file() && entry.path().extension() == ".yaml") {
            profiles.push_back(entry.path().stem().string());
        }
    }
    std::sort(profiles.begin(), profiles.end());
    return profiles;
}

material_info silicon_material_info() {
    return {
        .id                           = material_id::silicon,
        .name                         = "Silicon",
        .symbol                       = "Si",
        .lattice_constant_m           = 5.431e-10,
        .mass_density_kg_m3           = 2329.0,
        .static_relative_permittivity = 11.7,
    };
}

}  // namespace uepm::physics
