#include "materials.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <utility>

#include "yaml-cpp/yaml.h"

namespace uepm::physics {

namespace {

constexpr double angstrom_to_m = 1.0e-10;

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

double optional_double(const YAML::Node& node, std::string_view key) {
    const auto value = node[std::string(key)];
    return value ? value.as<double>() : 0.0;
}

double read_lattice_constant_m(const YAML::Node& node) {
    if (const auto value = node["lattice_constant_m"]) {
        return value.as<double>();
    }
    if (const auto value = node["lattice-constant"]) {
        return value.as<double>() * angstrom_to_m;
    }
    return 0.0;
}

double read_mass_density_kg_m3(const YAML::Node& node) {
    if (const auto value = node["mass_density_kg_m3"]) {
        return value.as<double>();
    }
    if (const auto value = node["mass-density-kg-m3"]) {
        return value.as<double>();
    }
    return 0.0;
}

double read_static_relative_permittivity(const YAML::Node& node) {
    if (const auto value = node["static_relative_permittivity"]) {
        return value.as<double>();
    }
    return optional_double(node, "dielectric-constant");
}

void validate_non_negative(const material_info& material) {
    if (material.lattice_constant_m < 0.0 || material.mass_density_kg_m3 < 0.0 ||
        material.static_relative_permittivity < 0.0) {
        throw std::runtime_error("epm_material '" + material.symbol + "' contains a negative physical property.");
    }
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
        throw std::invalid_argument("epm_material '" + std::string(name_or_symbol) + "' has no registered material_id.");
    }
    return it->id;
}

void material_database::load_from_file(const std::string& filename) {
    if (!std::filesystem::exists(filename)) {
        throw std::runtime_error("epm_material file '" + filename + "' does not exist.");
    }

    const YAML::Node config = YAML::LoadFile(filename);
    const YAML::Node nodes  = config["materials"] ? config["materials"] : config;
    if (!nodes.IsSequence()) {
        throw std::runtime_error("epm_material file '" + filename + "' must contain a material sequence.");
    }

    m_materials.clear();
    for (const auto& node : nodes) {
        if (!node["name"] || !node["symbol"]) {
            throw std::runtime_error("Every material must define both 'name' and 'symbol'.");
        }

        material_info material;
        material.name                         = node["name"].as<std::string>();
        material.symbol                       = node["symbol"].as<std::string>();
        material.id                           = material_id_from_name_or_symbol(material.symbol);
        material.lattice_constant_m           = read_lattice_constant_m(node);
        material.mass_density_kg_m3           = read_mass_density_kg_m3(node);
        material.static_relative_permittivity = read_static_relative_permittivity(node);
        add(std::move(material));
    }
}

void material_database::add(material_info material) {
    if (material.name.empty() || material.symbol.empty()) {
        throw std::invalid_argument("epm_material name and symbol must not be empty.");
    }
    validate_non_negative(material);
    if (contains(material.id)) {
        throw std::invalid_argument("epm_material '" + material.symbol + "' already exists in material database.");
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
        throw std::runtime_error("epm_material '" + std::string(to_string(id)) +
                                 "' is not available in material database.");
    }
    return *it;
}

const material_info& material_database::require(const std::string& name_or_symbol) const {
    const auto it = std::find_if(m_materials.begin(), m_materials.end(), [&name_or_symbol](const auto& material) {
        return material.name == name_or_symbol || material.symbol == name_or_symbol;
    });
    if (it == m_materials.end()) {
        throw std::runtime_error("epm_material '" + name_or_symbol + "' is not available in material database.");
    }
    return *it;
}

void material_database::print_materials() const {
    for (const auto& material : m_materials) {
        std::cout << "epm_material " << material.name << " (" << material.symbol << ")\n"
                  << "  lattice_constant_m = " << material.lattice_constant_m << '\n'
                  << "  mass_density_kg_m3 = " << material.mass_density_kg_m3 << '\n'
                  << "  static_relative_permittivity = " << material.static_relative_permittivity << '\n';
    }
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
