/**
 * @file materials.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-09-13
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "materials.hpp"

#include "impact_ionization.hpp"
#include "mobility_model.hpp"
#include "physical_constants.hpp"
#include "yaml-cpp/yaml.h"

namespace uepm {

namespace physic {
namespace material {

static const std::vector<std::string> list_required_parameters = {"dielectric-constant", "e-ionization-threshold", "h-ionization-threshold"};

static const std::vector<std::string> list_optional_parameters = {"lattice_constant", "electron_ionization_energy_threshold",
                                                                  "hole_ionization_energy_threshold"};

/**
 * @brief Load materials parameters from a YAML file.
 *
 * @param filename
 */
void list_materials::load_materials_from_file(const std::string& filename) {
    YAML::Node materials_file = YAML::LoadFile(filename);
    for (const auto& material_node : materials_file) {
        std::string material_name = material_node["name"].as<std::string>();
        std::string material_formula = material_node["symbol"].as<std::string>();
        std::map<std::string, double> material_parameters;
        for (const auto& required_parameter : list_required_parameters) {
            if (!material_node[required_parameter]) {
                throw std::runtime_error("Missing required parameter " + required_parameter + " for material " +
                                         material_node["name"].as<std::string>());
            } else {
                // std::cout << "Parameter " << required_parameter << " found for material " << material_node["name"].as<std::string>()
                //           << std::endl;
                material_parameters[required_parameter] = material_node[required_parameter].as<double>();
            }
        }

        material new_material(material_name, material_formula, material_parameters);
        m_materials.push_back(new_material);
    }
    // Add models for Silicon, that are already implemented in the code.
    for (auto& material : m_materials) {
        if (material.m_name == "Silicon") {
            material.set_mobility_model(model::Arora_Canali_Mobility_Silicon);
            material.set_impact_ionization_model(model::VanOverstratenDeManSilicon300K);
        }
    }
}

/**
 * @brief Check if a material is available in the list of materials.
 *
 * @param material_name
 * @return true
 * @return false
 */
bool list_materials::is_material_available(const std::string& material_name) const {
    auto it = std::find_if(m_materials.begin(), m_materials.end(), [&material_name](const material& arg_material) {
        return (arg_material.m_name == material_name || arg_material.m_formula == material_name);
    });
    return it != m_materials.end();
}

/**
 * @brief Get a material from the list of materials.
 *
 * @param material_name
 * @return const material&
 */
const material& list_materials::get_material(const std::string& material_name) const {
    if (!is_material_available(material_name)) {
        // std::cout << "Available materials are: " << std::endl;
        // print_materials();
        throw std::runtime_error("Material " + material_name + " is not available.");
    }
    auto it_material = std::find_if(m_materials.begin(), m_materials.end(), [&material_name](const material& arg_material) {
        return (arg_material.m_name == material_name || arg_material.m_formula == material_name);
    });
    return *it_material;
}

void list_materials::print_materials() const {
    for (const auto& material : m_materials) {
        std::cout << "Material " << material.m_name << " (" << material.m_formula << ") : " << std::endl;
        for (const auto& parameter : material.m_parameters) {
            std::cout << "\t" << parameter.first << " = " << parameter.second << std::endl;
        }
    }
}

}  // namespace material
}  // namespace physic
} // namespace uepm