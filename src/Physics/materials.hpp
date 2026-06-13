/**
 * @file materials.h
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-05
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once


#include <map>
#include <memory>
#include <string>
#include <vector>

#include "physical_constants.hpp"


namespace uepm {

namespace physic {
namespace material {



/**
 * @brief Class that represents a material.
 *
 */
struct material {
    std::string                    m_name;
    std::string                    m_formula;
    std::map<std::string, double>  m_parameters;

    material() : m_name(""), m_formula(""), m_parameters() {}

    material(const std::string &name, const std::string &formula) : m_name(name), m_formula(formula) {}

    material(const std::string &name, const std::string &formula, const std::map<std::string, double> &parameters)
        : m_name(name),
          m_formula(formula),
          m_parameters(parameters) {}
    double get_absolute_permittivity() const { return uepm::constants::eps_0 * m_parameters.at("dielectric-constant"); }
    double get_scalar_parameter(const std::string &parameter_name) const { return m_parameters.at(parameter_name); }
};

/**
 * @brief Class to store a list of materials.
 *
 */
class list_materials {
 private:
    std::vector<material> m_materials;

 public:
    list_materials() : m_materials() {}
    explicit list_materials(const std::string &filename) : m_materials() { load_materials_from_file(filename); }
    void            add_material(const material &arg_material) { m_materials.push_back(arg_material); }
    void            load_materials_from_file(const std::string &filename);
    bool            is_material_available(const std::string &material_name) const;
    const material &get_material(const std::string &material_name) const;
    void            print_materials() const;
};

}  // namespace material
}  // namespace physic

}  // namespace uepm