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
#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "impact_ionization.hpp"
#include "mobility_model.hpp"
#include "physical_constants.hpp"
#include "physical_model.hpp"

using sp_physical_model = std::shared_ptr<uepm::physic::model::physical_model>;

namespace uepm {

namespace physic {
namespace material {

struct ionization_energy_threshold {
    double m_electron_ionization_energy;
    double m_hole_ionization_energy;
    ionization_energy_threshold() : m_electron_ionization_energy(0.0), m_hole_ionization_energy(0.0) {}
    ionization_energy_threshold(double electron_energy, double hole_energy)
        : m_electron_ionization_energy(electron_energy),
          m_hole_ionization_energy(hole_energy) {}
};

static const ionization_energy_threshold null_ionization_energy_threshold{};

/**
 * @brief Class that represents a material.
 *
 */
struct material {
    std::string                    m_name;
    std::string                    m_formula;
    std::map<std::string, double>  m_parameters;
    model::mobility_model          m_mobility_model;
    model::impact_ionization_model m_impact_ionization_model;

    material() : m_name(""), m_formula(""), m_parameters() {}

    material(const std::string &name, const std::string &formula) : m_name(name), m_formula(formula) {}

    material(const std::string &name, const std::string &formula, const std::map<std::string, double> &parameters)
        : m_name(name),
          m_formula(formula),
          m_parameters(parameters) {}
    double get_absolute_permittivity() const { return constant::vacuum_permittivity * m_parameters.at("dielectric-constant"); }
    double get_scalar_parameter(const std::string &parameter_name) const { return m_parameters.at(parameter_name); }

    /**
     * @brief Get the mobility model.
     *
     * @return const model::mobility_model&
     */
    const model::mobility_model &get_mobility_model() const { return m_mobility_model; }

    /**
     * @brief Set the mobility model.
     *
     * @param mobility_model
     */
    void set_mobility_model(const model::mobility_model &mobility_model) { m_mobility_model = mobility_model; }

    /**
     * @brief Get the impact ionization model.
     *
     * @return const model::impact_ionization_model&
     */
    const model::impact_ionization_model &get_impact_ionization_model() const { return m_impact_ionization_model; }

    /**
     * @brief Set the impact ionization model.
     *
     * @param impact_ionization_model
     */
    void set_impact_ionization_model(const model::impact_ionization_model &impact_ionization_model) {
        m_impact_ionization_model = impact_ionization_model;
    }
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

// static material UnknownMaterial("Unknown", "", -1);
// static material Silicon("Silicon", "Si", 11.9, ionization_energy_threshold{1.8, 2.4});
// static material SiliconDiOxide("Oxide", "SiO2", 3.9);
// static material Oxide("Oxide", "Oxide", 3.9);
// static material Germanium("Germanium", "Ge", 16.0);
// static material Dioxohafnium("Dioxohafnium", "HfO2", 24);
// static material Dioxohafnium2("Hf02", "HfO2", 24);
// static material Silicide("Silicide", "Silicide", 24);
// static material Gas("Gas", "Gas", 1);

// static std::vector<material> list_materials = {UnknownMaterial, Silicon, SiliconDiOxide, Germanium, Gas, Oxide, Dioxohafnium2, Silicide};

// inline material get_material(const std::string &name) {
//     auto it_mat = std::find_if(list_materials.begin(), list_materials.end(),
//                                [&](const auto &myMaterial) { return (myMaterial.m_name == name || myMaterial.m_formula == name); });
//     if (it_mat == list_materials.end()) {
//         return material("Unknown", "", -1.0);
//     }
//     return *it_mat;
// }

}  // namespace material
}  // namespace physic

}  // namespace uepm