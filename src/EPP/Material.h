#pragma once

#include <cmath>
#include <map>
#include <string>

#include "NonLocalParameters.hpp"
#include "Pseudopotential.h"
#include "yaml-cpp/yaml.h"

namespace EmpiricalPseudopotential {

class Material {
 public:
    /**
     * @brief Name of the material. (Actually it is the symbole, eg.g., Si)
     *
     */
    std::string name;

    /**
     * @brief Lattice constant of the material, in Angrstom.
     *
     */
    double m_lattice_constant;

    /**
     * @brief Pseudopotential class of the material.
     *
     */
    Pseudopotential m_pseudopotential;

    NonLocalParameters m_non_local_parameters;

    bool m_is_non_local_parameters_populated = false;

    Material() : m_lattice_constant(0) {}
    Material(const std::string& Name, double a, double V3S, double V8S, double V11S, double V3A = 0, double V4A = 0, double V11A = 0);
    void populate_non_local_parameters(const YAML::Node& node);

    std::complex<double> compute_pseudopotential_non_local_correction(const Vector3D<double>& K1,
                                                                      const Vector3D<double>& K2,
                                                                      const Vector3D<double>& tau) const;

    double get_lattice_constant_meter() const {
        constexpr double angstrom_to_m = 1.0e-10;
        return angstrom_to_m * m_lattice_constant;
    }

    double get_atomic_volume_angstrom() const { return (1.0 / 8.0) * pow(m_lattice_constant, 3.0); }

    double get_atomic_volume() const { return (1.0 / 8.0) * pow(get_lattice_constant_meter(), 3.0); }

    bool get_is_non_local_parameters_populated() const { return m_is_non_local_parameters_populated; }

    void set_is_non_local_parameters_populated(bool is_non_local_parameters_populated) {
        m_is_non_local_parameters_populated = is_non_local_parameters_populated;
    }

    /**
     * @brief Return the fermi momentum of the material: (96*pi^2)**(1/3) / lattice_constant
     *
     */
    double get_fermi_momentum() const { return pow(96 * M_PI * M_PI, 1.0 / 3.0) / get_lattice_constant_meter(); }
};

class Materials {
 public:
    Materials() = default;
    void load_material_parameters(const std::string& filename);

    std::map<std::string, Material> materials;

    void print_materials_list() const;
    void print_material_parameters(const std::string& material_name) const;
    void print_material_parameters() const;
};

}  // namespace EmpiricalPseudopotential