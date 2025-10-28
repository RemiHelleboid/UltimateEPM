#pragma once

#include <cmath>
#include <map>
#include <string>

#include "NonLocalParameters.hpp"
#include "SpinOrbitParameters.hpp"
#include "Pseudopotential.h"
#include "yaml-cpp/yaml.h"

namespace uepm::pseudopotential {

class Material {
 protected:
    /**
     * @brief Name of the material. (Actually it is the symbol, eg.g., Si)
     *
     */
    std::string m_name;

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

    /**
     * @brief Non-local parameters strcut of the material.
     *
     */
    NonLocalParameters m_non_local_parameters;

    /**
     * @brief Spin-orbit coupling parameters struct of the material.
     *
     */
    SpinOrbitParameters m_spin_orbit_parameters;

    /**
     * @brief Flag to indicate if the non-local corrections parameters are set.
     *
     */
    bool m_is_non_local_parameters_populated = false;

    /**
     * @brief Flag to indicate if the spin-orbit coupling parameters are set.
     *
     */
    bool m_is_spin_orbit_parameters_populated = false;

 public:
    Material() : m_lattice_constant(0) {}
    Material(const std::string& Name,
             double             a,
             double             V3S,
             double             V4S,
             double             V8S,
             double             V11S,
             double             V3A  = 0,
             double             V4A  = 0,
             double             V8A  = 0,
             double             V11A = 0);
    Material(const Material&) = default;
    Material& operator=(const Material&) = default;
    Material(Material&&)                 = default;
    
    const std::string& get_name() const { return m_name; }

    /**
     * @brief Populate the non-local parameters of the material from a YAML parameter node.
     * 
     * @param node 
     */
    void populate_non_local_parameters(const YAML::Node& node) {
        m_non_local_parameters.populate_non_local_parameters(node);
    }

    void populate_spin_orbit_parameters(const YAML::Node& node) {
        m_spin_orbit_parameters.populate_from_yaml(node);
    }

    const SpinOrbitParameters& get_spin_orbit_parameters() const { return m_spin_orbit_parameters; }

    /**
     * @brief Compute the non-local correction for (K1, K2) element of the Hamiltonian. 
     * 
     * @param K1 
     * @param K2 
     * @param tau 
     * @return std::complex<double> 
     */
    std::complex<double> compute_pseudopotential_non_local_correction(const Vector3D<double>& K1,
                                                                      const Vector3D<double>& K2,
                                                                      const Vector3D<double>& tau) const;

    /**
     * @brief Get the lattice constant in meter. 
     * 
     * @return double 
     */
    double get_lattice_constant_meter() const {
        constexpr double angstrom_to_m = 1.0e-10;
        return angstrom_to_m * m_lattice_constant;
    }

    double get_fourier_factor() const { 
        constexpr double angstrom_to_m = 1.0e-10;
        return 2.0 * M_PI / (angstrom_to_m * m_lattice_constant);
    }

    /**
     * @brief Get the atomic volume of the material in angstrom^3.
     * 
     * @return double 
     */
    double get_atomic_volume_angstrom() const { return (1.0 / 8.0) * pow(m_lattice_constant, 3.0); }

    /**
     * @brief Get the atomic volume of the material in meter^3.
     * 
     * @return double 
     */
    double get_atomic_volume() const { return (1.0 / 8.0) * pow(get_lattice_constant_meter(), 3.0); }

    /**
     * @brief Check if the non-local parameters are populated.
     * 
     * @return true 
     * @return false 
     */
    bool is_non_local_parameters_populated() const { return m_is_non_local_parameters_populated; }

    /**
     * @brief Check if the spin-orbit parameters are populated.
     * 
     * @return true 
     * @return false 
     */
    bool is_spin_orbit_parameters_populated() const { return m_is_spin_orbit_parameters_populated; }

    /**
     * @brief Set the flag to indicate if the non-local parameters are populated.
     * 
     * @param is_non_local_parameters_populated 
     */
    void set_is_non_local_parameters_populated(bool is_non_local_parameters_populated) {
        m_is_non_local_parameters_populated = is_non_local_parameters_populated;
    }

    /**
     * @brief Set the flag to indicate if the spin-orbit parameters are populated.
     * 
     * @param is_spin_orbit_parameters_populated 
     */
    void set_is_spin_orbit_parameters_populated(bool is_spin_orbit_parameters_populated) {
        m_is_spin_orbit_parameters_populated = is_spin_orbit_parameters_populated;
    }

    /**
     * @brief Return the fermi momentum of the material: (96*pi^2)**(1/3) / lattice_constant
     *
     */
    double get_fermi_momentum() const { return pow(96 * M_PI * M_PI, 1.0 / 3.0) / get_lattice_constant_meter(); }

    const Pseudopotential& get_pseudopotential() const { return m_pseudopotential; }

    const NonLocalParameters& get_non_local_parameters() const { return m_non_local_parameters; }
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

}  // namespace uepm::pseudopotential