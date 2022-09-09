#include "Material.h"

#include <cmath>

#include "Constants.hpp"
#include "bessel_func.hpp"
#include "yaml-cpp/yaml.h"

namespace EmpiricalPseudopotential {

const double Bohr = 0.52917721092;  // in Angstroms

Material::Material(const std::string& Name, double a, double V3S, double V8S, double V11S, double V3A, double V4A, double V11A)
    : name(Name),
      m_lattice_constant(a),
      m_pseudopotential(V3S / 1.,
                        V8S / 1.,
                        V11S / 1.,
                        V3A / 1.,
                        V4A / 1.,
                        V11A / 1.)  // the values are in Rydbergs, convert them to Hartree, one Rydberg is half a Hartree
{}

/**
 * @brief Load material parameters from the passed filename.
 * The file is a YAML file containing the paramters for each material (lattice constant, pseudopotential parameters, etc.).
 *
 * @param filename
 */
void Materials::load_material_parameters(const std::string& filename) {
    YAML::Node config         = YAML::LoadFile(filename);
    auto       list_materials = config["materials"];
    for (const auto& material : list_materials) {
        std::string name                 = material["name"].as<std::string>();
        std::string symbol               = material["symbol"].as<std::string>();
        double      a                    = material["lattice-constant"].as<double>();
        auto        node_pseudopotential = material["pseudo-potential-parameters"];
        double      V3S                  = node_pseudopotential["V3S"].as<double>();
        double      V8S                  = node_pseudopotential["V8S"].as<double>();
        double      V11S                 = node_pseudopotential["V11S"].as<double>();
        double      V3A                  = node_pseudopotential["V3A"].as<double>();
        double      V4A                  = node_pseudopotential["V4A"].as<double>();
        double      V11A                 = node_pseudopotential["V11A"].as<double>();
        materials[symbol]                = Material(symbol, a, V3S, V8S, V11S, V3A, V4A, V11A);

        auto node_non_local_parameters = material["non-local-parameters"];
        if (node_non_local_parameters) {
            materials[symbol].m_non_local_parameters.populate_non_local_parameters(node_non_local_parameters);
            materials[symbol].set_is_non_local_parameters_populated(true);
            materials[symbol].m_non_local_parameters.print_non_local_parameters();
        }
    }
}

double compute_F_l_function(const Vector3D<int>& K1, const Vector3D<int>& K2, double atomic_radii, int l) {
    // std::cout << "atomic_radii = " << atomic_radii << std::endl;
    double norm_K1 = K1.Length();
    double norm_K2 = K2.Length();
    if (abs(norm_K1 * norm_K1 - norm_K2 * norm_K2) > 1.0e-9) {
        double pre_factor = (pow(atomic_radii, 2.0) / (norm_K1 * norm_K1 - norm_K2 * norm_K2));
        double F =
            pre_factor * (norm_K1 * generalized_bessel(l + 1, norm_K1 * atomic_radii) * generalized_bessel(l, norm_K2 * atomic_radii) -
                          norm_K2 * generalized_bessel(l + 1, norm_K1 * atomic_radii) * generalized_bessel(l, norm_K2 * atomic_radii));
        return F;
    } else if (fabs(norm_K1) > 1e-9) {
        double F = 0.5 * pow(atomic_radii, 3.0) *
                   (pow(generalized_bessel(l, norm_K1 * atomic_radii), 2.0) -
                    generalized_bessel(l - 1, norm_K1 * atomic_radii) * generalized_bessel(l + 1, norm_K1 * atomic_radii));
        return F;
    } else {
        return pow(atomic_radii, 3.0) / 3.0;
    }
}

/**
 * @brief Compute the non local correction to the EPM Hamiltonian.
 * It follows: Chelikowsky, J. R. & Cohen, M. L. Nonlocal pseudopotential calculations for the electronic structure of eleven diamond
 * and zinc-blende semiconductors. Phys. Rev. B 14, 556–582 (1976).
 * See also: Pötz, W. & Vogl, P. Theory of optical-phonon deformation
 * potentials in tetrahedral semiconductors. Phys. Rev. B 24, 2025–2037 (1981)
 *
 * K1 = (k + G)
 * K2 = (k + G')
 * tau = 1/8 * a * (1, 1, 1)
 *
 * @warning This function aims to be as close as possible to the original implementation of the authors.
 * It might not be the most efficient way, even though the compiler may optimize it for us.
 *
 * @warning Only square well pseudopotential are supported for now on.
 *
 * @param K1
 * @param K2
 * @param tau
 * @return std::complex<double>
 */
std::complex<double> Material::compute_pseudopotential_non_local_correction(const Vector3D<int>&    K1,
                                                                            const Vector3D<int>&    K2,
                                                                            const Vector3D<double>& tau) const {
    // std::cout << "K1: " << K1 << "   K2 " << K2 << std::endl;
    double norm_K1 = K1.Length();
    double norm_K2 = K2.Length();

    // G_diff = (K1 - K2) = (k + G) - (k + G') = G - G'
    const Vector3D<int> G_diff          = K1 - K2;
    double              cos_angle_K1_K2 = compte_cos_angle(K1, K2);
    double              V_pre_factor    = 4.0 * M_PI / get_atomic_volume_angstrom();
    // std::cout << "Cos angle: " << cos_angle_K1_K2 << std::endl;
    // std::cout << "V pre factor: " << V_pre_factor << std::endl;
    // First atomic species: anion
    double A_0_anion = m_non_local_parameters.m_alpha_0_anion +
                       m_non_local_parameters.m_beta_0_anion * (norm_K1 * norm_K2 - pow(this->get_fermi_momentum(), 2.0));
    double A_2_anion = m_non_local_parameters.m_A2_anion;
    double F_0_anion = compute_F_l_function(K1, K2, m_non_local_parameters.m_R0_anion, 0);
    double F_2_anion = compute_F_l_function(K1, K2, m_non_local_parameters.m_R2_anion, 2);

    if (norm_K1 == 0.0 && norm_K2 == 0.0) {
        std::cout << "F_0_anion: " << F_0_anion << std::endl;
        std::cout << "F_2_anion: " << F_2_anion << std::endl;
    }
    double V_anion = 0;
    // l = 0
    V_anion += A_0_anion * (2 * 0 + 1) * std::legendre(0, cos_angle_K1_K2) * F_0_anion;
    // std::cout << "A_0_anion: " << A_0_anion << std::endl;
    // std::cout << "Legendre : " << std::legendre(0, cos_angle_K1_K2) << std::endl;
    // std::cout << "V_anion 1: " << V_anion << std::endl;
    // l = 2
    V_anion += A_2_anion * (2 * 2 + 1) * std::legendre(2, cos_angle_K1_K2) * F_2_anion;
    V_anion *= V_pre_factor;
    // std::cout << "V_anion: " << V_anion << std::endl;

    // Second atomic species: cation
    double A_0_cation = m_non_local_parameters.m_alpha_0_cation +
                        m_non_local_parameters.m_beta_0_cation * (norm_K1 * norm_K2 - pow(this->get_fermi_momentum(), 2.0));
    double A_2_cation = m_non_local_parameters.m_A2_cation;
    double F_0_cation = compute_F_l_function(K1, K2, m_non_local_parameters.m_R0_cation, 0);
    double F_2_cation = compute_F_l_function(K1, K2, m_non_local_parameters.m_R2_cation, 2);
    double V_cation   = 0;
    // l = 0
    V_cation += A_0_cation * (2 * 0 + 1) * std::legendre(0, cos_angle_K1_K2) * F_0_cation;
    // l = 2
    V_cation += A_2_cation * (2 * 2 + 1) * std::legendre(2, cos_angle_K1_K2) * F_2_cation;
    double pre_factor = 4.0 * M_PI / get_atomic_volume_angstrom();
    V_cation *= V_pre_factor;
    // std::cout << "V_cation: " << V_cation << std::endl;

    double V_symmetric     = (V_anion + V_cation) / 2.0;
    double V_antisymmetric = (V_anion - V_cation) / 2.0;

    constexpr double const_two = 2.0;
    const double     Gtau      = const_two * M_PI * tau * G_diff;
    // std::cout << "Vsym: " << V_symmetric << std::endl;
    // std::cout << "Vanti: " << V_antisymmetric << std::endl;

    return std::complex<double>(cos(Gtau) * V_symmetric, sin(Gtau) * V_antisymmetric);
}

void Materials::print_materials_list() const {
    for (const auto& material : materials) {
        std::cout << material.first << std::endl;
    }
}

void Materials::print_material_parameters(const std::string& name) const {
    if (materials.find(name) == materials.end()) {
        std::cout << "Material " << name << " not found" << std::endl;
        return;
    }
    const Material& material = materials.at(name);
    std::cout << "Material: " << name << std::endl;
    std::cout << "Lattice constant: " << material.m_lattice_constant << " Bohr" << std::endl;
    std::cout << "Empirical Pseudopotential Parameters:" << std::endl;
    material.m_pseudopotential.print_parameters();
    std::cout << "-------------------------------------" << std::endl;
}

void Materials::print_material_parameters() const {
    for (const auto& material : materials) {
        print_material_parameters(material.first);
    }
}

}  // namespace EmpiricalPseudopotential
