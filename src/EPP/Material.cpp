#include "Material.h"

#include <cmath>

#include "Constants.hpp"
#include "bessel_func.hpp"
#include "yaml-cpp/yaml.h"

namespace EmpiricalPseudopotential {

const double Bohr = 0.52917721092;  // in Angstroms

Material::Material(const std::string& Name,
                   double             a,
                   double             V3S,
                   double             V4S,
                   double             V8S,
                   double             V11S,
                   double             V3A,
                   double             V4A,
                   double             V8A,
                   double             V11A)
    : m_name(Name),
      m_lattice_constant(a),
      m_pseudopotential(V3S, V4S, V8S, V11S, V3A, V4A, V8A, V11A) {
}

/**
 * @brief Load material parameters from the passed filename.
 * The file is a YAML file containing the parameters for each material (lattice constant, pseudopotential parameters, etc.).
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
        double      V3S                  = Constants::Ryd_to_eV * node_pseudopotential["V3S"].as<double>();
        double      V8S                  = Constants::Ryd_to_eV * node_pseudopotential["V8S"].as<double>();
        double      V11S                 = Constants::Ryd_to_eV * node_pseudopotential["V11S"].as<double>();
        double      V3A                  = Constants::Ryd_to_eV * node_pseudopotential["V3A"].as<double>();
        double      V4A                  = Constants::Ryd_to_eV * node_pseudopotential["V4A"].as<double>();
        double      V11A                 = Constants::Ryd_to_eV * node_pseudopotential["V11A"].as<double>();
        double      V4S                  = 0.0;
        double      V8A                  = 0.0;
        if (node_pseudopotential["V4S"]) {
            V4S = Constants::Ryd_to_eV * node_pseudopotential["V4S"].as<double>();
        }
        if (node_pseudopotential["V8A"]) {
            V8A = Constants::Ryd_to_eV * node_pseudopotential["V8A"].as<double>();
        }
        materials[symbol] = Material(symbol, a, V3S, V4S, V8S, V11S, V3A, V4A, V8A, V11A);

        auto node_non_local_parameters = material["non-local-parameters"];
        if (node_non_local_parameters) {
            materials[symbol].populate_non_local_parameters(node_non_local_parameters);
            materials[symbol].set_is_non_local_parameters_populated(true);
        }

        auto node_spin_orbit_parameters = material["spin-orbit-parameters"];
        if (node_spin_orbit_parameters) {
            materials[symbol].populate_spin_orbit_parameters(node_spin_orbit_parameters);
            materials[symbol].set_is_spin_orbit_parameters_populated(true);
            std::cout << "Spin-orbit parameters for " << symbol << " are set." << std::endl;
            // materials[symbol].get_spin_orbit_parameters().print_parameters();
            
        }
    }
}

/**
 * @brief Compute the so called F_l function, which is used in the non-local pseudopotential correction.
 * (See Chelikowsky, J. R. & Cohen, M. L. Nonlocal pseudopotential calculations for the electronic structure of eleven diamond
 * and zinc-blende semiconductors. Phys. Rev. B 14, 556–582 (1976).)
 *
 * For the values of function, see: Bloomfield, J. K., Face, S. H. P. & Moss, Z. Indefinite Integrals of Spherical Bessel
 * Functions. Preprint at http://arxiv.org/abs/1703.06428 (2017). Equations 49 and 59.
 *
 * @param K1
 * @param K2
 * @param atomic_radii
 * @param l
 * @return double
 */
double F_l_function(const Vector3D<double>& K1, const Vector3D<double>& K2, double atomic_radii, int l) {
    // This epsilon is used to avoid division by zero in the case of K1 == K2.
    // The value is quite big, but lower values lead to numerical instabilities (noisy bands).
    // Reason: K1 and K2 are of the order of  2PI / a_0 ~ 1e10 !
    constexpr double EPSILON = 1.0e-4;
    const double     norm_K1 = K1.Length();
    const double     norm_K2 = K2.Length();
    if (fabs(norm_K1 - norm_K2) > EPSILON) {
        const double pre_factor = pow(atomic_radii, 2.0) / (norm_K1 * norm_K1 - norm_K2 * norm_K2);
        const double F = norm_K1 * generalized_bessel(l + 1, norm_K1 * atomic_radii) * generalized_bessel(l, norm_K2 * atomic_radii) -
                         norm_K2 * generalized_bessel(l + 1, norm_K2 * atomic_radii) * generalized_bessel(l, norm_K1 * atomic_radii);
        return pre_factor * F;
    } else if (norm_K1 > EPSILON) {
        const double pre_factor = pow(atomic_radii, 3.0) / (2.0);
        const double F          = pow(generalized_bessel(l, norm_K1 * atomic_radii), 2.0) -
                         generalized_bessel(l - 1, norm_K1 * atomic_radii) * generalized_bessel(l + 1, norm_K1 * atomic_radii);
        return pre_factor * F;
    } else {
        return (l==0) ? pow(atomic_radii, 3.0) / (3.0) : 0.0;
    }
}

double F_2_function_gaussian(const Vector3D<double>& K1, const Vector3D<double>& K2, double atomic_radii) {
    const double norm_K1    = K1.Length();
    const double norm_K2    = K2.Length();
    const double bessel_arg = 0.5 * (atomic_radii * atomic_radii) * norm_K1 * norm_K2;
    return bessel_2nd_order_first_kind(bessel_arg) * exp(-0.25 * (norm_K1 * norm_K1 + norm_K2 * norm_K2) * atomic_radii * atomic_radii);
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
std::complex<double> Material::compute_pseudopotential_non_local_correction(const Vector3D<double>& K1_normalized,
                                                                            const Vector3D<double>& K2_normalized,
                                                                            const Vector3D<double>& tau) const {
    const double           diag_factor       = pow(Constants::h_bar, 2) / (2.0 * Constants::m_e * Constants::q_e);
    const double           fourier_factor    = 2.0 * M_PI / get_lattice_constant_meter();
    const Vector3D<double> G_diff_normalized = (K1_normalized - K2_normalized);
    const Vector3D<double> K1                = K1_normalized * fourier_factor;
    const Vector3D<double> K2                = K2_normalized * fourier_factor;
    const double           norm_K1           = K1.Length();
    const double           norm_K2           = K2.Length();
    const double           cos_angle_K1_K2   = compte_cos_angle(K1, K2);
    const double           V_pre_factor      = 4.0 * M_PI / get_atomic_volume();
    const double           legendre_0        = 1.0;
    const double           legendre_2        = 0.5 * (3 * cos_angle_K1_K2 * cos_angle_K1_K2 - 1);

    // First atomic species: anion
    double V_anion = 0;
    // l = 0
    const double A_0_anion = m_non_local_parameters.m_alpha_0_anion + diag_factor * m_non_local_parameters.m_beta_0_anion *
                                                                          (norm_K1 * norm_K2 - pow(this->get_fermi_momentum(), 2.0));
    const double F_0_anion = (m_non_local_parameters.m_R0_anion == 0.0) ? 0.0 : F_l_function(K1, K2, m_non_local_parameters.m_R0_anion, 0);
    V_anion += V_pre_factor * A_0_anion * (2 * 0 + 1) * 1.0 * F_0_anion;
    // l = 2
    double V_anion_2 = 0.0;
    if (m_non_local_parameters.m_A2_anion != 0) {
        const double A_2_anion = m_non_local_parameters.m_A2_anion;
        double       F_2_anion = 0.0;
        if (m_non_local_parameters.m_well_type == non_local_well_type::square) {
            F_2_anion = F_l_function(K1, K2, m_non_local_parameters.m_R2_anion, 2);
            V_anion_2 = V_pre_factor * A_2_anion * (2 * 2 + 1) * legendre_2 * F_2_anion;
        } else {
            F_2_anion = F_2_function_gaussian(K1, K2, m_non_local_parameters.m_R2_anion);
            V_anion_2 = 5.0 * pow(M_PI, 1.5) * (pow(m_non_local_parameters.m_R2_anion, 3.0) / get_atomic_volume()) * A_2_anion *
                        legendre_2 * F_2_anion;
        }
        V_anion += V_anion_2;
    }

    // Second atomic species: cation
    double V_cation = 0;
    // l = 0
    const double F_0_cation = (m_non_local_parameters.m_R0_cation == 0) ? 0.0 : F_l_function(K1, K2, m_non_local_parameters.m_R0_cation, 0);
    const double A_0_cation = m_non_local_parameters.m_alpha_0_cation + m_non_local_parameters.m_beta_0_cation * diag_factor *
                                                                            (norm_K1 * norm_K2 - pow(this->get_fermi_momentum(), 2.0));
    V_cation += V_pre_factor * A_0_cation * (2 * 0 + 1) * legendre_0 * F_0_cation;
    // l = 2
    double V_cation_2 = 0.0;
    if (m_non_local_parameters.m_A2_cation != 0.0) {
        const double A_2_cation = m_non_local_parameters.m_A2_cation;
        double       F_2_cation = 0.0;
        if (m_non_local_parameters.m_well_type == non_local_well_type::square) {
            F_2_cation = F_l_function(K1, K2, m_non_local_parameters.m_R2_cation, 2);
            V_cation_2 = V_pre_factor * A_2_cation * (2 * 2 + 1) * legendre_2 * F_2_cation;
        } else {
            F_2_cation = F_2_function_gaussian(K1, K2, m_non_local_parameters.m_R2_cation);
            V_cation_2 = 5.0 * pow(M_PI, 1.5) * (pow(m_non_local_parameters.m_R2_cation, 3.0) / get_atomic_volume()) * A_2_cation *
                         legendre_2 * F_2_cation;
        }
        V_cation += V_cation_2;
    }

    const double V_symmetric     = 1.0 * (V_anion + V_cation) / 2.0;
    const double V_antisymmetric = 1.0 * (V_anion - V_cation) / 2.0;

    constexpr double const_two        = 2.0;
    const double     lattice_constant = this->get_lattice_constant_meter();
    const double     Gtau             = (tau / lattice_constant) * (G_diff_normalized);

    return std::complex<double>(cos(const_two * M_PI * Gtau) * V_symmetric, sin(const_two * M_PI * Gtau) * V_antisymmetric);
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
    std::cout << "Lattice constant: " << material.get_lattice_constant_meter() << " Bohr" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
}

void Materials::print_material_parameters() const {
    for (const auto& material : materials) {
        print_material_parameters(material.first);
    }
}

}  // namespace EmpiricalPseudopotential
