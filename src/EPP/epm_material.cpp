#include "epm_material.hpp"

#include <cmath>
#include <filesystem>
#include <optional>

#include "bessel_func.hpp"
#include "physical_constants.hpp"
#include "yaml-cpp/yaml.h"

namespace uepm::pseudopotential {

const double Bohr = 0.52917721092;  // in Angstroms

epm_material::epm_material(const uepm::physics::material_info& material,
                           double                              V3S,
                           double                              V4S,
                           double                              V8S,
                           double                              V11S,
                           double                              V3A,
                           double                              V4A,
                           double                              V8A,
                           double                              V11A)
    : m_material_info(material),
      m_pseudopotential(V3S, V4S, V8S, V11S, V3A, V4A, V8A, V11A) {
    if (m_material_info.lattice_constant_m <= 0.0) {
        throw std::invalid_argument("EPP material '" + material.symbol +
                                    "' requires a positive common lattice constant.");
    }
}

namespace {
void load_material_from_yaml(std::map<std::string, epm_material>& materials,
                             const uepm::physics::material_info&  common_material,
                             const std::filesystem::path&         filename,
                             const YAML::Node&                    config,
                             const std::optional<std::string>&    expected_parameter_set = std::nullopt) {
    const auto material_node = config["material"];
    const auto model_node    = config["model"];
    if (!material_node || !model_node || material_node.as<std::string>() != common_material.symbol ||
        model_node.as<std::string>() != "epm") {
        throw std::runtime_error("Invalid EPM parameter file '" + filename.string() + "'.");
    }
    if (expected_parameter_set) {
        const auto parameter_set_node = config["parameter_set"];
        if (!parameter_set_node || parameter_set_node.as<std::string>() != *expected_parameter_set) {
            throw std::runtime_error("Invalid EPM parameter set in file '" + filename.string() + "'.");
        }
    }
    const auto node_pseudopotential = config["pseudo-potential-parameters"];
    if (!node_pseudopotential) {
        throw std::runtime_error("EPM parameter file '" + filename.string() +
                                 "' does not define pseudo-potential-parameters.");
    }

    const auto rydberg_value = [&](const char* key, double default_value = 0.0) {
        const auto value = node_pseudopotential[key];
        return value ? uepm::constants::Ryd_to_eV * value.as<double>() : default_value;
    };

    materials[common_material.symbol] = epm_material(common_material,
                                                     rydberg_value("V3S"),
                                                     rydberg_value("V4S"),
                                                     rydberg_value("V8S"),
                                                     rydberg_value("V11S"),
                                                     rydberg_value("V3A"),
                                                     rydberg_value("V4A"),
                                                     rydberg_value("V8A"),
                                                     rydberg_value("V11A"));

    if (const auto node = config["non-local-parameters"]) {
        materials[common_material.symbol].populate_non_local_parameters(node);
    }
    if (const auto node = config["spin-orbit-parameters"]) {
        materials[common_material.symbol].populate_spin_orbit_parameters(node);
    }
}
}  // namespace

void Materials::load_material(const uepm::physics::material_repository& repository,
                              const std::string&                        material_symbol,
                              const std::string&                        parameter_set) {
    const auto common_material = repository.load_material(material_symbol);
    const auto filename        = repository.parameter_file(material_symbol, "epm", parameter_set);
    const auto config          = YAML::LoadFile(filename.string());
    load_material_from_yaml(materials, common_material, filename, config, parameter_set);
}

void Materials::load_material_file(const uepm::physics::material_repository& repository,
                                   const std::string&                        material_symbol,
                                   const std::filesystem::path&              parameter_file) {
    const auto common_material = repository.load_material(material_symbol);
    const auto filename        = std::filesystem::absolute(parameter_file).lexically_normal();
    const auto config          = YAML::LoadFile(filename.string());
    load_material_from_yaml(materials, common_material, filename, config);
}

void Materials::load_parameter_set(const uepm::physics::material_repository& repository,
                                   const std::string&                        parameter_set) {
    materials.clear();
    for (const auto& symbol : repository.material_symbols()) {
        if (repository.has_parameter_set(symbol, "epm", parameter_set)) {
            load_material(repository, symbol, parameter_set);
        }
    }
}

/**
 * @brief Compute the so called F_l function, which is used in the non-local pseudopotential correction.
 * (See Chelikowsky, J. R. & Cohen, M. L. Nonlocal pseudopotential calculations for the electronic structure of eleven
 * diamond and zinc-blende semiconductors. Phys. Rev. B 14, 556–582 (1976).)
 *
 * For the values of function, see: Bloomfield, J. K., Face, S. H. P. & Moss, Z. Indefinite Integrals of Spherical
 * Bessel Functions. Preprint at http://arxiv.org/abs/1703.06428 (2017). Equations 49 and 59.
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
    // std::cout << "F_l_function: norm_K1 = " << norm_K1 << ", norm_K2 = " << norm_K2 << std::endl;
    if (fabs(norm_K1 - norm_K2) > EPSILON) {
        const double pre_factor = pow(atomic_radii, 2.0) / (norm_K1 * norm_K1 - norm_K2 * norm_K2);
        const double F =
            norm_K1 * generalized_bessel(l + 1, norm_K1 * atomic_radii) *
                generalized_bessel(l, norm_K2 * atomic_radii) -
            norm_K2 * generalized_bessel(l + 1, norm_K2 * atomic_radii) * generalized_bessel(l, norm_K1 * atomic_radii);
        return pre_factor * F;
    } else if (norm_K1 > EPSILON) {
        const double pre_factor = pow(atomic_radii, 3.0) / (2.0);
        // const double pre_factor = 1.0 / (2.0 * atomi c_radii * atomic_radii);
        const double F =
            pow(generalized_bessel(l, norm_K1 * atomic_radii), 2.0) -
            generalized_bessel(l - 1, norm_K1 * atomic_radii) * generalized_bessel(l + 1, norm_K1 * atomic_radii);
        return pre_factor * F;
    } else {
        return (l == 0) ? pow(atomic_radii, 3.0) / (3.0) : 0.0;
    }
}

double F_2_function_gaussian(const Vector3D<double>& K1, const Vector3D<double>& K2, double atomic_radii) {
    const double norm_K1    = K1.Length();
    const double norm_K2    = K2.Length();
    const double bessel_arg = 0.5 * (atomic_radii * atomic_radii) * norm_K1 * norm_K2;
    return bessel_2nd_order_first_kind(bessel_arg) *
           exp(-0.25 * (norm_K1 * norm_K1 + norm_K2 * norm_K2) * atomic_radii * atomic_radii);
}

/**
 * @brief Compute the non local correction to the EPM Hamiltonian.
 * It follows: Chelikowsky, J. R. & Cohen, M. L. Nonlocal pseudopotential calculations for the electronic structure of
 * eleven diamond and zinc-blende semiconductors. Phys. Rev. B 14, 556–582 (1976). See also: Pötz, W. & Vogl, P. Theory
 * of optical-phonon deformation potentials in tetrahedral semiconductors. Phys. Rev. B 24, 2025–2037 (1981)
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
std::complex<double> epm_material::compute_pseudopotential_non_local_correction(const Vector3D<double>& K1_normalized,
                                                                                const Vector3D<double>& K2_normalized,
                                                                                const Vector3D<double>& tau) const {
    const double diag_factor    = pow(uepm::constants::h_bar, 2) / (2.0 * uepm::constants::m_e * uepm::constants::q_e);
    const double fourier_factor = 2.0 * M_PI / get_lattice_constant_meter();
    const Vector3D<double> G_diff_normalized = (K1_normalized - K2_normalized);
    const Vector3D<double> K1                = K1_normalized * fourier_factor;
    const Vector3D<double> K2                = K2_normalized * fourier_factor;
    const double           norm_K1           = K1.Length();
    const double           norm_K2           = K2.Length();
    const double           cos_angle_K1_K2   = compute_cos_angle(K1, K2);
    const double           V_pre_factor      = 4.0 * M_PI / get_atomic_volume();
    const double           legendre_0        = 1.0;
    const double           legendre_2        = 0.5 * (3 * cos_angle_K1_K2 * cos_angle_K1_K2 - 1);

    // First atomic species: anion
    double V_anion = 0;
    // l = 0
    const double A_0_anion =
        m_non_local_parameters.m_alpha_0_anion + diag_factor * m_non_local_parameters.m_beta_0_anion *
                                                     (norm_K1 * norm_K2 - pow(this->get_fermi_momentum(), 2.0));
    const double F_0_anion =
        (m_non_local_parameters.m_R0_anion == 0.0) ? 0.0 : F_l_function(K1, K2, m_non_local_parameters.m_R0_anion, 0);
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
            V_anion_2 = 5.0 * pow(M_PI, 1.5) * (pow(m_non_local_parameters.m_R2_anion, 3.0) / get_atomic_volume()) *
                        A_2_anion * legendre_2 * F_2_anion;
        }
        V_anion += V_anion_2;
    }

    // Second atomic species: cation
    double V_cation = 0;
    // l = 0
    const double F_0_cation =
        (m_non_local_parameters.m_R0_cation == 0) ? 0.0 : F_l_function(K1, K2, m_non_local_parameters.m_R0_cation, 0);
    const double A_0_cation =
        m_non_local_parameters.m_alpha_0_cation + m_non_local_parameters.m_beta_0_cation * diag_factor *
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
            V_cation_2 = 5.0 * pow(M_PI, 1.5) * (pow(m_non_local_parameters.m_R2_cation, 3.0) / get_atomic_volume()) *
                         A_2_cation * legendre_2 * F_2_cation;
        }
        V_cation += V_cation_2;
    }

    const double V_symmetric     = 1.0 * (V_anion + V_cation) / 2.0;
    const double V_antisymmetric = 1.0 * (V_anion - V_cation) / 2.0;

    constexpr double const_two        = 2.0;
    const double     lattice_constant = this->get_lattice_constant_meter();
    const double     Gtau             = (tau / lattice_constant) * (G_diff_normalized);

    return std::complex<double>(cos(const_two * M_PI * Gtau) * V_symmetric,
                                sin(const_two * M_PI * Gtau) * V_antisymmetric);
}

void Materials::print_materials_list() const {
    for (const auto& material : materials) {
        std::cout << material.first << std::endl;
    }
}

void Materials::print_material_parameters(const std::string& name) const {
    if (materials.find(name) == materials.end()) {
        std::cout << "epm_material " << name << " not found" << std::endl;
        return;
    }
    const epm_material& material = materials.at(name);
    std::cout << "epm_material: " << name << std::endl;
    std::cout << "Lattice constant: " << material.get_lattice_constant_meter() << " Bohr" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
}

void Materials::print_material_parameters() const {
    for (const auto& material : materials) {
        print_material_parameters(material.first);
    }
}

}  // namespace uepm::pseudopotential
