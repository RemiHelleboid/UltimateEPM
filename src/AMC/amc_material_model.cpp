/**
 * @file amc_material_model.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-26
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "amc_material_model.hpp"

#include <stdexcept>

#include "physical_constants.hpp"

namespace uepm::amc {

constexpr double m0 = uepm::constants::m_e;

namespace {
// Local z -> global z
valley_model::mat3 rotation_local_z_to_global_z() { return valley_model::identity_matrix(); }

// Local z -> global x
valley_model::mat3 rotation_local_z_to_global_x() { return {{{{0.0, 0.0, 1.0}}, {{0.0, 1.0, 0.0}}, {{-1.0, 0.0, 0.0}}}}; }

// Local z -> global y
valley_model::mat3 rotation_local_z_to_global_y() { return {{{{1.0, 0.0, 0.0}}, {{0.0, 0.0, 1.0}}, {{0.0, -1.0, 0.0}}}}; }

} // namespace

dielectric_properties make_silicon_dielectric_properties() {
    dielectric_properties props;
    props.epsilon_r = 11.7;
    return props;
}

carrier_impurity_mobility_parameters make_silicon_impurity_mobility_parameters() {
    carrier_impurity_mobility_parameters parameters;

    parameters.m_electron.m_mu0_cm2_per_V_s    = 1417.0;
    parameters.m_electron.m_mu_min_cm2_per_V_s = 52.2;
    parameters.m_electron.m_n_ref_cm_3         = 9.68e16;
    parameters.m_electron.m_alpha              = 0.68;

    parameters.m_hole.m_mu0_cm2_per_V_s    = 470.5;
    parameters.m_hole.m_mu_min_cm2_per_V_s = 44.9;
    parameters.m_hole.m_n_ref_cm_3         = 2.23e17;
    parameters.m_hole.m_alpha              = 0.70;

    return parameters;
}

carrier_impact_ionization_parameters make_silicon_impact_ionization_parameters() {
    carrier_impact_ionization_parameters parameters;

    parameters.m_electron.m_threshold_eV  = 1.12;
    // parameters.m_electron.m_prefactor_s_1 = 1.0e11;
    parameters.m_electron.m_prefactor_s_1 = 1.1e14;
    parameters.m_electron.m_exponent      = 2.5;

    parameters.m_hole.m_threshold_eV  = 1.49;
    parameters.m_hole.m_prefactor_s_1 = 1.4e12;
    parameters.m_hole.m_exponent      = 3.4;

    return parameters;
}

// Data from Annexe HDR Dollfus
std::vector<valley_model> make_silicon_delta_valleys() {
    std::vector<valley_model> valleys;
    valleys.reserve(6);

    valley_model::parameters p;
    p.transverse_effective_mass   = 0.1905 * m0;
    p.longitudinal_effective_mass = 0.9163 * m0;

    p.dispersion = valley_model::band_type::kane;

    p.non_parabolicity        = 0.5;
    p.energy_offset           = 0.0;
    p.phonon_reference_energy = 0.0;
    p.degeneracy              = 1;

    p.name     = "Delta_x_plus";
    p.rotation = rotation_local_z_to_global_x();
    valleys.emplace_back(p);

    p.name     = "Delta_x_minus";
    p.rotation = rotation_local_z_to_global_x();
    valleys.emplace_back(p);

    p.name     = "Delta_y_plus";
    p.rotation = rotation_local_z_to_global_y();
    valleys.emplace_back(p);

    p.name     = "Delta_y_minus";
    p.rotation = rotation_local_z_to_global_y();
    valleys.emplace_back(p);

    p.name     = "Delta_z_plus";
    p.rotation = rotation_local_z_to_global_z();
    valleys.emplace_back(p);

    p.name     = "Delta_z_minus";
    p.rotation = rotation_local_z_to_global_z();
    valleys.emplace_back(p);

    return valleys;
}

// IDEM
std::vector<valley_model> make_silicon_hole_bands() {
    std::vector<valley_model> bands;
    bands.reserve(2);

    valley_model::parameters p;
    p.dispersion              = valley_model::band_type::parabolic;
    p.non_parabolicity        = 0.0;
    p.energy_offset           = 0.0;
    p.phonon_reference_energy = 0.0;
    p.degeneracy              = 1;
    p.rotation                = valley_model::identity_matrix();

    // Heavy-hole band: isotropic parabolic approximation
    p.name                        = "heavy_hole";
    p.transverse_effective_mass   = 0.87 * m0;
    p.longitudinal_effective_mass = 0.87 * m0;
    bands.emplace_back(p);

    // Light-hole band: isotropic parabolic approximation
    p.name                        = "light_hole";
    p.transverse_effective_mass   = 0.25 * m0;
    p.longitudinal_effective_mass = 0.25 * m0;
    bands.emplace_back(p);

    return bands;
}

std::vector<hole_optical_transition> make_silicon_hole_optical_transitions() {
    constexpr double eV_per_cm_to_eV_per_m = 100.0;

    const double dop_eV_per_m = 8.0e8 * eV_per_cm_to_eV_per_m;
    const double eop_eV       = 63.0e-3;

    std::vector<hole_optical_transition> transitions;
    transitions.reserve(4);

    transitions.push_back({"hh_to_hh", 0, 0, eop_eV, dop_eV_per_m, 0.5});
    transitions.push_back({"hh_to_lh", 0, 1, eop_eV, dop_eV_per_m, 1.0});
    transitions.push_back({"lh_to_hh", 1, 0, eop_eV, dop_eV_per_m, 1.0});
    transitions.push_back({"lh_to_lh", 1, 1, eop_eV, dop_eV_per_m, 0.5});

    return transitions;
}

void amc_material_model::validate() const {
    if (m_name.empty() || m_symbol.empty()) {
        throw std::invalid_argument("AMC material name and symbol must not be empty");
    }
    if (m_dielectric.epsilon_r <= 0.0) {
        throw std::invalid_argument("AMC material relative permittivity must be positive");
    }
    const auto validate_acoustic = [](const acoustic_scattering_parameters& parameters) {
        if (parameters.mass_density_kg_per_m3 <= 0.0) {
            throw std::invalid_argument("AMC material mass density must be positive");
        }
        if (parameters.sound_velocity_m_per_s <= 0.0) {
            throw std::invalid_argument("AMC material sound velocity must be positive");
        }
        if (parameters.deformation_potential_eV < 0.0) {
            throw std::invalid_argument("AMC acoustic deformation potential must be non-negative");
        }
        if (parameters.overlap_factor < 0.0) {
            throw std::invalid_argument("AMC acoustic overlap factor must be non-negative");
        }
    };
    validate_acoustic(m_electron_acoustic);
    validate_acoustic(m_hole_acoustic);
    if (m_electron_valleys.empty() || m_hole_bands.empty()) {
        throw std::invalid_argument("AMC material must define electron valleys and hole bands");
    }
    for (const auto& transition : m_electron_intervalley_transitions) {
        if (transition.m_name.empty()) {
            throw std::invalid_argument("AMC electron intervalley transition name must not be empty");
        }
    }
    for (const auto& transition : m_hole_optical_transitions) {
        if (transition.name.empty()) {
            throw std::invalid_argument("AMC hole optical transition name must not be empty");
        }
        if (transition.initial_band >= m_hole_bands.size() ||
            transition.final_band >= m_hole_bands.size()) {
            throw std::invalid_argument("AMC hole optical transition references an invalid band");
        }
    }
}

amc_material_model make_silicon_amc_material_model() {
    constexpr double silicon_mass_density_kg_per_m3 = 2.329e3;
    constexpr double electron_longitudinal_sound_velocity_m_per_s = 9.002e3;
    constexpr double electron_transverse_sound_velocity_m_per_s   = 5.409e3;

    amc_material_model material;
    material.m_name       = "Silicon";
    material.m_symbol     = "Si";
    material.m_dielectric = make_silicon_dielectric_properties();

    material.m_electron_acoustic = acoustic_scattering_parameters{
        .mass_density_kg_per_m3   = silicon_mass_density_kg_per_m3,
        .sound_velocity_m_per_s   = (electron_longitudinal_sound_velocity_m_per_s +
                                   2.0 * electron_transverse_sound_velocity_m_per_s) /
                                  3.0,
        .deformation_potential_eV = 6.55,
        .overlap_factor           = 1.0,
    };
    material.m_hole_acoustic = acoustic_scattering_parameters{
        .mass_density_kg_per_m3   = silicon_mass_density_kg_per_m3,
        .sound_velocity_m_per_s   = 6.606e3,
        .deformation_potential_eV = 5.5,
        .overlap_factor           = 0.5,
    };

    material.m_electron_valleys                = make_silicon_delta_valleys();
    material.m_electron_intervalley_transitions = make_silicon_intervalley_phonon_branches();
    material.m_hole_bands                      = make_silicon_hole_bands();
    material.m_hole_optical_transitions        = make_silicon_hole_optical_transitions();
    material.m_impurity_mobility               = make_silicon_impurity_mobility_parameters();
    material.m_impact_ionization               = make_silicon_impact_ionization_parameters();
    material.validate();
    return material;
}

}  // namespace uepm::amc
