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

carrier_impact_ionization_parameters make_silicon_impact_ionization_parameters() {
    carrier_impact_ionization_parameters parameters;

    parameters.m_electron.m_threshold_eV  = 1.2;
    parameters.m_electron.m_prefactor_s_1 = 1.0e11;
    parameters.m_electron.m_exponent      = 6.0;

    parameters.m_hole.m_threshold_eV  = 1.49;
    parameters.m_hole.m_prefactor_s_1 = 1.4e12;
    parameters.m_hole.m_exponent      = 3.4;

    return parameters;
}

std::vector<valley_model> make_silicon_delta_valleys() {
    std::vector<valley_model> valleys;
    valleys.reserve(6);

    valley_model::parameters p;
    p.transverse_effective_mass   = 0.19 * m0;
    p.longitudinal_effective_mass = 0.916 * m0;

    // Keep this simple for now.
    // Use parabolic first; switch to Kane once the loop works.
    p.dispersion = valley_model::band_type::parabolic;

    // If you want Kane immediately, set:
    // p.dispersion = valley_model::band_type::kane;
    // p.non_parabolicity = ...; // in 1/eV or 1/J depending on your convention

    p.non_parabolicity        = 0.0;
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
    p.transverse_effective_mass   = 0.24 * m0;
    p.longitudinal_effective_mass = 0.24 * m0;
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

}  // namespace uepm::amc
