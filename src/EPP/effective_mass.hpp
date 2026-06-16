#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "Vector3D.h"

namespace uepm::pseudopotential {

enum class band_edge_kind { minimum, maximum };

struct valley_fit_sample {
    Vector3D<double> k_reduced;
    double           energy_eV = 0.0;
};

struct effective_mass_fit_result {
    Vector3D<double> k0_reduced;
    double           edge_energy_eV          = 0.0;
    band_edge_kind   edge_kind               = band_edge_kind::minimum;
    double           non_parabolicity_eV_inv = 0.0;
    double           rms_error_meV           = 0.0;
    double           mass_rms_error_meV      = 0.0;
    double           alpha_rms_error_meV     = 0.0;
    std::size_t      sample_count            = 0;
    std::size_t      mass_sample_count       = 0;
    std::size_t      alpha_sample_count      = 0;

    std::array<std::array<double, 3>, 3> hessian_eV_m2{};
    std::array<std::array<double, 3>, 3> principal_axes{};
    std::array<double, 3>                principal_masses_m0{};
};

effective_mass_fit_result fit_effective_mass_and_nonparabolicity(
    const std::vector<valley_fit_sample>& samples,
    const Vector3D<double>&               k0_reduced,
    double                                edge_energy_eV,
    double                                lattice_constant_m,
    band_edge_kind                        edge_kind);

effective_mass_fit_result fit_effective_mass_tensor(const std::vector<valley_fit_sample>& samples,
                                                    const Vector3D<double>&               k0_reduced,
                                                    double                                edge_energy_eV,
                                                    double                                lattice_constant_m,
                                                    band_edge_kind                        edge_kind);

effective_mass_fit_result fit_nonparabolicity_with_fixed_mass(
    const std::vector<valley_fit_sample>& samples,
    const effective_mass_fit_result&      mass_fit,
    double                                lattice_constant_m,
    double                                max_kinetic_energy_eV,
    bool                                  clamp_nonnegative);

effective_mass_fit_result fit_effective_mass_then_nonparabolicity(
    const std::vector<valley_fit_sample>& mass_samples,
    const std::vector<valley_fit_sample>& alpha_samples,
    const Vector3D<double>&               k0_reduced,
    double                                edge_energy_eV,
    double                                lattice_constant_m,
    band_edge_kind                        edge_kind,
    double                                max_kinetic_energy_eV,
    bool                                  clamp_nonnegative);

std::string to_string(band_edge_kind edge_kind);
band_edge_kind band_edge_kind_from_string(const std::string& value);

}  // namespace uepm::pseudopotential
