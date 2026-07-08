/**
 * @file impact_ionization_matrix_element.hpp
 * @brief Plane-wave two-body matrix elements for ab initio impact ionization.
 */

#pragma once

#include <Eigen/Dense>

#include <array>
#include <complex>
#include <functional>
#include <vector>

#include "vector_bz.hpp"

namespace uepm::mesh_bz {

using complex_d = std::complex<double>;
using impact_screened_interaction = std::function<complex_d(const vector3& q_SI, double energy_transfer_eV)>;

struct ImpactIonizationPlaneWaveState {
    vector3          k_SI{};
    double           energy_eV = 0.0;
    Eigen::VectorXcd coefficients{};
};

struct ImpactIonizationMatrixElementConfig {
    double momentum_tolerance_SI = 1.0e6;
    double normalization_volume_m3 = 1.0;
};

complex_d compute_screened_two_body_matrix_element(
    const ImpactIonizationPlaneWaveState& initial_a,
    const ImpactIonizationPlaneWaveState& initial_b,
    const ImpactIonizationPlaneWaveState& final_a,
    const ImpactIonizationPlaneWaveState& final_b,
    const std::vector<vector3>&           basis_vectors_SI,
    const impact_screened_interaction&    screened_interaction,
    const ImpactIonizationMatrixElementConfig& config = {});

std::array<complex_d, 2> compute_direct_exchange_impact_ionization_matrix_element(
    const ImpactIonizationPlaneWaveState& initial_hot_electron,
    const ImpactIonizationPlaneWaveState& initial_valence_electron,
    const ImpactIonizationPlaneWaveState& final_electron_1,
    const ImpactIonizationPlaneWaveState& final_electron_2,
    const std::vector<vector3>&           basis_vectors_SI,
    const impact_screened_interaction&    screened_interaction,
    const ImpactIonizationMatrixElementConfig& config = {});

double antisymmetrized_impact_ionization_strength(complex_d direct, complex_d exchange);

}  // namespace uepm::mesh_bz
