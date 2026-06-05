/**
 * @file amc_scattering_model.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-26
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include "amc_material_model.hpp"
#include "intervalley_phonon.hpp"
#include "particle_amc.hpp"
#include "valley_model.hpp"

namespace uepm::amc {

double bose_einstein_occupation(double phonon_energy_eV, double temperature_K);

double acoustic_scattering_rate_silicon(const valley_model& valley, double energy_eV, double temperature_K);
double acoustic_scattering_rate_silicon_holes(const valley_model& band, double energy_eV, double temperature_K);
double silicon_density_of_states_mass(const valley_model& valley);

double silicon_nonparabolicity_per_joule(const valley_model& valley);

double intervalley_zeroth_order_rate(const valley_model&              valley,
                                     const intervalley_phonon_branch& branch,
                                     double                           initial_energy_eV,
                                     bool                             absorption,
                                     double                           temperature_K);

double intervalley_first_order_rate(const valley_model&              valley,
                                    const intervalley_phonon_branch& branch,
                                    double                           initial_energy_eV,
                                    bool                             absorption,
                                    double                           temperature_K);

double intervalley_scattering_rate(const valley_model&              valley,
                                   const intervalley_phonon_branch& branch,
                                   double                           initial_energy_eV,
                                   bool                             absorption,
                                   double                           temperature_K);

double optical_scattering_rate_silicon_holes(const valley_model&            final_band,
                                             const hole_optical_transition& transition,
                                             double                         initial_energy_eV,
                                             bool                           absorption,
                                             double                         temperature_K);

double caughey_thomas_mobility_cm2_per_V_s(double                              impurity_density_cm_3,
                                           const impurity_mobility_parameters& parameters);

double impurity_momentum_relaxation_rate_silicon(const valley_model&                         band_or_valley,
                                                 particle_type                               carrier_type,
                                                 double                                      impurity_density_cm_3,
                                                 const carrier_impurity_mobility_parameters& parameters);

double screened_coulomb_impurity_momentum_relaxation_rate_silicon(const valley_model& band_or_valley,
                                                                  particle_type       carrier_type,
                                                                  double              energy_eV,
                                                                  double              impurity_density_cm_3,
                                                                  double              screening_density_cm_3,
                                                                  double              temperature_K);

};  // namespace uepm::amc