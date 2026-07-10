/**
 * @file pbmc_transport_kernel.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-26
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "pbmc_transport_kernel.hpp"

#include <fmt/core.h>

#include <cmath>
#include <stdexcept>
#include <utility>

#include "pbmc_scattering_model.hpp"
#include "unit_conversion.hpp"

namespace uepm::PBMC {

namespace {

double sample_thermal_energy_eV(double temperature_K, std::mt19937_64& rng) {
    const double                    kT_eV        = uepm::constants::k_b_eV * temperature_K;
    constexpr double                three_halves = 1.5;
    std::gamma_distribution<double> dist(three_halves, kT_eV);
    return dist(rng);
}

std::size_t valley_axis(std::size_t valley_index) {
    if (valley_index < 2) {
        return 0;
    }
    if (valley_index < 4) {
        return 1;
    }
    if (valley_index < 6) {
        return 2;
    }
    throw std::out_of_range("invalid valley index in valley_axis");
}

std::size_t draw_g_destination_valley(std::size_t current_valley_index) {
    switch (current_valley_index) {
        case 0:
            return 1;
        case 1:
            return 0;
        case 2:
            return 3;
        case 3:
            return 2;
        case 4:
            return 5;
        case 5:
            return 4;
        default:
            throw std::out_of_range("invalid valley index in draw_g_destination_valley");
    }
}

std::size_t draw_f_destination_valley(std::size_t current_valley_index, std::mt19937_64& rng) {
    const std::size_t current_axis = valley_axis(current_valley_index);

    std::array<std::size_t, 4> candidates{};
    std::size_t                count = 0;

    for (std::size_t valley_index = 0; valley_index < 6; ++valley_index) {
        if (valley_axis(valley_index) != current_axis) {
            candidates[count++] = valley_index;
        }
    }

    if (count != 4) {
        throw std::runtime_error("unexpected number of f-type destination valleys");
    }

    std::uniform_int_distribution<std::size_t> dist(0, 3);
    return candidates[dist(rng)];
}

std::size_t draw_intervalley_destination_valley(const intervalley_phonon_branch& branch,
                                                std::size_t                      current_valley_index,
                                                std::mt19937_64&                 rng) {
    switch (branch.m_family) {
        case intervalley_family::g:
            return draw_g_destination_valley(current_valley_index);
        case intervalley_family::f:
            return draw_f_destination_valley(current_valley_index, rng);
        default:
            throw std::runtime_error("unknown intervalley family");
    }
}

double impurity_density_cm_3(const pbmc_particle& particle, const pbmc_transport_config& cfg) {
    switch (cfg.m_impurity_density_source) {
        case impurity_density_source::background:
            return std::isfinite(cfg.m_background_impurity_density_cm_3) && cfg.m_background_impurity_density_cm_3 > 0.0
                       ? cfg.m_background_impurity_density_cm_3
                       : 0.0;
        case impurity_density_source::particle_local: {
            const double local_density = particle.state().impurity_concentration_cm_3;
            return std::isfinite(local_density) && local_density > 0.0 ? local_density : 0.0;
        }
    }
    return 0.0;
}

}  // namespace

pbmc_transport_kernel::pbmc_transport_kernel()
    : m_cfg{},
      m_material_model(make_silicon_pbmc_material_model()),
      m_rng(std::random_device{}()) {}

pbmc_transport_kernel::pbmc_transport_kernel(const pbmc_transport_config& cfg)
    : m_cfg(cfg),
      m_material_model(make_silicon_pbmc_material_model()),
      m_rng(std::random_device{}()) {}

pbmc_transport_kernel::pbmc_transport_kernel(const pbmc_transport_config& cfg, std::uint64_t seed)
    : m_cfg(cfg),
      m_material_model(make_silicon_pbmc_material_model()),
      m_rng(seed) {}

pbmc_transport_kernel::pbmc_transport_kernel(const pbmc_transport_config& cfg, pbmc_material_model material)
    : m_cfg(cfg),
      m_material_model(std::move(material)),
      m_rng(std::random_device{}()) {}

pbmc_transport_kernel::pbmc_transport_kernel(const pbmc_transport_config& cfg,
                                             pbmc_material_model          material,
                                             std::uint64_t                seed)
    : m_cfg(cfg),
      m_material_model(std::move(material)),
      m_rng(seed) {}

void pbmc_transport_kernel::initialize() {
    m_material_model.validate();
    if (m_cfg.m_carrier_type == particle_type::electron) {
        m_valleys              = m_material_model.m_electron_valleys;
        m_intervalley_branches = m_material_model.m_electron_intervalley_transitions;
        m_hole_optical_transitions.clear();
    } else {
        m_valleys = m_material_model.m_hole_bands;
        m_intervalley_branches.clear();
        m_hole_optical_transitions = m_material_model.m_hole_optical_transitions;
    }

    m_impact_ionization_parameters = m_material_model.m_impact_ionization;
    m_impurity_mobility_parameters = m_material_model.m_impurity_mobility;

    m_gamma_max_by_valley_s_1.assign(m_valleys.size(), 0.0);
    m_gamma_max_s_1 = 0.0;
    for (std::size_t valley_index = 0; valley_index < m_valleys.size(); ++valley_index) {
        double valley_gamma_max = 0.0;
        for (std::size_t i = 0; i < m_cfg.m_gamma_max_energy_samples; ++i) {
            const double x         = static_cast<double>(i) / static_cast<double>(m_cfg.m_gamma_max_energy_samples - 1);
            const double energy_eV = x * m_cfg.m_max_energy_eV;
            valley_gamma_max =
                std::max(valley_gamma_max,
                         total_scattering_rate_for_energy(valley_index, energy_eV, m_cfg.m_lattice_temperature));
        }
        valley_gamma_max *= m_cfg.m_self_scattering_safety_factor;
        m_gamma_max_by_valley_s_1[valley_index] = valley_gamma_max;
        m_gamma_max_s_1                         = std::max(m_gamma_max_s_1, valley_gamma_max);
    }
}

void pbmc_transport_kernel::initialize_particle_state(pbmc_particle& p) {
    initialize_particle_state(p, m_cfg.m_lattice_temperature);
}

void pbmc_transport_kernel::initialize_particle_state(pbmc_particle& p, double temperature_K) {
    if (m_valleys.empty()) {
        throw std::runtime_error("transport kernel is not initialized");
    }

    if (p.state().valley_index >= m_valleys.size()) {
        throw std::out_of_range("invalid valley index in initialize_particle_state");
    }

    if (!std::isfinite(temperature_K) || temperature_K <= 0.0) {
        throw std::invalid_argument("particle initialization temperature must be positive and finite");
    }

    const auto&  valley    = m_valleys[p.state().valley_index];
    const double energy_eV = sample_thermal_energy_eV(temperature_K, m_rng);

    p.state().local_k        = valley.draw_random_k_valley_at_energy(energy_eV, m_rng);
    p.state().gamma          = valley.gamma_from_k_valley(p.state().local_k);
    p.state().kinetic_energy = valley.kinetic_energy_from_gamma(p.state().gamma);
    p.state().velocity       = valley.to_global_frame(valley.velocity_from_k_valley(p.state().local_k));
}

double pbmc_transport_kernel::uniform01() {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(m_rng);
}

void pbmc_transport_kernel::drift_particle(pbmc_particle& p, const mesh::vector3& electric_field_Vm, double dt) {
    if (dt < 0.0) {
        throw std::invalid_argument("drift time step must be non-negative");
    }

    const auto valley_index = p.state().valley_index;
    if (valley_index >= m_valleys.size()) {
        throw std::out_of_range("invalid valley index in drift_particle");
    }

    const auto&   valley                = m_valleys[valley_index];
    const vector3 electric_field_valley = valley.to_valley_frame(electric_field_Vm);

    const double  prefactor    = p.get_signed_charge() / uepm::constants::h_bar;
    const vector3 old_velocity = p.state().velocity;

    p.state().local_k += electric_field_valley * (prefactor * dt);

    p.state().gamma          = valley.gamma_from_k_valley(p.state().local_k);
    p.state().kinetic_energy = valley.kinetic_energy_from_gamma(p.state().gamma);
    const vector3 velocity_valley =
        valley.velocity_from_k_valley_and_energy(p.state().local_k, p.state().kinetic_energy);
    p.state().velocity = valley.to_global_frame(velocity_valley);

    const vector3 avg_velocity = 0.5 * (old_velocity + p.state().velocity);

    // Note: the mesh is in microns, while the velocity is in m/s.
    p.state().position += avg_velocity * dt * uepm::units::meter_to_micron;
    p.state().time += dt;
}

void pbmc_transport_kernel::set_particle_velocity_direction_preserving_energy(
    pbmc_particle& p,
    const vector3& desired_global_direction) const {
    const auto valley_index = p.state().valley_index;
    if (valley_index >= m_valleys.size()) {
        throw std::out_of_range("invalid valley index in boundary reflection");
    }

    if (desired_global_direction.norm_squared() == 0.0) {
        return;
    }

    const auto&   valley                   = m_valleys[valley_index];
    const vector3 desired_direction_valley = valley.to_valley_frame(desired_global_direction);
    p.state().local_k =
        valley.k_valley_from_energy_velocity_direction(p.state().kinetic_energy, desired_direction_valley);
    p.state().gamma          = valley.gamma_from_k_valley(p.state().local_k);
    p.state().kinetic_energy = valley.kinetic_energy_from_gamma(p.state().gamma);
    p.state().velocity =
        valley.to_global_frame(valley.velocity_from_k_valley_and_energy(p.state().local_k, p.state().kinetic_energy));
}

double pbmc_transport_kernel::sample_free_flight_time() {
    if (m_gamma_max_s_1 <= 0.0) {
        throw std::invalid_argument("max self-scattering rate must be > 0");
    }

    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    double                                 u = 0.0;

    do {
        u = unif01(m_rng);
    } while (u <= 0.0);

    return -std::log(u) / m_gamma_max_s_1;
}

double pbmc_transport_kernel::gamma_max(const pbmc_particle& p) const {
    const auto valley_index = p.state().valley_index;
    if (valley_index >= m_gamma_max_by_valley_s_1.size()) {
        throw std::out_of_range("invalid valley index in gamma_max");
    }
    return m_gamma_max_by_valley_s_1[valley_index];
}

double pbmc_transport_kernel::sample_free_flight_time(const pbmc_particle& p) {
    const double local_gamma_max = gamma_max(p);
    if (local_gamma_max <= 0.0) {
        throw std::invalid_argument("valley max self-scattering rate must be > 0");
    }

    double u = 0.0;
    do {
        u = uniform01();
    } while (u <= 0.0);
    return -std::log(u) / local_gamma_max;
}

const impact_ionization_parameters& pbmc_transport_kernel::impact_ionization_parameters_for_carrier() const {
    if (m_cfg.m_carrier_type == particle_type::electron) {
        return m_impact_ionization_parameters.m_electron;
    }

    return m_impact_ionization_parameters.m_hole;
}

double pbmc_transport_kernel::impact_ionization_rate(double energy_eV) const {
    if (!m_cfg.m_enable_impact_ionization) {
        return 0.0;
    }

    const auto& parameters = impact_ionization_parameters_for_carrier();
    if (parameters.m_threshold_eV <= 0.0) {
        throw std::invalid_argument("Impact ionization threshold must be positive.");
    }
    if (parameters.m_prefactor_s_1 < 0.0) {
        throw std::invalid_argument("Impact ionization prefactor must be non-negative.");
    }
    if (parameters.m_exponent <= 0.0) {
        throw std::invalid_argument("Impact ionization exponent must be positive.");
    }
    if (energy_eV <= parameters.m_threshold_eV) {
        return 0.0;
    }
    const double excess_ratio = (energy_eV - parameters.m_threshold_eV) / parameters.m_threshold_eV;
    return parameters.m_prefactor_s_1 * std::pow(excess_ratio, parameters.m_exponent);
}

scattering_channel pbmc_transport_kernel::select_scattering_channel(const pbmc_particle& p) {
    const auto channels = build_scattering_channels(p);

    double total_rate = 0.0;
    for (const auto& channel : channels) {
        total_rate += channel.rate_s_1;
    }

    return select_scattering_channel(channels, total_rate);
}

scattering_channel pbmc_transport_kernel::select_scattering_channel(const scattering_channel_list& channels,
                                                                    double                         total_rate) {
    if (total_rate <= 0.0) {
        throw std::runtime_error("cannot select scattering channel with zero total rate");
    }

    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    const double                           r_select = unif01(m_rng) * total_rate;

    double cumulative = 0.0;
    for (const auto& channel : channels) {
        cumulative += channel.rate_s_1;
        if (r_select < cumulative) {
            return channel;
        }
    }

    throw std::runtime_error("failed to select a real scattering channel");
}
double pbmc_transport_kernel::impurity_rate_for_particle(const pbmc_particle& p,
                                                         const valley_model&  current_band,
                                                         double               energy_eV) const {
    const double impurity_density = impurity_density_cm_3(p, m_cfg);
    if (impurity_density <= 0.0) {
        return 0.0;
    }
    switch (m_cfg.m_impurity_scattering_model) {
        case impurity_scattering_model::mobility_empirical:
            return impurity_momentum_relaxation_rate(current_band,
                                                     p.type(),
                                                     impurity_density,
                                                     m_impurity_mobility_parameters);
        case impurity_scattering_model::screened_coulomb:
            return screened_coulomb_impurity_momentum_relaxation_rate(current_band,
                                                                      m_material_model.m_dielectric.epsilon_r,
                                                                      energy_eV,
                                                                      impurity_density,
                                                                      impurity_density,
                                                                      p.get_lattice_temperature(),
                                                                      m_cfg.m_impurity_screening_model);
    }
    throw std::runtime_error("unknown impurity scattering model");
}

double pbmc_transport_kernel::impurity_rate_for_energy(const valley_model& band_or_valley,
                                                       double              energy_eV,
                                                       double              temperature_K) const {
    const double impurity_density_cm_3 = m_cfg.m_background_impurity_density_cm_3;
    if (impurity_density_cm_3 <= 0.0) {
        return 0.0;
    }
    switch (m_cfg.m_impurity_scattering_model) {
        case impurity_scattering_model::mobility_empirical:
            return impurity_momentum_relaxation_rate(band_or_valley,
                                                     m_cfg.m_carrier_type,
                                                     impurity_density_cm_3,
                                                     m_impurity_mobility_parameters);

        case impurity_scattering_model::screened_coulomb:
            return screened_coulomb_impurity_momentum_relaxation_rate(band_or_valley,
                                                                      m_material_model.m_dielectric.epsilon_r,
                                                                      energy_eV,
                                                                      impurity_density_cm_3,
                                                                      impurity_density_cm_3,
                                                                      temperature_K,
                                                                      m_cfg.m_impurity_screening_model);
    }

    throw std::runtime_error("unknown impurity scattering model");
}

scattering_channel_list pbmc_transport_kernel::build_scattering_channels(const pbmc_particle& p) const {
    const auto current_band_index = p.state().valley_index;
    if (current_band_index >= m_valleys.size()) {
        throw std::out_of_range("invalid band/valley index in build_scattering_channels");
    }

    const auto&  current_band = m_valleys[current_band_index];
    const double energy_eV    = p.state().kinetic_energy;

    scattering_channel_list channels;

    if (p.type() == particle_type::hole) {
        const double acoustic_rate = acoustic_scattering_rate(current_band,
                                                              m_material_model.m_hole_acoustic,
                                                              energy_eV,
                                                              p.get_lattice_temperature());

        if (acoustic_rate > 0.0) {
            channels.push_back(scattering_channel{.mechanism         = scattering_mechanism::acoustic,
                                                  .rate_s_1          = acoustic_rate,
                                                  .final_energy_eV   = energy_eV,
                                                  .destination_index = current_band_index,
                                                  .branch            = nullptr,
                                                  .process           = intervalley_process::none,
                                                  .transition_name   = "acoustic"});
        }

        for (const auto& transition : m_hole_optical_transitions) {
            if (transition.initial_band != current_band_index) {
                continue;
            }

            const auto& final_band = m_valleys[transition.final_band];

            const double rate_abs =
                optical_scattering_rate_holes(final_band,
                                              transition,
                                              m_material_model.m_hole_acoustic.mass_density_kg_per_m3,
                                              energy_eV,
                                              true,
                                              p.get_lattice_temperature());

            if (rate_abs > 0.0) {
                channels.push_back(scattering_channel{.mechanism         = scattering_mechanism::intervalley,
                                                      .rate_s_1          = rate_abs,
                                                      .final_energy_eV   = energy_eV + transition.phonon_energy_eV,
                                                      .destination_index = transition.final_band,
                                                      .branch            = nullptr,
                                                      .process           = intervalley_process::absorption,
                                                      .transition_name   = transition.name});
            }

            const double final_energy_emission_eV = energy_eV - transition.phonon_energy_eV;
            const double rate_em =
                optical_scattering_rate_holes(final_band,
                                              transition,
                                              m_material_model.m_hole_acoustic.mass_density_kg_per_m3,
                                              energy_eV,
                                              false,
                                              p.get_lattice_temperature());

            if (rate_em > 0.0 && final_energy_emission_eV >= 0.0) {
                channels.push_back(scattering_channel{.mechanism         = scattering_mechanism::intervalley,
                                                      .rate_s_1          = rate_em,
                                                      .final_energy_eV   = final_energy_emission_eV,
                                                      .destination_index = transition.final_band,
                                                      .branch            = nullptr,
                                                      .process           = intervalley_process::emission,
                                                      .transition_name   = transition.name});
            }
        }

        if (m_cfg.m_enable_impurity_scattering) {
            const double rate_impurity = impurity_rate_for_particle(p, current_band, energy_eV);

            if (rate_impurity > 0.0) {
                channels.push_back(scattering_channel{
                    .mechanism         = scattering_mechanism::impurity,
                    .rate_s_1          = rate_impurity,
                    .final_energy_eV   = energy_eV,
                    .destination_index = current_band_index,
                    .transition_name   = "impurity",
                });
            }
        }

        if (m_cfg.m_enable_impact_ionization) {
            const double impact_rate = impact_ionization_rate(energy_eV);

            if (impact_rate > 0.0) {
                const auto& impact_parameters = impact_ionization_parameters_for_carrier();

                channels.push_back(scattering_channel{.mechanism         = scattering_mechanism::impact_ionization,
                                                      .rate_s_1          = impact_rate,
                                                      .final_energy_eV   = energy_eV - impact_parameters.m_threshold_eV,
                                                      .destination_index = current_band_index,
                                                      .branch            = nullptr,
                                                      .process           = intervalley_process::none,
                                                      .transition_name   = "impact_ionization"});
            }
        }

        return channels;
    }

    const double acoustic_rate = acoustic_scattering_rate(current_band,
                                                          m_material_model.m_electron_acoustic,
                                                          energy_eV,
                                                          p.get_lattice_temperature());

    if (acoustic_rate > 0.0) {
        channels.push_back(scattering_channel{.mechanism       = scattering_mechanism::acoustic,
                                              .rate_s_1        = acoustic_rate,
                                              .final_energy_eV = energy_eV,
                                              .branch          = nullptr,
                                              .process         = intervalley_process::none,
                                              .transition_name = "acoustic"});
    }
    const double optical_boost = 1.0;
    for (const auto& branch : m_intervalley_branches) {
        const double rate_abs = intervalley_scattering_rate(current_band,
                                                            branch,
                                                            m_material_model.m_electron_acoustic.mass_density_kg_per_m3,
                                                            energy_eV,
                                                            true,
                                                            p.get_lattice_temperature());

        if (rate_abs > 0.0) {
            channels.push_back(scattering_channel{.mechanism       = scattering_mechanism::intervalley,
                                                  .rate_s_1        = rate_abs * optical_boost,
                                                  .final_energy_eV = energy_eV + branch.m_phonon_energy_eV,
                                                  .branch          = &branch,
                                                  .process         = intervalley_process::absorption,
                                                  .transition_name = branch.m_name});
        }

        const double rate_em = intervalley_scattering_rate(current_band,
                                                           branch,
                                                           m_material_model.m_electron_acoustic.mass_density_kg_per_m3,
                                                           energy_eV,
                                                           false,
                                                           p.get_lattice_temperature());

        const double final_energy_emission_eV = energy_eV - branch.m_phonon_energy_eV;
        if (rate_em > 0.0 && final_energy_emission_eV >= 0.0) {
            channels.push_back(scattering_channel{.mechanism       = scattering_mechanism::intervalley,
                                                  .rate_s_1        = rate_em * optical_boost,
                                                  .final_energy_eV = final_energy_emission_eV,
                                                  .branch          = &branch,
                                                  .process         = intervalley_process::emission,
                                                  .transition_name = branch.m_name});
        }
    }

    if (m_cfg.m_enable_impurity_scattering) {
        const double rate_impurity = impurity_rate_for_particle(p, current_band, energy_eV);

        if (rate_impurity > 0.0) {
            channels.push_back(scattering_channel{
                .mechanism         = scattering_mechanism::impurity,
                .rate_s_1          = rate_impurity,
                .final_energy_eV   = energy_eV,
                .destination_index = current_band_index,
                .transition_name   = "impurity",
            });
        }
    }

    if (m_cfg.m_enable_impact_ionization) {
        const double impact_rate = impact_ionization_rate(energy_eV);

        if (impact_rate > 0.0) {
            const auto& impact_parameters = impact_ionization_parameters_for_carrier();

            channels.push_back(scattering_channel{.mechanism         = scattering_mechanism::impact_ionization,
                                                  .rate_s_1          = impact_rate,
                                                  .final_energy_eV   = energy_eV - impact_parameters.m_threshold_eV,
                                                  .destination_index = current_band_index,
                                                  .branch            = nullptr,
                                                  .process           = intervalley_process::none,
                                                  .transition_name   = "impact_ionization"});
        }
    }

    return channels;
}

double pbmc_transport_kernel::total_scattering_rate(const pbmc_particle& p) const {
    const auto current_band_index = p.state().valley_index;
    if (current_band_index >= m_valleys.size()) {
        throw std::out_of_range("invalid band/valley index in total_scattering_rate");
    }

    const auto&  current_band = m_valleys[current_band_index];
    const double energy_eV    = p.state().kinetic_energy;
    const double temperature  = p.get_lattice_temperature();
    double       total_rate   = 0.0;

    if (p.type() == particle_type::hole) {
        total_rate += acoustic_scattering_rate(current_band, m_material_model.m_hole_acoustic, energy_eV, temperature);

        for (const auto& transition : m_hole_optical_transitions) {
            if (transition.initial_band != current_band_index) {
                continue;
            }

            const auto& final_band = m_valleys[transition.final_band];
            total_rate += optical_scattering_rate_holes(final_band,
                                                        transition,
                                                        m_material_model.m_hole_acoustic.mass_density_kg_per_m3,
                                                        energy_eV,
                                                        true,
                                                        temperature);

            if (energy_eV >= transition.phonon_energy_eV) {
                total_rate += optical_scattering_rate_holes(final_band,
                                                            transition,
                                                            m_material_model.m_hole_acoustic.mass_density_kg_per_m3,
                                                            energy_eV,
                                                            false,
                                                            temperature);
            }
        }

        if (m_cfg.m_enable_impurity_scattering) {
            total_rate += impurity_rate_for_particle(p, current_band, energy_eV);
        }
        if (m_cfg.m_enable_impact_ionization) {
            total_rate += impact_ionization_rate(energy_eV);
        }
        return total_rate;
    }

    total_rate += acoustic_scattering_rate(current_band, m_material_model.m_electron_acoustic, energy_eV, temperature);

    for (const auto& branch : m_intervalley_branches) {
        total_rate += intervalley_scattering_rate(current_band,
                                                  branch,
                                                  m_material_model.m_electron_acoustic.mass_density_kg_per_m3,
                                                  energy_eV,
                                                  true,
                                                  temperature);

        if (energy_eV >= branch.m_phonon_energy_eV) {
            total_rate += intervalley_scattering_rate(current_band,
                                                      branch,
                                                      m_material_model.m_electron_acoustic.mass_density_kg_per_m3,
                                                      energy_eV,
                                                      false,
                                                      temperature);
        }
    }

    if (m_cfg.m_enable_impurity_scattering) {
        total_rate += impurity_rate_for_particle(p, current_band, energy_eV);
    }
    if (m_cfg.m_enable_impact_ionization) {
        total_rate += impact_ionization_rate(energy_eV);
    }
    return total_rate;
}

double pbmc_transport_kernel::total_scattering_rate_for_energy(std::size_t band_or_valley_index,
                                                               double      energy_eV,
                                                               double      max_temperature_K) const {
    if (band_or_valley_index >= m_valleys.size()) {
        throw std::out_of_range("invalid band/valley index in total_scattering_rate_for_energy");
    }

    const auto& band_or_valley = m_valleys[band_or_valley_index];
    double      total_rate     = 0.0;

    if (m_cfg.m_carrier_type == particle_type::hole) {
        total_rate +=
            acoustic_scattering_rate(band_or_valley, m_material_model.m_hole_acoustic, energy_eV, max_temperature_K);

        for (const auto& transition : m_hole_optical_transitions) {
            if (transition.initial_band != band_or_valley_index) {
                continue;
            }

            const auto& final_band = m_valleys[transition.final_band];

            total_rate += optical_scattering_rate_holes(final_band,
                                                        transition,
                                                        m_material_model.m_hole_acoustic.mass_density_kg_per_m3,
                                                        energy_eV,
                                                        true,
                                                        max_temperature_K);

            total_rate += optical_scattering_rate_holes(final_band,
                                                        transition,
                                                        m_material_model.m_hole_acoustic.mass_density_kg_per_m3,
                                                        energy_eV,
                                                        false,
                                                        max_temperature_K);
        }

        if (m_cfg.m_enable_impurity_scattering && m_cfg.m_background_impurity_density_cm_3 > 0.0) {
            total_rate += impurity_rate_for_energy(band_or_valley, energy_eV, max_temperature_K);
        }

        if (m_cfg.m_enable_impact_ionization) {
            total_rate += impact_ionization_rate(energy_eV);
        }

        return total_rate;
    }

    total_rate +=
        acoustic_scattering_rate(band_or_valley, m_material_model.m_electron_acoustic, energy_eV, max_temperature_K);

    for (const auto& branch : m_intervalley_branches) {
        total_rate += intervalley_scattering_rate(band_or_valley,
                                                  branch,
                                                  m_material_model.m_electron_acoustic.mass_density_kg_per_m3,
                                                  energy_eV,
                                                  true,
                                                  max_temperature_K);
        total_rate += intervalley_scattering_rate(band_or_valley,
                                                  branch,
                                                  m_material_model.m_electron_acoustic.mass_density_kg_per_m3,
                                                  energy_eV,
                                                  false,
                                                  max_temperature_K);
    }

    if (m_cfg.m_enable_impurity_scattering && m_cfg.m_background_impurity_density_cm_3 > 0.0) {
        total_rate += impurity_rate_for_energy(band_or_valley, energy_eV, max_temperature_K);
    }
    if (m_cfg.m_enable_impact_ionization) {
        total_rate += impact_ionization_rate(energy_eV);
    }
    return total_rate;
}

void pbmc_transport_kernel::ensure_gamma_max_covers(double total_rate) {
    if (total_rate <= m_gamma_max_s_1) {
        return;
    }

    const double old_gamma_max = m_gamma_max_s_1;
    m_gamma_max_s_1            = total_rate * m_cfg.m_self_scattering_safety_factor;
    for (double& valley_gamma_max : m_gamma_max_by_valley_s_1) {
        valley_gamma_max = std::max(valley_gamma_max, m_gamma_max_s_1);
    }

    fmt::print(stderr,
               "Warning: gamma_max increased from {:.6e} to {:.6e} s^-1 "
               "after observing total scattering rate {:.6e} s^-1. "
               "Consider increasing --max-energy or --gamma-safety.\n",
               old_gamma_max,
               m_gamma_max_s_1,
               total_rate);
}

void pbmc_transport_kernel::ensure_gamma_max_covers(std::size_t band_or_valley_index, double total_rate) {
    if (band_or_valley_index >= m_gamma_max_by_valley_s_1.size()) {
        throw std::out_of_range("invalid valley index in ensure_gamma_max_covers");
    }
    if (total_rate <= m_gamma_max_by_valley_s_1[band_or_valley_index]) {
        return;
    }

    m_gamma_max_by_valley_s_1[band_or_valley_index] = total_rate * m_cfg.m_self_scattering_safety_factor;
    m_gamma_max_s_1 = std::max(m_gamma_max_s_1, m_gamma_max_by_valley_s_1[band_or_valley_index]);
}

double pbmc_transport_kernel::compute_max_self_scattering_rate(double      max_energy_eV,
                                                               double      max_temperature_K,
                                                               std::size_t n_samples) const {
    if (max_energy_eV <= 0.0) {
        throw std::invalid_argument("max energy must be > 0");
    }

    if (n_samples < 2) {
        throw std::invalid_argument("number of gamma-max energy samples must be >= 2");
    }

    double gamma_max = 0.0;

    for (std::size_t valley_index = 0; valley_index < m_valleys.size(); ++valley_index) {
        for (std::size_t i = 0; i < n_samples; ++i) {
            const double x         = static_cast<double>(i) / static_cast<double>(n_samples - 1);
            const double energy_eV = x * max_energy_eV;

            const double total_rate = total_scattering_rate_for_energy(valley_index, energy_eV, max_temperature_K);

            gamma_max = std::max(gamma_max, total_rate);
        }
    }

    return gamma_max * m_cfg.m_self_scattering_safety_factor;
}

scattering_event pbmc_transport_kernel::apply_scattering_channel(pbmc_particle& p, const scattering_channel& channel) {
    if (channel.rate_s_1 < 0.0) {
        throw std::invalid_argument("negative scattering channel rate");
    }

    switch (channel.mechanism) {
        case scattering_mechanism::acoustic: {
            const auto band_or_valley_index = p.state().valley_index;
            if (band_or_valley_index >= m_valleys.size()) {
                throw std::out_of_range("invalid band/valley index in apply_scattering_channel acoustic");
            }

            const auto& band_or_valley = m_valleys[band_or_valley_index];

            p.state().local_k        = band_or_valley.draw_random_k_valley_at_energy(channel.final_energy_eV, m_rng);
            p.state().gamma          = band_or_valley.gamma_from_k_valley(p.state().local_k);
            p.state().kinetic_energy = band_or_valley.kinetic_energy_from_gamma(p.state().gamma);
            p.state().velocity       = band_or_valley.to_global_frame(
                band_or_valley.velocity_from_k_valley_and_energy(p.state().local_k, p.state().kinetic_energy));

            p.increment_scattering_event_count();
            p.add_scattering_event(scattering_event::acoustic);
            return scattering_event::acoustic;
        }

        case scattering_mechanism::intervalley: {
            if (p.type() == particle_type::hole) {
                if (channel.destination_index >= m_valleys.size()) {
                    throw std::out_of_range("invalid destination band in apply_scattering_channel for holes");
                }
                const auto& dst_band     = m_valleys[channel.destination_index];
                p.state().valley_index   = channel.destination_index;
                p.state().local_k        = dst_band.draw_random_k_valley_at_energy(channel.final_energy_eV, m_rng);
                p.state().gamma          = dst_band.gamma_from_k_valley(p.state().local_k);
                p.state().kinetic_energy = dst_band.kinetic_energy_from_gamma(p.state().gamma);
                p.state().velocity       = dst_band.to_global_frame(
                    dst_band.velocity_from_k_valley_and_energy(p.state().local_k, p.state().kinetic_energy));
                p.increment_scattering_event_count();
                p.add_transition_event(channel.transition_name);
                if (channel.process == intervalley_process::absorption) {
                    p.add_scattering_event(scattering_event::intervalley_absorption);
                    return scattering_event::intervalley_absorption;
                }
                if (channel.process == intervalley_process::emission) {
                    p.add_scattering_event(scattering_event::intervalley_emission);
                    return scattering_event::intervalley_emission;
                }
                throw std::runtime_error("hole optical channel missing absorption/emission tag");
            }

            if (channel.branch == nullptr) {
                throw std::runtime_error("electron intervalley channel missing branch");
            }
            const auto current_valley_index = p.state().valley_index;
            if (current_valley_index >= m_valleys.size()) {
                throw std::out_of_range("invalid current valley index in apply_scattering_channel for electrons");
            }
            const std::size_t destination_valley =
                draw_intervalley_destination_valley(*channel.branch, current_valley_index, m_rng);
            if (destination_valley >= m_valleys.size()) {
                throw std::out_of_range("invalid destination valley in apply_scattering_channel for electrons");
            }

            const auto& dst_valley   = m_valleys[destination_valley];
            p.state().valley_index   = destination_valley;
            p.state().local_k        = dst_valley.draw_random_k_valley_at_energy(channel.final_energy_eV, m_rng);
            p.state().gamma          = dst_valley.gamma_from_k_valley(p.state().local_k);
            p.state().kinetic_energy = dst_valley.kinetic_energy_from_gamma(p.state().gamma);
            p.state().velocity       = dst_valley.to_global_frame(
                dst_valley.velocity_from_k_valley_and_energy(p.state().local_k, p.state().kinetic_energy));
            p.increment_scattering_event_count();
            p.add_transition_event(channel.transition_name);
            if (channel.process == intervalley_process::absorption) {
                p.add_scattering_event(scattering_event::intervalley_absorption);
                return scattering_event::intervalley_absorption;
            }
            if (channel.process == intervalley_process::emission) {
                p.add_scattering_event(scattering_event::intervalley_emission);
                return scattering_event::intervalley_emission;
            }
            throw std::runtime_error("electron intervalley channel missing absorption/emission tag");
        }
        case scattering_mechanism::impurity: {
            const auto band_or_valley_index = p.state().valley_index;

            if (band_or_valley_index >= m_valleys.size()) {
                throw std::out_of_range("invalid band/valley index in apply_scattering_channel impurity");
            }
            const auto& band_or_valley = m_valleys[band_or_valley_index];
            // Elastic impurity scattering:
            // same band/valley, same kinetic energy, randomized direction.
            p.state().local_k        = band_or_valley.draw_random_k_valley_at_energy(channel.final_energy_eV, m_rng);
            p.state().gamma          = band_or_valley.gamma_from_k_valley(p.state().local_k);
            p.state().kinetic_energy = band_or_valley.kinetic_energy_from_gamma(p.state().gamma);
            p.state().velocity       = band_or_valley.to_global_frame(
                band_or_valley.velocity_from_k_valley_and_energy(p.state().local_k, p.state().kinetic_energy));
            p.increment_scattering_event_count();
            p.add_scattering_event(scattering_event::impurity);

            return scattering_event::impurity;
        }
        case scattering_mechanism::impact_ionization: {
            const auto band_or_valley_index = p.state().valley_index;
            if (band_or_valley_index >= m_valleys.size()) {
                throw std::out_of_range("invalid band/valley index in apply_scattering_channel impact ionization");
            }
            const auto&  band_or_valley  = m_valleys[band_or_valley_index];
            const double final_energy_eV = std::max(0.0, channel.final_energy_eV);
            p.state().local_k            = band_or_valley.draw_random_k_valley_at_energy(final_energy_eV, m_rng);
            p.state().gamma              = band_or_valley.gamma_from_k_valley(p.state().local_k);
            p.state().kinetic_energy     = band_or_valley.kinetic_energy_from_gamma(p.state().gamma);
            p.state().velocity           = band_or_valley.to_global_frame(
                band_or_valley.velocity_from_k_valley_and_energy(p.state().local_k, p.state().kinetic_energy));
            p.increment_scattering_event_count();
            p.add_scattering_event(scattering_event::impact_ionization);
            return scattering_event::impact_ionization;
        }
    }

    throw std::runtime_error("unknown scattering mechanism");
}

std::optional<scattering_event> pbmc_transport_kernel::scatter_particle(pbmc_particle& p, double dt) {
    if (dt < 0.0) {
        throw std::invalid_argument("scatter time step must be non-negative");
    }

    const double total_rate = total_scattering_rate(p);

    if (total_rate <= 0.0) {
        return std::nullopt;
    }

    const double scatter_probability = 1.0 - std::exp(-total_rate * dt);

    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    if (unif01(m_rng) >= scatter_probability) {
        return std::nullopt;
    }

    const auto   channels = build_scattering_channels(p);
    const double r_select = unif01(m_rng) * total_rate;

    double cumulative = 0.0;
    for (const auto& channel : channels) {
        cumulative += channel.rate_s_1;
        if (r_select < cumulative) {
            auto event = apply_scattering_channel(p, channel);
            return event;
        }
    }

    throw std::runtime_error("failed to select a scattering channel");
}

}  // namespace uepm::PBMC
