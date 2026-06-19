/**
 * @file single_part_fbmc.cpp
 * @author
 * @brief Single-particle FBMC with null-collision (Option A)
 * @version 0.1
 * @date 2025-09-19
 */

#include "single_part_fbmc.hpp"

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

#include "keldysh_impactio.hpp"
#include "physical_constants.hpp"

namespace uepm::fbmc {

namespace {

double sample_thermal_energy_eV(double temperature_K, std::mt19937& rng) {
    const double                    kT_eV = uepm::constants::k_b_eV * temperature_K;
    std::gamma_distribution<double> distribution(1.5, kT_eV);
    return distribution(rng);
}

std::uint32_t particle_seed(std::uint64_t base_seed, std::size_t particle_index) {
    std::uint64_t value = base_seed + 0x9e3779b97f4a7c15ULL * (particle_index + 1);
    value               = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value               = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    value ^= value >> 31U;
    return static_cast<std::uint32_t>(value ^ (value >> 32U));
}

double sampled_time_after_warmup(double time_before_s, double time_after_s, double warmup_time_s, double final_time_s) {
    const double sample_start = std::max(time_before_s, warmup_time_s);
    const double sample_stop  = std::min(time_after_s, final_time_s);
    return std::max(0.0, sample_stop - sample_start);
}

uepm::mesh_bz::MeshParticleType mesh_particle_type(particle_type type) {
    return type == particle_type::electron ? uepm::mesh_bz::MeshParticleType::conduction
                                           : uepm::mesh_bz::MeshParticleType::valence;
}

std::pair<vector3, std::size_t> draw_carrier_state_at_energy(uepm::mesh_bz::ElectronPhonon& mesh,
                                                             particle_type                    type,
                                                             double                           energy_eV,
                                                             std::mt19937&                    rng) {
    std::vector<std::size_t> candidate_bands;
    std::vector<double>      weights;
    for (const std::size_t band : mesh.get_band_indices(mesh_particle_type(type))) {
        const auto [minimum, maximum] = mesh.get_min_max_energy_at_band(static_cast<int>(band));
        if (energy_eV < minimum || energy_eV > maximum) {
            continue;
        }
        const double dos = mesh.compute_dos_at_energy_and_band(energy_eV, static_cast<int>(band));
        if (dos > 0.0 && std::isfinite(dos)) {
            candidate_bands.push_back(band);
            weights.push_back(dos);
        }
    }
    if (candidate_bands.empty()) {
        throw std::runtime_error("No carrier band has non-zero DOS at the requested energy");
    }
    std::discrete_distribution<std::size_t> band_distribution(weights.begin(), weights.end());
    const std::size_t band = candidate_bands[band_distribution(rng)];
    return {mesh.draw_random_k_point_at_energy(energy_eV, band, rng), band};
}

}  // namespace

double impact_ionization_coefficient_statistics::event_rate_per_carrier_s_1() const {
    return m_carrier_time_s > 0.0 ? static_cast<double>(m_events) / m_carrier_time_s : 0.0;
}

double impact_ionization_coefficient_statistics::average_drift_velocity_m_per_s() const {
    return m_carrier_time_s > 0.0 ? m_drift_velocity_time_integral_m / m_carrier_time_s : 0.0;
}

double impact_ionization_coefficient_statistics::ionization_coefficient_cm_1() const {
    const double drift_velocity_cm_per_s = average_drift_velocity_m_per_s() * 1.0e2;
    return drift_velocity_cm_per_s > 0.0 ? event_rate_per_carrier_s_1() / drift_velocity_cm_per_s : 0.0;
}

Single_particle_simulation::Single_particle_simulation(uepm::mesh_bz::ElectronPhonon* ptr_mesh_bz,
                                                       const Bulk_environment&        bulk_env,
                                                       const Simulation_parameters&   sim_params,
                                                       std::size_t                    nb_particles)
    : m_ptr_mesh_bz(ptr_mesh_bz),
      m_bulk_env(bulk_env),
      m_sim_params(sim_params),
      m_nb_particles(nb_particles) {
    if (m_ptr_mesh_bz == nullptr) {
        throw std::invalid_argument("Single_particle_simulation: BZ mesh pointer must not be null");
    }
    if (!(m_bulk_env.m_temperature > 0.0) || !std::isfinite(m_bulk_env.m_temperature)) {
        throw std::invalid_argument("Single_particle_simulation: temperature must be finite and positive");
    }
    if (!(m_sim_params.m_simulation_time > 0.0) || !std::isfinite(m_sim_params.m_simulation_time)) {
        throw std::invalid_argument("Single_particle_simulation: simulation time must be finite and positive");
    }
    if (!(m_sim_params.m_warmup_fraction >= 0.0 && m_sim_params.m_warmup_fraction < 1.0) ||
        !std::isfinite(m_sim_params.m_warmup_fraction)) {
        throw std::invalid_argument("Single_particle_simulation: warmup fraction must be finite and in [0, 1)");
    }
    if (m_sim_params.m_nb_openmp_threads == 0) {
        throw std::invalid_argument("Single_particle_simulation: thread count must be positive");
    }
    if (!(m_sim_params.m_max_energy_eV > 0.0) || !std::isfinite(m_sim_params.m_max_energy_eV)) {
        throw std::invalid_argument("Single_particle_simulation: maximum energy must be finite and positive");
    }
    if (!(m_sim_params.m_self_scattering_safety_factor >= 1.0) ||
        !std::isfinite(m_sim_params.m_self_scattering_safety_factor)) {
        throw std::invalid_argument(
            "Single_particle_simulation: self-scattering safety factor must be finite and at least one");
    }
    if (m_nb_particles == 0) {
        throw std::invalid_argument("Single_particle_simulation: particle count must be positive");
    }
    if (m_ptr_mesh_bz->get_particle_type() != mesh_particle_type(m_sim_params.m_particle_type)) {
        throw std::invalid_argument("Single_particle_simulation: mesh and FBMC carrier types do not match");
    }

    m_list_particle.reserve(m_nb_particles);
    for (std::size_t i = 0; i < m_nb_particles; ++i) {
        particle new_particle(i, m_sim_params.m_particle_type, m_ptr_mesh_bz);
        new_particle.seed_random_generator(particle_seed(m_sim_params.m_random_seed, i));
        new_particle.state().m_position = {0.0, 0.0, 0.0};
        m_list_particle.emplace_back(std::move(new_particle));
    }

    for (auto& current_particle : m_list_particle) {
        bool initialized = false;
        for (std::size_t attempt = 0; attempt < 10000 && !initialized; ++attempt) {
            const double thermal_energy =
                sample_thermal_energy_eV(m_bulk_env.m_temperature, current_particle.get_random_generator());
            if (thermal_energy > m_sim_params.m_max_energy_eV) {
                continue;
            }
            try {
                auto [initial_k, initial_band] = draw_carrier_state_at_energy(*m_ptr_mesh_bz,
                                                                              current_particle.get_type(),
                                                                              thermal_energy,
                                                                              current_particle.get_random_generator());
                current_particle.state().m_k_vector   = initial_k;
                current_particle.state().m_band_index = static_cast<int>(initial_band);
                initialized                           = true;
            } catch (const std::runtime_error&) {
                // The thermal distribution may sample an energy outside every loaded band.
            }
        }
        if (!initialized) {
            throw std::runtime_error("Could not draw an initial thermal state inside the loaded FBMC bands");
        }

        uepm::mesh_bz::Tetra* containing_tetra =
            m_ptr_mesh_bz->find_tetra_at_location(current_particle.state().m_k_vector);
        if (containing_tetra == nullptr) {
            throw std::runtime_error("Initial k-point is out of the Brillouin zone mesh.");
        }
        current_particle.set_containing_bz_mesh_tetra(containing_tetra);
        current_particle.update_energy();
        current_particle.update_group_velocity();
        const auto rates =
            current_particle.interpolate_phonon_scattering_rate_at_location(current_particle.state().m_k_vector);
        current_particle.state().m_gamma = std::accumulate(rates.begin(), rates.end(), 0.0);
        if (m_sim_params.m_record_history) {
            current_particle.update_history();
        }
    }

    fmt::print("Initialized {} thermal FBMC particles at {:.1f} K with seed {}\n",
               m_list_particle.size(),
               m_bulk_env.m_temperature,
               m_sim_params.m_random_seed);
}

Single_particle_simulation::Single_particle_simulation(uepm::mesh_bz::ElectronPhonon*     ptr_mesh_bz,
                                                       const bulk_fbmc_simulation_config& config)
    : Single_particle_simulation(
          ptr_mesh_bz,
          Bulk_environment{config.m_lattice_temperature, config.m_electric_field, config.m_doping_concentration_cm_3},
          Simulation_parameters{.m_simulation_time               = config.m_final_time,
                                .m_warmup_fraction               = config.m_warmup_fraction,
                                .m_nb_openmp_threads             = config.m_nb_threads,
                                .m_max_energy_eV                 = config.m_max_energy_eV,
                                .m_self_scattering_safety_factor = config.m_self_scattering_safety_factor,
                                .m_enable_impact_ionization      = config.m_enable_impact_ionization,
                                .m_record_history                = config.m_record_history,
                                .m_random_seed                   = config.m_random_seed,
                                .m_history_export_prefix         = config.m_history_export_prefix,
                                .m_particle_type                 = config.m_particle_type},
          config.m_number_of_particles) {}

void Single_particle_simulation::run_simulation() {
    if (m_sim_params.m_particle_type == particle_type::hole && m_sim_params.m_enable_impact_ionization) {
        throw std::invalid_argument("FBMC hole impact ionization is not implemented");
    }
    const double p_gamma_elph = m_ptr_mesh_bz->compute_P_Gamma();
    if (!(p_gamma_elph > 0.0) || !std::isfinite(p_gamma_elph)) {
        throw std::runtime_error("Invalid maximum carrier-phonon scattering rate");
    }

    KeldyshImpactIonization keldysh_impactio;
    const double            p_gamma_impact =
        m_sim_params.m_enable_impact_ionization ? keldysh_impactio.compute_rate(m_sim_params.m_max_energy_eV) : 0.0;
    if (!(p_gamma_impact >= 0.0) || !std::isfinite(p_gamma_impact)) {
        throw std::runtime_error("Invalid maximum impact-ionization scattering rate");
    }

    m_gamma_max_s_1 = (p_gamma_elph + p_gamma_impact) * m_sim_params.m_self_scattering_safety_factor;
    if (!(m_gamma_max_s_1 > 0.0) || !std::isfinite(m_gamma_max_s_1)) {
        throw std::runtime_error("Invalid total self-scattering rate");
    }

    fmt::print("Using gamma_max = {:.3e} s^-1 (el-ph: {:.3e}, impact: {:.3e}, safety: {:.3f})\n",
               m_gamma_max_s_1,
               p_gamma_elph,
               p_gamma_impact,
               m_sim_params.m_self_scattering_safety_factor);

    const double  final_time_s  = m_sim_params.m_simulation_time;
    const double  warmup_time_s = m_sim_params.m_warmup_fraction * final_time_s;
    const double  field_norm    = m_bulk_env.m_electric_field.norm();
    const vector3 field_direction =
        field_norm > 0.0 ? m_bulk_env.m_electric_field / field_norm : vector3{0.0, 0.0, 0.0};

    m_observables                          = {};
    m_observables.m_electric_field_V_per_m = field_norm;
    m_impact_ionization_statistics         = {};

    double                   reduced_weighted_velocity_x  = 0.0;
    double                   reduced_weighted_velocity_y  = 0.0;
    double                   reduced_weighted_velocity_z  = 0.0;
    double                   reduced_weighted_energy      = 0.0;
    double                   reduced_accumulated_time     = 0.0;
    std::size_t              reduced_ii_events            = 0;
    double                   reduced_ii_carrier_time      = 0.0;
    double                   reduced_ii_velocity_integral = 0.0;
    double                   max_observed_total_rate      = 0.0;
    std::exception_ptr       parallel_exception;
    std::atomic<std::size_t> completed_particles{0};
    std::atomic<int>         last_reported_progress_bucket{0};
    const auto               simulation_start = std::chrono::steady_clock::now();

    fmt::print("FBMC progress: 0% (0/{})\n", m_list_particle.size());

#pragma omp parallel for schedule(dynamic) num_threads(m_sim_params.m_nb_openmp_threads) \
    reduction(+ : reduced_weighted_velocity_x,                                           \
                  reduced_weighted_velocity_y,                                           \
                  reduced_weighted_velocity_z,                                           \
                  reduced_weighted_energy,                                               \
                  reduced_accumulated_time,                                              \
                  reduced_ii_events,                                                     \
                  reduced_ii_carrier_time,                                               \
                  reduced_ii_velocity_integral) reduction(max : max_observed_total_rate)
    for (std::size_t idx = 0; idx < m_list_particle.size(); ++idx) {
        try {
            auto&                                  current_particle = m_list_particle[idx];
            std::uniform_real_distribution<double> U01(0.0, 1.0);

            while (current_particle.state().m_time < final_time_s) {
                const double time_before_flight = current_particle.state().m_time;
                current_particle.draw_free_flight_time(m_gamma_max_s_1);

                const double  sampled_dt             = current_particle.state().m_free_flight_time;
                const double  remaining_time         = final_time_s - current_particle.state().m_time;
                const double  dt                     = std::min(sampled_dt, remaining_time);
                const vector3 velocity_before_flight = current_particle.state().m_velocity;

                current_particle.state().m_free_flight_time = dt;
                current_particle.update_k_vector(m_bulk_env.m_electric_field);

                uepm::mesh_bz::Tetra* containing_tetra =
                    m_ptr_mesh_bz->find_tetra_at_location(current_particle.state().m_k_vector);
                if (containing_tetra == nullptr) {
                    current_particle.state().m_k_vector =
                        m_ptr_mesh_bz->retrieve_k_inside_mesh_geometry(current_particle.state().m_k_vector);
                    containing_tetra = m_ptr_mesh_bz->find_tetra_at_location(current_particle.state().m_k_vector);
                    if (containing_tetra == nullptr) {
                        throw std::runtime_error("Particle k-point is out of the Brillouin zone mesh after folding.");
                    }
                }
                current_particle.set_containing_bz_mesh_tetra(containing_tetra);

                current_particle.update_energy();
                current_particle.update_group_velocity();

                const vector3 velocity_after_flight = current_particle.state().m_velocity;
                current_particle.update_position(0.5 * (velocity_before_flight + velocity_after_flight), dt);
                current_particle.advance_time_by_free_flight();

                const double sampled_dt_after_warmup = sampled_time_after_warmup(time_before_flight,
                                                                                 current_particle.state().m_time,
                                                                                 warmup_time_s,
                                                                                 final_time_s);
                if (sampled_dt_after_warmup > 0.0) {
                    reduced_weighted_velocity_x += current_particle.state().m_velocity.x() * sampled_dt_after_warmup;
                    reduced_weighted_velocity_y += current_particle.state().m_velocity.y() * sampled_dt_after_warmup;
                    reduced_weighted_velocity_z += current_particle.state().m_velocity.z() * sampled_dt_after_warmup;
                    reduced_weighted_energy += current_particle.state().m_energy * sampled_dt_after_warmup;
                    reduced_accumulated_time += sampled_dt_after_warmup;
                    reduced_ii_carrier_time += sampled_dt_after_warmup;
                    if (field_norm > 0.0) {
                        reduced_ii_velocity_integral +=
                            std::abs(current_particle.state().m_velocity.dot(field_direction)) *
                            sampled_dt_after_warmup;
                    }
                }
                if (m_sim_params.m_record_history) {
                    current_particle.update_history();
                }

                if (sampled_dt > remaining_time) {
                    break;
                }

                std::array<double, 8> rates_elph = current_particle.interpolate_phonon_scattering_rate_at_location(
                    current_particle.state().m_k_vector);
                const double Gamma_elph = std::accumulate(rates_elph.begin(), rates_elph.end(), 0.0);

                const double energy_particle = current_particle.state().m_energy;
                const double rate_impactio =
                    m_sim_params.m_enable_impact_ionization ? keldysh_impactio.compute_rate(energy_particle) : 0.0;
                const double Gamma               = Gamma_elph + rate_impactio;
                max_observed_total_rate          = std::max(max_observed_total_rate, Gamma);
                current_particle.state().m_gamma = Gamma;

                if (!(Gamma > 0.0) || !std::isfinite(Gamma)) {
                    current_particle.add_scattering_event_to_history(9);
                    if (m_sim_params.m_record_history) {
                        current_particle.update_history();
                    }
                    continue;
                }

                double accept = Gamma / m_gamma_max_s_1;

                constexpr double eps = 1e-12;
                if (accept > 1.0 + eps) {
                    throw std::runtime_error(
                        "FBMC self-scattering bound violated. Increase --maxenergy or --gamma-safety.");
                }
                accept = std::min(accept, 1.0);
                if (U01(current_particle.get_random_generator()) > accept) {
                    current_particle.add_scattering_event_to_history(9);
                    if (m_sim_params.m_record_history) {
                        current_particle.update_history();
                    }
                    continue;
                }

                const double rsel = U01(current_particle.get_random_generator()) * Gamma;
                if (rsel > Gamma_elph) {
                    current_particle.select_final_state_after_impact_ionization(keldysh_impactio.m_E_threshold);
                    current_particle.add_scattering_event_to_history(8);
                    if (current_particle.state().m_time >= warmup_time_s) {
                        ++reduced_ii_events;
                    }
                } else {
                    double      cum       = 0.0;
                    std::size_t event_idx = 0;
                    for (; event_idx < rates_elph.size(); ++event_idx) {
                        cum += rates_elph[event_idx];
                        if (rsel <= cum) {
                            break;
                        }
                    }
                    if (event_idx >= rates_elph.size()) {
                        event_idx = rates_elph.size() - 1;
                    }

                    current_particle.select_final_state_after_phonon_scattering(event_idx);
                    current_particle.add_scattering_event_to_history(event_idx);
                }
                if (m_sim_params.m_record_history) {
                    current_particle.update_history();
                }
            }

            if (m_sim_params.m_record_history && !m_sim_params.m_history_export_prefix.empty()) {
                const std::string history_filename = fmt::format("{}_particle_{}.csv",
                                                                 m_sim_params.m_history_export_prefix,
                                                                 current_particle.get_index());
                current_particle.export_history_to_csv(history_filename);
                current_particle.reset_history();
            }

            const std::size_t completed        = completed_particles.fetch_add(1) + 1;
            const int         progress_percent = static_cast<int>((100 * completed) / m_list_particle.size());
            const int         progress_bucket  = std::min(10, progress_percent / 10);
            int               previous_bucket  = last_reported_progress_bucket.load();
            if (progress_bucket > previous_bucket &&
                last_reported_progress_bucket.compare_exchange_strong(previous_bucket, progress_bucket)) {
                const double elapsed_s =
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - simulation_start).count();
                fmt::print("FBMC progress: {}% ({}/{}, {:.1f} s)\n",
                           progress_bucket * 10,
                           completed,
                           m_list_particle.size(),
                           elapsed_s);
            }
        } catch (...) {
#pragma omp critical(fbmc_parallel_exception)
            {
                if (parallel_exception == nullptr) {
                    parallel_exception = std::current_exception();
                }
            }
        }
    }

    if (parallel_exception != nullptr) {
        std::rethrow_exception(parallel_exception);
    }

    m_observables.m_weighted_velocity_x_m                           = reduced_weighted_velocity_x;
    m_observables.m_weighted_velocity_y_m                           = reduced_weighted_velocity_y;
    m_observables.m_weighted_velocity_z_m                           = reduced_weighted_velocity_z;
    m_observables.m_weighted_kinetic_energy_eV_s                    = reduced_weighted_energy;
    m_observables.m_accumulated_time_s                              = reduced_accumulated_time;
    m_impact_ionization_statistics.m_events                         = reduced_ii_events;
    m_impact_ionization_statistics.m_carrier_time_s                 = reduced_ii_carrier_time;
    m_impact_ionization_statistics.m_drift_velocity_time_integral_m = reduced_ii_velocity_integral;

    if (!(m_observables.m_accumulated_time_s > 0.0)) {
        throw std::runtime_error("No FBMC steady-state samples were accumulated");
    }

    fmt::print("Completed self-scattering FBMC run with {} particles\n", m_list_particle.size());
    fmt::print("Maximum observed total scattering rate: {:.6e} s^-1\n", max_observed_total_rate);
    fmt::print("Steady-state average velocity: ({:.6e}, {:.6e}, {:.6e}) m/s\n",
               m_observables.m_weighted_velocity_x_m / m_observables.m_accumulated_time_s,
               m_observables.m_weighted_velocity_y_m / m_observables.m_accumulated_time_s,
               m_observables.m_weighted_velocity_z_m / m_observables.m_accumulated_time_s);
    fmt::print("Steady-state average energy: {:.6f} eV\n",
               m_observables.m_weighted_kinetic_energy_eV_s / m_observables.m_accumulated_time_s);
    fmt::print("Impact ionization coefficient: {:.6e} cm^-1\n",
               m_impact_ionization_statistics.ionization_coefficient_cm_1());
}

void Single_particle_simulation::extract_stats_and_export(const std::string& filename) {
    export_observables_to_csv(filename);
}

void Single_particle_simulation::export_observables_to_csv(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error(fmt::format("Could not open {} for writing", filename));
    }
    if (!(m_observables.m_accumulated_time_s > 0.0)) {
        throw std::runtime_error("No FBMC observables have been accumulated");
    }

    file << "charge_C,"
            "temperature_K,"
            "electric_field_V_per_m,"
            "impurity_density_cm_3,"
            "mean_velocity_x_m_per_s,"
            "mean_velocity_y_m_per_s,"
            "mean_velocity_z_m_per_s,"
            "mean_kinetic_energy_eV,"
            "sample_count,"
            "impact_ionization_events,"
            "impact_ionization_rate_per_carrier_s_1,"
            "impact_ionization_drift_velocity_m_per_s,"
            "impact_ionization_coefficient_cm_1\n";

    const double carrier_charge_C =
        m_sim_params.m_particle_type == particle_type::electron ? -uepm::constants::q_e : uepm::constants::q_e;
    file << fmt::format("{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
                        carrier_charge_C,
                        m_bulk_env.m_temperature,
                        m_observables.m_electric_field_V_per_m,
                        m_bulk_env.m_doping_concentration,
                        m_observables.m_weighted_velocity_x_m / m_observables.m_accumulated_time_s,
                        m_observables.m_weighted_velocity_y_m / m_observables.m_accumulated_time_s,
                        m_observables.m_weighted_velocity_z_m / m_observables.m_accumulated_time_s,
                        m_observables.m_weighted_kinetic_energy_eV_s / m_observables.m_accumulated_time_s,
                        m_observables.m_accumulated_time_s,
                        m_impact_ionization_statistics.m_events,
                        m_impact_ionization_statistics.event_rate_per_carrier_s_1(),
                        m_impact_ionization_statistics.average_drift_velocity_m_per_s(),
                        m_impact_ionization_statistics.ionization_coefficient_cm_1());
    fmt::print("Exported FBMC observables to {}\n", filename);
}

void Single_particle_simulation::export_history(const std::string& filename) {
    if (!m_sim_params.m_record_history) {
        throw std::runtime_error("FBMC history recording is disabled");
    }
    for (auto& particle : m_list_particle) {
        std::string filename_particle = fmt::format("{}_particle_{}.csv", filename, particle.get_index());
        particle.export_history_to_csv(filename_particle);
        fmt::print("Exported history of particle {} to {}\n", particle.get_index(), filename_particle);
    }
}

}  // namespace uepm::fbmc
