/**
 * @file pbmc_device_runner.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-09
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "pbmc_device_runner.hpp"

#include <fmt/chrono.h>
#include <fmt/core.h>
#include <fmt/format.h>

#include <chrono>
#include <filesystem>
#include <stdexcept>
#include <thread>

#include "pbmc_device_setup.hpp"
#include "pbmc_run_manifest.hpp"
#include "device.hpp"
#include "materials.hpp"
#include "msh_file.hpp"

namespace uepm::PBMC {

void run_self_consistent_device_pbmc_simulation(const self_consistent_device_pbmc_run_config& config) {
    validate_material_symbol(config.material_symbol);
    auto device_options               = config.device_options;
    device_options.m_simulation_name  = config.simulation_name;
    device_options.m_output_directory = config.output_dir;

    auto self_consistent_options_2d = config.self_consistent_options_2d;
    auto self_consistent_options_3d = config.self_consistent_options_3d;

    device_options.validate();
    self_consistent_options_2d.validate();

    const std::string output_dir =
        config.output_dir.empty() ? make_default_output_directory(config.mesh_file) : config.output_dir;
    const std::string trajectory_dir  = fmt::format("{}/trajectory", output_dir);
    device_options.m_output_directory = output_dir;

    std::filesystem::create_directories(output_dir);
    if (config.device_options.m_export_time_step || config.device_options.m_keep_particles_history) {
        std::filesystem::create_directories(trajectory_dir);
    }

    fmt::print("Loading mesh: {}\n", config.mesh_file);

    uepm::file::msh_file msh_file(config.mesh_file);
    msh_file.read_mesh();
    msh_file.read_states();

    uepm::mesh::mesh* mesh = msh_file.get_p_mesh();
    if (mesh == nullptr) {
        throw std::runtime_error("Mesh loading failed.");
    }

    const int mesh_dimension = mesh->get_dimension();
    if (mesh_dimension != 2 && mesh_dimension != 3) {
        throw std::runtime_error("Only 2D and 3D meshes are supported.");
    }

    fmt::print("Mesh dimension: {}D\n", mesh_dimension);
    const uepm::physics::material_repository material_repository =
        config.material_root.empty() ? uepm::physics::material_repository{}
                                     : uepm::physics::material_repository{config.material_root};
    fmt::print("Loading materials from repository: {}\n", material_repository.root().string());

    uepm::physics::material_database material_database = material_repository.load_all_materials();
    const auto&                      common_material   = material_database.require(config.material_symbol);
    device_options.m_material_model                    = load_pbmc_material_model(material_repository, common_material);

    uepm::device::device simulation_device(mesh);
    add_collecting_contacts(simulation_device, *mesh, config.collecting_contacts);

    const std::string material_symbol = config.material_symbol;

    fmt::print("Self-consistent PBMC {}D simulation\n", mesh_dimension);
    fmt::print("  mesh vertices: {}\n", mesh->get_nb_vertices());
    fmt::print("  output directory: {}\n", output_dir);
    fmt::print("  material: {}\n", material_symbol);
    fmt::print("  final time: {:.6e} s\n", config.device_options.m_t_max);
    fmt::print("  time step: {:.6e} s\n", config.device_options.m_time_step);
    fmt::print("  Poisson frequency: {}\n", config.self_consistent_options_2d.m_common.m_poisson_frequency);
    fmt::print("  contact voltages:\n");
    for (const auto& [contact_name, voltage_V] :
         config.self_consistent_options_2d.m_common.m_contact_voltages_V) {
        fmt::print("    {}: {:.6e} V\n", contact_name, voltage_V);
    }
    fmt::print("  Ramo electrode: {}\n", config.self_consistent_options_2d.m_common.m_ramo_electrode);
    fmt::print("  collecting contacts:");
    for (const auto& contact_name : config.collecting_contacts) {
        fmt::print(" {}", contact_name);
    }
    fmt::print("\n");
    fmt::print("  built-in potential: {}\n",
               config.self_consistent_options_2d.m_common.m_enable_built_in_potential ? "enabled" : "disabled");
    if (config.self_consistent_options_2d.m_common.m_enable_built_in_potential) {
        fmt::print("  intrinsic concentration: {:.6e} cm^-3\n",
                   config.self_consistent_options_2d.m_common.m_intrinsic_concentration_cm_3);
    }
    fmt::print("  export time steps: {}\n", config.device_options.m_export_time_step ? "enabled" : "disabled");
    if (mesh_dimension == 2) {
        fmt::print("  effective depth: {:.6e} um\n", config.self_consistent_options_2d.m_effective_depth_um);
        fmt::print("  particle z period: {:.6e} um\n", config.self_consistent_options_2d.m_particle_z_period_um);
    }
    fmt::print("  initial electrons: {}\n", config.number_electrons_start);
    fmt::print("  initial holes: {}\n", config.number_holes_start);
    fmt::print("  initial position: ({:.6e}, {:.6e}, {:.6e})\n",
               config.starting_position.x(),
               config.starting_position.y(),
               config.starting_position.z());

    const auto& common_options =
        mesh_dimension == 2 ? self_consistent_options_2d.m_common : self_consistent_options_3d.m_common;
    simulation_manifest manifest;
    manifest.add("run", "simulation_type", "self_consistent_device_PBMC");
    manifest.add("run", "simulation_name", config.simulation_name);
    manifest.add("run", "status", "completed");
    manifest.add("run", "started_at_utc", current_utc_timestamp());
    manifest.add("run", "command_line", config.command_line);
    manifest.add("run", "working_directory", std::filesystem::current_path().string());

    manifest.add("build", "project_version", pbmc_project_version());
    manifest.add("build", "build_type", pbmc_build_type());
    manifest.add("build", "compiler", pbmc_compiler());
    manifest.add("build", "hardware_concurrency", static_cast<std::size_t>(std::thread::hardware_concurrency()));

    manifest.add("input", "device_mesh", std::filesystem::absolute(config.mesh_file).string());
    manifest.add("input", "material_root", material_repository.root().string());
    manifest.add("input", "material_file", material_repository.material_file(config.material_symbol).string());
    manifest.add("input", "material", config.material_symbol);
    manifest.add("input", "output_directory", std::filesystem::absolute(output_dir).string());
    manifest.add("input", "mesh_dimension", mesh_dimension);
    manifest.add("input", "mesh_vertices", static_cast<std::size_t>(mesh->get_nb_vertices()));
    manifest.add("input", "mesh_regions", static_cast<std::size_t>(mesh->get_nb_regions()));

    manifest.add("simulation", "requested_threads", device_options.m_nb_threads);
    manifest.add("simulation", "random_seed", config.seed_random_generator);
    manifest.add("simulation", "final_time_s", device_options.m_t_max);
    manifest.add("simulation", "time_step_s", device_options.m_time_step);
    manifest.add("simulation", "lattice_temperature_K", device_options.m_lattice_temperature);
    manifest.add("simulation", "initial_electrons", config.number_electrons_start);
    manifest.add("simulation", "initial_holes", config.number_holes_start);
    manifest.add("simulation", "initial_x_um", config.starting_position.x());
    manifest.add("simulation", "initial_y_um", config.starting_position.y());
    manifest.add("simulation", "initial_z_um", config.starting_position.z());
    manifest.add("simulation", "max_particles", device_options.m_max_number_particle);
    manifest.add("simulation", "stop_when_no_electrons", device_options.m_stop_simu_when_no_electron_remaining);
    manifest.add("simulation", "keep_particle_history", device_options.m_keep_particles_history);
    manifest.add("simulation", "export_time_steps", device_options.m_export_time_step);
    manifest.add("simulation", "export_frequency", device_options.m_frequency_export_trajectory);
    if (mesh_dimension == 2) {
        manifest.add("simulation", "effective_depth_um", self_consistent_options_2d.m_effective_depth_um);
        manifest.add("simulation", "particle_z_period_um", self_consistent_options_2d.m_particle_z_period_um);
    }

    manifest.add("transport", "impact_ionization_enabled", device_options.m_activate_impact_ionization);
    manifest.add("transport", "particle_creation_enabled", device_options.m_particle_creation_activated);
    manifest.add("transport", "impurity_scattering_enabled", device_options.m_enable_impurity_scattering);
    manifest.add("transport",
                 "impurity_model",
                 device_options.m_impurity_scattering_model == impurity_scattering_model::mobility_empirical
                     ? "mobility-empirical"
                     : "screened-coulomb");
    manifest.add("transport",
                 "impurity_screening",
                 impurity_screening_model_name(device_options.m_impurity_screening_model));
    manifest.add("transport", "gamma_max_energy_eV", device_options.m_max_energy_eV);
    manifest.add("transport", "gamma_safety_factor", device_options.m_self_scattering_safety_factor);
    manifest.add("transport", "gamma_energy_samples", device_options.m_gamma_max_energy_samples);

    manifest.add("self_consistent", "poisson_frequency", common_options.m_poisson_frequency);
    manifest.add("self_consistent", "frozen_field_mode", common_options.m_frozen_field_mode);
    manifest.add("self_consistent", "ramo_electrode", common_options.m_ramo_electrode);
    for (const auto& [contact_name, voltage_V] : common_options.m_contact_voltages_V) {
        manifest.add("contact_voltages_V", contact_name, voltage_V);
    }
    for (const auto& contact_name : config.collecting_contacts) {
        manifest.add("collecting_contacts", contact_name, true);
    }
    manifest.add("self_consistent", "built_in_potential_enabled", common_options.m_enable_built_in_potential);
    manifest.add("self_consistent",
                 "intrinsic_concentration_cm_3",
                 common_options.m_intrinsic_concentration_cm_3);
    manifest.add("self_consistent",
                 "built_in_contact_voltage_scale",
                 common_options.m_built_in_contact_voltage_scale);
    manifest.add("self_consistent",
                 "initialize_particles_from_doping",
                 common_options.m_initialize_particles_from_doping);
    manifest.add("self_consistent", "initial_particle_weight", common_options.m_initial_particle_weight);
    manifest.add("self_consistent",
                 "contact_injection_particle_weight",
                 common_options.m_contact_injection_particle_weight);
    manifest.add("self_consistent", "background_ramo_current_A", common_options.m_background_ramo_current_A);
    manifest.add("self_consistent",
                 "auto_background_ramo_current",
                 common_options.m_auto_background_ramo_current);

    const auto& quench_options = common_options.m_passive_quench_circuit;
    manifest.add("quench_circuit", "enabled", quench_options.m_enabled);
    manifest.add("quench_circuit",
                 "biased_contact",
                 common_options.m_quench_biased_contact);
    manifest.add("quench_circuit", "bias_voltage_V", quench_options.m_bias_voltage_V);
    manifest.add("quench_circuit", "initial_device_voltage_V", quench_options.m_initial_device_voltage_V);
    manifest.add("quench_circuit", "resistance_ohm", quench_options.m_resistance_ohm);
    manifest.add("quench_circuit", "capacitance_F", quench_options.m_capacitance_F);
    manifest.add("quench_circuit",
                 "ramo_current_to_quench_current_sign",
                 common_options.m_ramo_current_to_quench_current_sign);
    manifest.add("avalanche_detection",
                 "voltage_drop_threshold_V",
                 common_options.m_avalanche_voltage_drop_threshold_V);
    manifest.add("successful_quench_detection",
                 "high_field_threshold_V_per_cm",
                 common_options.m_quench_high_field_threshold_V_per_cm);
    manifest.add("successful_quench_detection", "quiet_time_s", common_options.m_quench_quiet_time_s);

    const auto& injection = device_options.m_scheduled_particle_injection;
    manifest.add("scheduled_injection", "enabled", device_options.m_enable_scheduled_particle_injection);
    manifest.add("scheduled_injection", "time_s", injection.m_time_s);
    manifest.add("scheduled_injection", "x_um", injection.m_position_um.x());
    manifest.add("scheduled_injection", "y_um", injection.m_position_um.y());
    manifest.add("scheduled_injection", "z_um", injection.m_position_um.z());
    manifest.add("scheduled_injection",
                 "particle_type",
                 injection.m_particle_type == particle_type::electron ? "electron" : "hole");
    manifest.add("scheduled_injection", "weight", injection.m_weight);

    std::size_t remaining_electrons        = 0;
    std::size_t remaining_holes            = 0;
    std::size_t impact_events              = 0;
    double      final_time_s               = 0.0;
    double      total_electron_weight      = 0.0;
    double      total_hole_weight          = 0.0;
    double      final_ramo_current_A       = 0.0;
    bool        avalanche_detected         = false;
    double      avalanche_time_s           = 0.0;
    double      avalanche_voltage_drop_V   = 0.0;
    bool        successful_quench_detected = false;
    double      successful_quench_time_s   = 0.0;

    const auto start = std::chrono::high_resolution_clock::now();

    if (mesh_dimension == 2) {
        self_consistent_device_pbmc_simulation_2d simulation(simulation_device,
                                                            device_options,
                                                            config.self_consistent_options_2d,
                                                            material_database,
                                                            config.starting_position,
                                                            config.number_electrons_start,
                                                            config.number_holes_start,
                                                            config.seed_random_generator);

        simulation.set_prefix_export_trajectory_filename(trajectory_dir);
        simulation.run_self_consistent_transport_simulation();

        const std::string history_file = fmt::format("{}/device_history.csv", output_dir);
        simulation.export_history_to_csv(history_file);
        fmt::print("Wrote {}\n", history_file);
        if (device_options.m_keep_particles_history) {
            const std::string trajectory_prefix = fmt::format("{}/particle_", trajectory_dir);
            simulation.export_all_trajectories_as_csv(trajectory_prefix);
            fmt::print("Wrote particle trajectories to {}\n", trajectory_dir);
        }
        fmt::print("Remaining electrons: {}\n", simulation.get_number_electrons());
        fmt::print("Remaining holes: {}\n", simulation.get_number_holes());
        remaining_electrons   = simulation.get_number_electrons();
        remaining_holes       = simulation.get_number_holes();
        total_electron_weight = simulation.get_total_electron_weight();
        total_hole_weight     = simulation.get_total_hole_weight();
        final_time_s          = simulation.get_current_time().value_or(0.0);
        const auto& history   = simulation.get_simulation_history();
        impact_events         = history.m_impact_ionization_positions.size();
        if (!history.m_list_ramo_current.empty()) {
            final_ramo_current_A = history.m_list_ramo_current.back();
        }
        const auto& avalanche = simulation.avalanche_detection();
        avalanche_detected    = avalanche.m_detected;
        if (avalanche.m_time_s.has_value()) {
            avalanche_time_s = *avalanche.m_time_s;
        }
        if (avalanche.m_voltage_drop_V.has_value()) {
            avalanche_voltage_drop_V = *avalanche.m_voltage_drop_V;
        }
        const auto& successful_quench = simulation.successful_quench_detection();
        successful_quench_detected    = successful_quench.m_detected;
        if (successful_quench.m_time_s.has_value()) {
            successful_quench_time_s = *successful_quench.m_time_s;
        }
    } else {
        self_consistent_device_pbmc_simulation_3d simulation(simulation_device,
                                                            device_options,
                                                            self_consistent_options_3d,
                                                            material_database,
                                                            config.starting_position,
                                                            config.number_electrons_start,
                                                            config.number_holes_start,
                                                            config.seed_random_generator);

        simulation.set_prefix_export_trajectory_filename(fmt::format("{}/time_step", trajectory_dir));
        simulation.run_self_consistent_transport_simulation();

        const std::string history_file = fmt::format("{}/device_history.csv", output_dir);
        simulation.export_history_to_csv(history_file);
        fmt::print("Wrote {}\n", history_file);
        if (device_options.m_keep_particles_history) {
            const std::string trajectory_prefix = fmt::format("{}/particle_", trajectory_dir);
            simulation.export_all_trajectories_as_csv(trajectory_prefix);
            fmt::print("Wrote particle trajectories to {}\n", trajectory_dir);
        }
        fmt::print("Remaining electrons: {}\n", simulation.get_number_electrons());
        fmt::print("Remaining holes: {}\n", simulation.get_number_holes());
        remaining_electrons   = simulation.get_number_electrons();
        remaining_holes       = simulation.get_number_holes();
        total_electron_weight = simulation.get_total_electron_weight();
        total_hole_weight     = simulation.get_total_hole_weight();
        final_time_s          = simulation.get_current_time().value_or(0.0);
        const auto& history   = simulation.get_simulation_history();
        impact_events         = history.m_impact_ionization_positions.size();
        if (!history.m_list_ramo_current.empty()) {
            final_ramo_current_A = history.m_list_ramo_current.back();
        }
        const auto& avalanche = simulation.avalanche_detection();
        avalanche_detected    = avalanche.m_detected;
        if (avalanche.m_time_s.has_value()) {
            avalanche_time_s = *avalanche.m_time_s;
        }
        if (avalanche.m_voltage_drop_V.has_value()) {
            avalanche_voltage_drop_V = *avalanche.m_voltage_drop_V;
        }
        const auto& successful_quench = simulation.successful_quench_detection();
        successful_quench_detected    = successful_quench.m_detected;
        if (successful_quench.m_time_s.has_value()) {
            successful_quench_time_s = *successful_quench.m_time_s;
        }
        if (device_options.m_enable_scheduled_particle_injection) {
            const auto& injection = device_options.m_scheduled_particle_injection;
            fmt::print("  scheduled injection: enabled\n");
            fmt::print("    time: {:.6e} s\n", injection.m_time_s);
            fmt::print("    position: ({:.6e}, {:.6e}, {:.6e}) um\n",
                       injection.m_position_um.x(),
                       injection.m_position_um.y(),
                       injection.m_position_um.z());
            fmt::print("    type: {}\n", injection.m_particle_type == particle_type::electron ? "electron" : "hole");
            fmt::print("    weight: {:.6e}\n", injection.m_weight);
        }
    }

    const auto                          stop    = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> elapsed = stop - start;
    fmt::print("Simulation completed in {:.3f} s\n", elapsed.count());

    manifest.add("run", "finished_at_utc", current_utc_timestamp());
    manifest.add("run", "elapsed_seconds", elapsed.count());
    manifest.add("results", "final_time_s", final_time_s);
    manifest.add("results", "remaining_electrons", remaining_electrons);
    manifest.add("results", "remaining_holes", remaining_holes);
    manifest.add("results", "total_electron_weight", total_electron_weight);
    manifest.add("results", "total_hole_weight", total_hole_weight);
    manifest.add("results", "impact_ionization_events", impact_events);
    manifest.add("results", "final_ramo_current_A", final_ramo_current_A);
    manifest.add("results", "avalanche_detected", avalanche_detected);
    if (avalanche_detected) {
        manifest.add("results", "avalanche_time_s", avalanche_time_s);
        manifest.add("results", "avalanche_voltage_drop_V", avalanche_voltage_drop_V);
    }
    manifest.add("results", "successful_quench_detected", successful_quench_detected);
    if (successful_quench_detected) {
        manifest.add("results", "successful_quench_time_s", successful_quench_time_s);
        manifest.add("results", "avalanche_to_successful_quench_s", successful_quench_time_s - avalanche_time_s);
    }
    manifest.add("outputs", "device_history_csv", fmt::format("{}/device_history.csv", output_dir));
    manifest.add("outputs", "trajectory_directory", trajectory_dir);
    manifest.add("outputs", "particle_trajectories_exported", device_options.m_keep_particles_history);
    manifest.add("outputs", "time_steps_exported", device_options.m_export_time_step);

    const std::filesystem::path manifest_file = std::filesystem::path(output_dir) / "simulation_manifest.txt";
    manifest.add("outputs", "simulation_manifest", manifest_file.string());
    manifest.write(manifest_file);
    fmt::print("Wrote {}\n", manifest_file.string());
}

}  // namespace uepm::PBMC
