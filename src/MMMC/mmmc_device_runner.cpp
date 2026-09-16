/**
 * @file mmmc_device_runner.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-07-10
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "mmmc_device_runner.hpp"

#include <fmt/core.h>

#include <chrono>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <thread>

#include "materials.hpp"
#include "msh_file.hpp"
#include "pbmc_device_setup.hpp"
#include "pbmc_material_model.hpp"
#include "pbmc_run_manifest.hpp"

namespace uepm::MMMC {
namespace {

std::string quote_command_argument(std::string_view argument) {
    if (argument.find_first_of(" \t\"'\\") == std::string_view::npos) {
        return std::string(argument);
    }
    std::string quoted = "\"";
    for (const char character : argument) {
        if (character == '"' || character == '\\') {
            quoted.push_back('\\');
        }
        quoted.push_back(character);
    }
    quoted.push_back('"');
    return quoted;
}

std::string make_default_mmmc_output_directory(const std::string& mesh_file) {
    return fmt::format("self_consistent_mmmc_{}", std::filesystem::path(mesh_file).stem().string());
}

}  // namespace

std::string mmmc_command_line_from_arguments(int argc, const char* const* argv) {
    std::string command_line;
    for (int index = 0; index < argc; ++index) {
        if (!command_line.empty()) {
            command_line.push_back(' ');
        }
        command_line += quote_command_argument(argv[index]);
    }
    return command_line;
}

void run_self_consistent_device_mmmc_simulation(const self_consistent_device_mmmc_run_config& config) {
    PBMC::validate_material_symbol(config.material_symbol);

    const std::string output_dir =
        config.output_dir.empty() ? make_default_mmmc_output_directory(config.mesh_file) : config.output_dir;
    const std::string trajectory_dir = (std::filesystem::path(output_dir) / "trajectory").string();
    std::filesystem::create_directories(output_dir);
    if (config.device_options.m_pbmc.m_export_time_step || config.device_options.m_pbmc.m_keep_particles_history) {
        std::filesystem::create_directories(trajectory_dir);
    }

    fmt::print("Loading mesh: {}\n", config.mesh_file);
    file::msh_file msh_file(config.mesh_file);
    msh_file.read_mesh();
    msh_file.read_states();
    auto* mesh = msh_file.get_p_mesh();
    if (mesh == nullptr) {
        throw std::runtime_error("MMMC mesh loading failed.");
    }
    if (mesh->get_dimension() != 2) {
        throw std::runtime_error("MMMC self-consistent device mode currently supports 2D meshes.");
    }

    const physics::material_repository material_repository = config.material_root.empty()
                                                                 ? physics::material_repository{}
                                                                 : physics::material_repository{config.material_root};
    physics::material_database         material_database   = material_repository.load_all_materials();
    const auto&                        common_material     = material_database.require(config.material_symbol);

    auto device_options                      = config.device_options;
    device_options.m_pbmc.m_simulation_name  = config.simulation_name;
    device_options.m_pbmc.m_output_directory = output_dir;
    device_options.m_pbmc.m_material_model   = PBMC::load_pbmc_material_model(material_repository, common_material);
    device_options.synchronize_from_pbmc();
    device_options.validate();

    device::device simulation_device(mesh);
    PBMC::add_collecting_contacts(simulation_device, *mesh, config.collecting_contacts);

    fmt::print("Self-consistent MMMC 2D simulation\n");
    fmt::print("  output directory: {}\n", output_dir);
    fmt::print("  material: {}\n", config.material_symbol);
    fmt::print("  final time: {:.6e} s\n", device_options.m_pbmc.m_t_max);
    fmt::print("  time step: {:.6e} s\n", device_options.m_pbmc.m_time_step);
    fmt::print("  Poisson frequency: {}\n", config.self_consistent_options_2d.m_common.m_poisson_frequency);
    fmt::print("  PBMC bbox: x=[{:.6e},{:.6e}] y=[{:.6e},{:.6e}] z=[{:.6e},{:.6e}] um\n",
               config.self_consistent_options_2d.m_policy.m_pbmc_region_um.get_x_min(),
               config.self_consistent_options_2d.m_policy.m_pbmc_region_um.get_x_max(),
               config.self_consistent_options_2d.m_policy.m_pbmc_region_um.get_y_min(),
               config.self_consistent_options_2d.m_policy.m_pbmc_region_um.get_y_max(),
               config.self_consistent_options_2d.m_policy.m_pbmc_region_um.get_z_min(),
               config.self_consistent_options_2d.m_policy.m_pbmc_region_um.get_z_max());
    fmt::print("  PBMC bbox buffer width: {:.6e} um\n",
               config.self_consistent_options_2d.m_policy.m_buffer_width_um);

    const auto&               common_options = config.self_consistent_options_2d.m_common;
    const auto&               policy_box     = config.self_consistent_options_2d.m_policy.m_pbmc_region_um;
    PBMC::simulation_manifest manifest;
    manifest.add("run", "simulation_type", "self_consistent_device_MMMC");
    manifest.add("run", "simulation_name", config.simulation_name);
    manifest.add("run", "status", "completed");
    manifest.add("run", "started_at_utc", PBMC::current_utc_timestamp());
    manifest.add("run", "command_line", config.command_line);
    manifest.add("run", "working_directory", std::filesystem::current_path().string());
    manifest.add("build", "project_version", PBMC::pbmc_project_version());
    manifest.add("build", "build_type", PBMC::pbmc_build_type());
    manifest.add("build", "compiler", PBMC::pbmc_compiler());
    manifest.add("build", "hardware_concurrency", static_cast<std::size_t>(std::thread::hardware_concurrency()));
    manifest.add("input", "device_mesh", std::filesystem::absolute(config.mesh_file).string());
    manifest.add("input", "material_root", material_repository.root().string());
    manifest.add("input", "material_file", material_repository.material_file(config.material_symbol).string());
    manifest.add("input", "material", config.material_symbol);
    manifest.add("input", "output_directory", std::filesystem::absolute(output_dir).string());
    manifest.add("input", "mesh_dimension", static_cast<int>(mesh->get_dimension()));
    manifest.add("input", "mesh_vertices", static_cast<std::size_t>(mesh->get_nb_vertices()));
    manifest.add("input", "mesh_regions", static_cast<std::size_t>(mesh->get_nb_regions()));
    manifest.add("simulation", "requested_threads", device_options.m_pbmc.m_nb_threads);
    manifest.add("simulation", "random_seed", config.seed_random_generator);
    manifest.add("simulation", "final_time_s", device_options.m_pbmc.m_t_max);
    manifest.add("simulation", "time_step_s", device_options.m_pbmc.m_time_step);
    manifest.add("simulation", "lattice_temperature_K", device_options.m_pbmc.m_lattice_temperature);
    manifest.add("simulation", "max_particles", device_options.m_pbmc.m_max_number_particle);
    manifest.add("simulation", "effective_depth_um", config.self_consistent_options_2d.m_effective_depth_um);
    manifest.add("simulation", "contact_current_window_s", device_options.m_pbmc.m_contact_current_window_s);
    manifest.add("mmmc_policy", "type", "pbmc_bbox_um");
    manifest.add("mmmc_policy", "x_min_um", policy_box.get_x_min());
    manifest.add("mmmc_policy", "x_max_um", policy_box.get_x_max());
    manifest.add("mmmc_policy", "y_min_um", policy_box.get_y_min());
    manifest.add("mmmc_policy", "y_max_um", policy_box.get_y_max());
    manifest.add("mmmc_policy", "z_min_um", policy_box.get_z_min());
    manifest.add("mmmc_policy", "z_max_um", policy_box.get_z_max());
    manifest.add("mmmc_policy",
                 "bbox_buffer_width_um",
                 config.self_consistent_options_2d.m_policy.m_buffer_width_um);
    manifest.add("self_consistent", "poisson_frequency", common_options.m_poisson_frequency);
    manifest.add("self_consistent",
                 "poisson_mode",
                 common_options.m_nonlinear_steady_state_poisson ? "nonlinear_steady_state" : "linear_transient");
    manifest.add("self_consistent", "poisson_mixing_enabled", common_options.m_enable_poisson_mixing);
    manifest.add("self_consistent", "ramo_electrode", common_options.m_ramo_electrode);
    manifest.add("self_consistent",
                 "contact_injection_distribution",
                 PBMC::contact_injection_distribution_name(common_options.m_contact_injection_distribution));
    manifest.add("quench_circuit", "enabled", common_options.m_passive_quench_circuit.m_enabled);
    manifest.add("current_probe", "enabled", device_options.m_pbmc.m_current_probe.m_enabled);
    for (const auto& [contact_name, voltage_V] : common_options.m_contact_voltages_V) {
        manifest.add("contact_voltages_V", contact_name, voltage_V);
    }
    for (const auto& contact_name : config.collecting_contacts) {
        manifest.add("collecting_contacts", contact_name, true);
    }

    const auto                                start = std::chrono::high_resolution_clock::now();
    self_consistent_device_mmmc_simulation_2d simulation(simulation_device,
                                                         device_options,
                                                         config.self_consistent_options_2d,
                                                         material_database,
                                                         config.starting_position,
                                                         config.number_electrons_start,
                                                         config.number_holes_start,
                                                         config.seed_random_generator);
    simulation.set_prefix_export_trajectory_filename(trajectory_dir);
    simulation.run_self_consistent_transport_simulation();
    const auto                          stop    = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> elapsed = stop - start;

    const std::string final_particle_state_file =
        (std::filesystem::path(output_dir) / "final_particle_state.csv").string();
    simulation.export_mmmc_particle_state_csv(final_particle_state_file);
    fmt::print("Wrote {}\n", final_particle_state_file);

    const std::string history_file = (std::filesystem::path(output_dir) / "device_history.csv").string();
    simulation.export_history_to_csv(history_file);
    fmt::print("Wrote {}\n", history_file);
    if (device_options.m_pbmc.m_keep_particles_history) {
        simulation.export_all_mmmc_trajectories_as_csv(trajectory_dir);
        fmt::print("Wrote MMMC particle trajectories to {}\n", trajectory_dir);
    }
    fmt::print("Remaining electrons: {}\n", simulation.get_total_number_electrons());
    fmt::print("Remaining holes: {}\n", simulation.get_total_number_holes());
    fmt::print("Final PBMC particles: {}\n", simulation.get_number_pbmc_particles());
    fmt::print("Final ADMC particles: {}\n", simulation.get_number_admc_particles());
    fmt::print("Total transfers PBMC->ADMC: {}\n", simulation.total_pbmc_to_admc_transfers());
    fmt::print("Total transfers ADMC->PBMC: {}\n", simulation.total_admc_to_pbmc_transfers());
    fmt::print("Simulation completed in {:.3f} s\n", elapsed.count());

    const auto& history = simulation.get_simulation_history();
    manifest.add("run", "finished_at_utc", PBMC::current_utc_timestamp());
    manifest.add("run", "elapsed_seconds", elapsed.count());
    manifest.add("results", "final_time_s", simulation.get_current_time().value_or(0.0));
    manifest.add("results", "remaining_electrons", simulation.get_total_number_electrons());
    manifest.add("results", "remaining_holes", simulation.get_total_number_holes());
    manifest.add("results", "total_electron_weight", simulation.get_total_electron_weight());
    manifest.add("results", "total_hole_weight", simulation.get_total_hole_weight());
    manifest.add("results", "final_pbmc_particles", simulation.get_number_pbmc_particles());
    manifest.add("results", "final_admc_particles", simulation.get_number_admc_particles());
    manifest.add("results", "total_pbmc_to_admc_transfers", simulation.total_pbmc_to_admc_transfers());
    manifest.add("results", "total_admc_to_pbmc_transfers", simulation.total_admc_to_pbmc_transfers());
    manifest.add("results", "impact_ionization_events", history.m_impact_ionization_positions.size());
    if (!history.m_list_ramo_current.empty()) {
        manifest.add("results", "final_ramo_current_A", history.m_list_ramo_current.back());
    }
    manifest.add("results", "avalanche_detected", simulation.avalanche_detection().m_detected);
    manifest.add("results", "successful_quench_detected", simulation.successful_quench_detection().m_detected);
    manifest.add("outputs", "device_history_csv", history_file);
    manifest.add("outputs", "final_particle_state_csv", final_particle_state_file);
    manifest.add("outputs", "trajectory_directory", trajectory_dir);
    manifest.add("outputs", "particle_trajectories_exported", device_options.m_pbmc.m_keep_particles_history);
    manifest.add("outputs", "time_steps_exported", device_options.m_pbmc.m_export_time_step);
    const std::filesystem::path manifest_file = std::filesystem::path(output_dir) / "simulation_manifest.txt";
    manifest.add("outputs", "simulation_manifest", manifest_file.string());
    manifest.write(manifest_file);
    fmt::print("Wrote {}\n", manifest_file.string());
}

}  // namespace uepm::MMMC
