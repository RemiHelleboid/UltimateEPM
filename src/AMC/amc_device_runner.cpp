/**
 * @file amc_device_runner.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-09
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "amc_device_runner.hpp"

#include <fmt/chrono.h>
#include <fmt/core.h>
#include <fmt/format.h>

#include <chrono>
#include <filesystem>
#include <stdexcept>

#include "amc_device_setup.hpp"
#include "device.hpp"
#include "materials.hpp"
#include "msh_file.hpp"

namespace uepm::amc {

void run_self_consistent_device_amc_simulation(const self_consistent_device_amc_run_config& config) {
    validate_material_symbol(config.material_symbol);
    auto device_options = config.device_options;
    device_options.m_simulation_name  = config.simulation_name;
    device_options.m_output_directory = config.output_dir;

    auto self_consistent_options_2d = config.self_consistent_options_2d;
    auto self_consistent_options_3d = config.self_consistent_options_3d;

    device_options.validate();
    self_consistent_options_2d.validate();

    const std::string output_dir = config.output_dir.empty() ? make_default_output_directory(config.mesh_file)
                                                             : config.output_dir;
    const std::string trajectory_dir = fmt::format("{}/trajectory", output_dir);

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
    fmt::print("Loading materials: {}\n", config.material_file);

    uepm::physic::material::list_materials list_of_materials;
    list_of_materials.load_materials_from_file(config.material_file);

    uepm::device::device simulation_device(mesh);
    add_default_contacts(simulation_device, *mesh);

    const std::string material_symbol = config.material_symbol;

    fmt::print("Self-consistent AMC {}D simulation\n", mesh_dimension);
    fmt::print("  mesh vertices: {}\n", mesh->get_nb_vertices());
    fmt::print("  output directory: {}\n", output_dir);
    fmt::print("  material: {}\n", material_symbol);
    fmt::print("  final time: {:.6e} s\n", config.device_options.m_t_max);
    fmt::print("  time step: {:.6e} s\n", config.device_options.m_time_step);
    fmt::print("  Poisson frequency: {}\n", config.self_consistent_options_2d.m_common.m_poisson_frequency);
    fmt::print("  anode voltage: {:.6e} V\n", config.self_consistent_options_2d.m_common.m_anode_voltage);
    fmt::print("  cathode voltage: {:.6e} V\n", config.self_consistent_options_2d.m_common.m_cathode_voltage);
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

    const auto start = std::chrono::high_resolution_clock::now();

    if (mesh_dimension == 2) {
        self_consistent_device_amc_simulation_2d simulation(simulation_device,
                                    device_options,
                                                            config.self_consistent_options_2d,
                                                            list_of_materials,
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
    } else {
        self_consistent_device_amc_simulation_3d simulation(simulation_device,
                                    device_options,
                                    self_consistent_options_3d,
                                                            list_of_materials,
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
        if (device_options.m_enable_scheduled_particle_injection) {
            const auto& injection = device_options.m_scheduled_particle_injection;
            fmt::print("  scheduled injection: enabled\n");
            fmt::print("    time: {:.6e} s\n", injection.m_time_s);
            fmt::print("    position: ({:.6e}, {:.6e}, {:.6e}) um\n",
                       injection.m_position_um.x(),
                       injection.m_position_um.y(),
                       injection.m_position_um.z());
            fmt::print("    type: {}\n",
                       injection.m_particle_type == particle_type::electron ? "electron" : "hole");
            fmt::print("    weight: {:.6e}\n", injection.m_weight);
        }
    }

    const auto                          stop    = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> elapsed = stop - start;
    fmt::print("Simulation completed in {:.3f} s\n", elapsed.count());
}

}  // namespace uepm::amc