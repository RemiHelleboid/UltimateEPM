/**
 * @file mmmc_device_runner.cpp
 * @brief Runner support for self-consistent MMMC device simulations.
 */

#include "mmmc_device_runner.hpp"

#include <fmt/core.h>

#include <chrono>
#include <filesystem>
#include <sstream>
#include <stdexcept>

#include "materials.hpp"
#include "msh_file.hpp"
#include "pbmc_device_setup.hpp"
#include "pbmc_material_model.hpp"

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

    const physics::material_repository material_repository =
        config.material_root.empty() ? physics::material_repository{}
                                     : physics::material_repository{config.material_root};
    physics::material_database material_database = material_repository.load_all_materials();
    const auto& common_material = material_database.require(config.material_symbol);

    auto device_options = config.device_options;
    device_options.m_pbmc.m_simulation_name = config.simulation_name;
    device_options.m_pbmc.m_output_directory = output_dir;
    device_options.m_pbmc.m_material_model = PBMC::load_pbmc_material_model(material_repository, common_material);
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

    const auto start = std::chrono::high_resolution_clock::now();
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
    const auto stop = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> elapsed = stop - start;

    const std::string final_particle_state_file =
        (std::filesystem::path(output_dir) / "final_particle_state.csv").string();
    simulation.export_mmmc_particle_state_csv(final_particle_state_file);
    fmt::print("Wrote {}\n", final_particle_state_file);

    const std::string history_file = (std::filesystem::path(output_dir) / "device_history.csv").string();
    simulation.export_history_to_csv(history_file);
    fmt::print("Wrote {}\n", history_file);
    fmt::print("Remaining electrons: {}\n", simulation.get_total_number_electrons());
    fmt::print("Remaining holes: {}\n", simulation.get_total_number_holes());
    fmt::print("Final PBMC particles: {}\n", simulation.get_number_pbmc_particles());
    fmt::print("Final ADMC particles: {}\n", simulation.get_number_admc_particles());
    fmt::print("Total transfers PBMC->ADMC: {}\n", simulation.total_pbmc_to_admc_transfers());
    fmt::print("Total transfers ADMC->PBMC: {}\n", simulation.total_admc_to_pbmc_transfers());
    fmt::print("Simulation completed in {:.3f} s\n", elapsed.count());
}

}  // namespace uepm::MMMC
