/**
 * @file admc_device_runner.cpp
 * @brief Runner support for self-consistent ADMC device simulations.
 */

#include "admc_device_runner.hpp"

#include <fmt/core.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>

#include "msh_file.hpp"

namespace uepm::ADMC {
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

std::string current_utc_timestamp() {
    const auto now    = std::chrono::system_clock::now();
    const auto time_t = std::chrono::system_clock::to_time_t(now);
    std::tm    utc_time{};
#if defined(_WIN32)
    gmtime_s(&utc_time, &time_t);
#else
    gmtime_r(&time_t, &utc_time);
#endif
    std::ostringstream stream;
    stream << std::put_time(&utc_time, "%Y-%m-%dT%H:%M:%SZ");
    return stream.str();
}

std::string make_default_output_directory(const std::string& mesh_file) {
    return fmt::format("self_consistent_admc_{}", std::filesystem::path(mesh_file).stem().string());
}

void validate_material_symbol(const std::string& material_symbol) {
    if (material_symbol != "Si") {
        throw std::invalid_argument("Only Si is currently supported by the ADMC mobility model.");
    }
}

void add_collecting_contacts(uepm::device::device&           simulation_device,
                             uepm::mesh::mesh&               mesh,
                             const std::vector<std::string>& contact_names) {
    constexpr double contact_collection_depth = 0.001;
    constexpr double contact_margin           = 10.0;
    constexpr double ohmic_resistance         = 0.0;

    const mesh::bbox device_bbox   = mesh.get_bounding_box();
    const double     tolerance     = 1.0e-8 * std::max(device_bbox.get_diagonal_size(), 1.0);
    const auto       bulk_elements = mesh.get_list_bulk_element();

    for (const auto& contact_name : contact_names) {
        const auto* region = mesh.get_p_region(contact_name);
        if (region == nullptr) {
            throw std::runtime_error("Collecting contact region '" + contact_name + "' does not exist in the mesh.");
        }
        if (region->get_region_type() != mesh::RegionType::contact) {
            throw std::runtime_error("Collecting contact '" + contact_name + "' is not a mesh contact region.");
        }

        const mesh::bbox            region_bbox = region->compute_bounding_box();
        const std::array<double, 3> region_min{region_bbox.get_x_min(),
                                               region_bbox.get_y_min(),
                                               region_bbox.get_z_min()};
        const std::array<double, 3> region_max{region_bbox.get_x_max(),
                                               region_bbox.get_y_max(),
                                               region_bbox.get_z_max()};
        std::array<double, 3>       box_min{};
        std::array<double, 3>       box_max{};
        std::optional<std::size_t>  normal_axis;

        for (std::size_t axis = 0; axis < 3; ++axis) {
            const double region_size = region_max[axis] - region_min[axis];
            const double device_size = axis == 0   ? device_bbox.get_x_size()
                                       : axis == 1 ? device_bbox.get_y_size()
                                                   : device_bbox.get_z_size();
            if (device_size > tolerance && region_size <= tolerance) {
                if (normal_axis.has_value()) {
                    throw std::runtime_error("Cannot determine a unique normal for collecting contact '" +
                                             contact_name + "'.");
                }
                normal_axis = axis;
            } else if (device_size <= tolerance) {
                box_min[axis] = region_min[axis] - contact_margin;
                box_max[axis] = region_max[axis] + contact_margin;
            } else {
                box_min[axis] = region_min[axis] - tolerance;
                box_max[axis] = region_max[axis] + tolerance;
            }
        }

        if (!normal_axis.has_value()) {
            throw std::runtime_error("Cannot determine the normal of collecting contact '" + contact_name + "'.");
        }

        const auto adjacent_indices = mesh.get_idx_bulk_elements_adjacent_to_contact_region(contact_name);
        if (adjacent_indices.empty()) {
            throw std::runtime_error("Collecting contact '" + contact_name + "' has no adjacent bulk element.");
        }

        double mean_adjacent_coordinate = 0.0;
        for (const auto element_index : adjacent_indices) {
            if (element_index >= bulk_elements.size()) {
                throw std::runtime_error("Contact-adjacent bulk element index is out of range.");
            }
            const auto barycenter = bulk_elements[element_index]->get_barycenter();
            mean_adjacent_coordinate += (*normal_axis == 0   ? barycenter.x()
                                         : *normal_axis == 1 ? barycenter.y()
                                                             : barycenter.z());
        }
        mean_adjacent_coordinate /= static_cast<double>(adjacent_indices.size());

        const double contact_plane = 0.5 * (region_min[*normal_axis] + region_max[*normal_axis]);
        if (mean_adjacent_coordinate < contact_plane) {
            box_min[*normal_axis] = contact_plane - contact_collection_depth;
            box_max[*normal_axis] = contact_plane + contact_margin;
        } else {
            box_min[*normal_axis] = contact_plane - contact_margin;
            box_max[*normal_axis] = contact_plane + contact_collection_depth;
        }

        simulation_device.add_contact(contact_name,
                                      mesh::vector3{box_min[0], box_min[1], box_min[2]},
                                      mesh::vector3{box_max[0], box_max[1], box_max[2]},
                                      ohmic_resistance);
    }
}

void write_manifest(const std::filesystem::path&                     filename,
                    const self_consistent_device_admc_run_config&    config,
                    const std::string&                               output_dir,
                    const uepm::physics::material_repository&        material_repository,
                    int                                              mesh_dimension,
                    std::size_t                                      mesh_vertices,
                    std::size_t                                      mesh_regions,
                    double                                           elapsed_seconds,
                    const self_consistent_device_admc_simulation_2d& simulation) {
    std::ofstream stream(filename);
    if (!stream.is_open()) {
        throw std::runtime_error("Could not open ADMC manifest for writing: " + filename.string());
    }
    stream << "# UltimateEPM ADMC simulation manifest\n";
    stream << "\n[run]\n";
    stream << "simulation_type = self_consistent_device_ADMC\n";
    stream << "simulation_name = " << config.simulation_name << "\n";
    stream << "finished_at_utc = " << current_utc_timestamp() << "\n";
    stream << "command_line = " << config.command_line << "\n";
    stream << "working_directory = " << std::filesystem::current_path().string() << "\n";
    stream << "elapsed_seconds = " << elapsed_seconds << "\n";
    stream << "\n[input]\n";
    stream << "device_mesh = " << std::filesystem::absolute(config.mesh_file).string() << "\n";
    stream << "material_root = " << material_repository.root().string() << "\n";
    stream << "material_file = " << material_repository.material_file(config.material_symbol).string() << "\n";
    stream << "material = " << config.material_symbol << "\n";
    stream << "output_directory = " << std::filesystem::absolute(output_dir).string() << "\n";
    stream << "mesh_dimension = " << mesh_dimension << "\n";
    stream << "mesh_vertices = " << mesh_vertices << "\n";
    stream << "mesh_regions = " << mesh_regions << "\n";
    stream << "\n[simulation]\n";
    stream << "random_seed = " << config.random_seed << "\n";
    stream << "final_time_s = " << config.device_options.m_final_time_s << "\n";
    stream << "time_step_s = " << config.device_options.m_time_step_s << "\n";
    stream << "lattice_temperature_K = " << config.device_options.m_lattice_temperature_K << "\n";
    stream << "initial_electrons = " << config.number_electrons_start << "\n";
    stream << "initial_holes = " << config.number_holes_start << "\n";
    stream << "initial_x_um = " << config.starting_position_um.x() << "\n";
    stream << "initial_y_um = " << config.starting_position_um.y() << "\n";
    stream << "initial_z_um = " << config.starting_position_um.z() << "\n";
    stream << "max_particles = " << config.device_options.m_max_number_particles << "\n";
    stream << "stop_when_no_electrons = " << (config.device_options.m_stop_when_no_electrons ? "true" : "false")
           << "\n";
    stream << "export_time_steps = " << (config.device_options.m_export_time_step ? "true" : "false") << "\n";
    stream << "export_frequency = " << config.device_options.m_frequency_export << "\n";
    stream << "mesh_particle_local_averages = "
           << (config.device_options.m_export_mesh_particle_local_averages ? "true" : "false") << "\n";
    stream << "boundary_reflection = "
           << mesh::boundary_reflection_model_name(config.device_options.m_boundary_reflection_model) << "\n";
    stream << "current_probe_enabled = " << (config.device_options.m_current_probe.m_enabled ? "true" : "false")
           << "\n";
    stream << "current_probe_x_min_um = " << config.device_options.m_current_probe.m_box_um.get_x_min() << "\n";
    stream << "current_probe_x_max_um = " << config.device_options.m_current_probe.m_box_um.get_x_max() << "\n";
    stream << "current_probe_y_min_um = " << config.device_options.m_current_probe.m_box_um.get_y_min() << "\n";
    stream << "current_probe_y_max_um = " << config.device_options.m_current_probe.m_box_um.get_y_max() << "\n";
    stream << "current_probe_z_min_um = " << config.device_options.m_current_probe.m_box_um.get_z_min() << "\n";
    stream << "current_probe_z_max_um = " << config.device_options.m_current_probe.m_box_um.get_z_max() << "\n";
    stream << "effective_depth_um = " << config.self_consistent_options_2d.m_effective_depth_um << "\n";
    stream << "particle_z_period_um = " << config.self_consistent_options_2d.m_particle_z_period_um << "\n";
    stream << "poisson_frequency = " << config.self_consistent_options_2d.m_common.m_poisson_frequency << "\n";
    stream << "\n[self_consistent]\n";
    stream << "ramo_electrode = " << config.self_consistent_options_2d.m_common.m_ramo_electrode << "\n";
    stream << "built_in_potential_enabled = "
           << (config.self_consistent_options_2d.m_common.m_enable_built_in_potential ? "true" : "false") << "\n";
    stream << "intrinsic_concentration_model = silicon_varshni_normalized_1e10_cm-3_at_300K\n";
    stream << "intrinsic_concentration_cm_3 = "
           << silicon_intrinsic_concentration_admc_cm_3(config.device_options.m_lattice_temperature_K) << "\n";
    stream << "built_in_contact_voltage_scale = "
           << config.self_consistent_options_2d.m_common.m_built_in_contact_voltage_scale << "\n";
    stream << "poisson_mixing_enabled = "
           << (config.self_consistent_options_2d.m_common.m_enable_poisson_mixing ? "true" : "false") << "\n";
    stream << "poisson_mixing_old_solution_fraction = "
           << config.self_consistent_options_2d.m_common.m_poisson_mixing_old_solution_fraction << "\n";
    stream << "initialize_particles_from_doping = "
           << (config.self_consistent_options_2d.m_common.m_initialize_particles_from_doping ? "true" : "false")
           << "\n";
    stream << "initial_particle_weight = " << config.self_consistent_options_2d.m_common.m_initial_particle_weight
           << "\n";
    stream << "initial_particle_state_file = "
           << config.self_consistent_options_2d.m_common.m_initial_particle_state_file << "\n";
    stream << "contact_injection_particle_weight = "
           << config.self_consistent_options_2d.m_common.m_contact_injection_particle_weight << "\n";
    stream << "\n[contact_voltages_V]\n";
    for (const auto& [contact_name, voltage_V] : config.self_consistent_options_2d.m_common.m_contact_voltages_V) {
        stream << contact_name << " = " << voltage_V << "\n";
    }
    stream << "\n[collecting_contacts]\n";
    for (const auto& contact_name : config.collecting_contacts) {
        stream << contact_name << " = true\n";
    }
    stream << "\n[scheduled_injection]\n";
    stream << "enabled = " << (config.device_options.m_enable_scheduled_particle_injection ? "true" : "false") << "\n";
    stream << "time_s = " << config.device_options.m_scheduled_injection_time_s << "\n";
    stream << "x_um = " << config.device_options.m_scheduled_injection_position_um.x() << "\n";
    stream << "y_um = " << config.device_options.m_scheduled_injection_position_um.y() << "\n";
    stream << "z_um = " << config.device_options.m_scheduled_injection_position_um.z() << "\n";
    stream << "particle_type = "
           << (config.device_options.m_scheduled_injection_type == carrier_type::electron ? "electron" : "hole")
           << "\n";
    stream << "weight = " << config.device_options.m_scheduled_injection_weight << "\n";
    stream << "\n[results]\n";
    stream << "final_time_s = " << simulation.current_time_s() << "\n";
    stream << "remaining_electrons = " << simulation.get_number_electrons() << "\n";
    stream << "remaining_holes = " << simulation.get_number_holes() << "\n";
    stream << "total_electron_weight = " << simulation.get_total_electron_weight() << "\n";
    stream << "total_hole_weight = " << simulation.get_total_hole_weight() << "\n";
    stream << "\n[outputs]\n";
    stream << "device_history_csv = " << (std::filesystem::path(output_dir) / "device_history.csv").string() << "\n";
    stream << "trajectory_directory = " << (std::filesystem::path(output_dir) / "trajectory").string() << "\n";
    stream << "particle_trajectories_exported = false\n";
    stream << "time_steps_exported = " << (config.device_options.m_export_time_step ? "true" : "false") << "\n";
    stream << "simulation_manifest = " << filename.string() << "\n";
}

}  // namespace

std::string admc_command_line_from_arguments(int argc, const char* const* argv) {
    std::string command_line;
    for (int index = 0; index < argc; ++index) {
        if (!command_line.empty()) {
            command_line.push_back(' ');
        }
        command_line += quote_command_argument(argv[index]);
    }
    return command_line;
}

void run_self_consistent_device_admc_simulation(const self_consistent_device_admc_run_config& config) {
    validate_material_symbol(config.material_symbol);
    config.device_options.validate();
    config.self_consistent_options_2d.validate();

    const std::string output_dir =
        config.output_dir.empty() ? make_default_output_directory(config.mesh_file) : config.output_dir;
    std::filesystem::create_directories(output_dir);
    auto device_options = config.device_options;
    device_options.m_output_directory = output_dir;

    fmt::print("Loading mesh: {}\n", config.mesh_file);
    uepm::file::msh_file msh_file(config.mesh_file);
    msh_file.read_mesh();
    msh_file.read_states();
    auto* mesh = msh_file.get_p_mesh();
    if (mesh == nullptr) {
        throw std::runtime_error("Mesh loading failed.");
    }
    if (mesh->get_dimension() != 2) {
        throw std::runtime_error("ADMC self-consistent device mode currently supports 2D meshes.");
    }

    const uepm::physics::material_repository material_repository =
        config.material_root.empty() ? uepm::physics::material_repository{}
                                     : uepm::physics::material_repository{config.material_root};
    auto material_database = material_repository.load_all_materials();
    (void)material_database.require(config.material_symbol);

    uepm::device::device simulation_device(mesh);
    add_collecting_contacts(simulation_device, *mesh, config.collecting_contacts);

    fmt::print("Self-consistent ADMC 2D simulation\n");
    fmt::print("  output directory: {}\n", output_dir);
    fmt::print("  material: {}\n", config.material_symbol);
    fmt::print("  final time: {:.6e} s\n", device_options.m_final_time_s);
    fmt::print("  time step: {:.6e} s\n", device_options.m_time_step_s);
    fmt::print("  Poisson frequency: {}\n", config.self_consistent_options_2d.m_common.m_poisson_frequency);
    fmt::print("  contact voltages:\n");
    for (const auto& [contact_name, voltage_V] : config.self_consistent_options_2d.m_common.m_contact_voltages_V) {
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
    fmt::print("  Poisson mixing: {}\n",
               config.self_consistent_options_2d.m_common.m_enable_poisson_mixing ? "enabled" : "disabled");
    if (config.self_consistent_options_2d.m_common.m_enable_poisson_mixing) {
        fmt::print("    old solution fraction: {:.6e}\n",
                   config.self_consistent_options_2d.m_common.m_poisson_mixing_old_solution_fraction);
    }
    if (config.self_consistent_options_2d.m_common.m_enable_built_in_potential) {
        fmt::print("  intrinsic concentration at {:.3f} K: {:.6e} cm^-3\n",
                   device_options.m_lattice_temperature_K,
                   silicon_intrinsic_concentration_admc_cm_3(device_options.m_lattice_temperature_K));
    }
    fmt::print("  export time steps: {}\n", device_options.m_export_time_step ? "enabled" : "disabled");
    fmt::print("  current probe: {}\n", device_options.m_current_probe.m_enabled ? "enabled" : "disabled");
    if (device_options.m_current_probe.m_enabled) {
        const auto& box = device_options.m_current_probe.m_box_um;
        fmt::print("    x: [{:.6e}, {:.6e}] um\n", box.get_x_min(), box.get_x_max());
        fmt::print("    y: [{:.6e}, {:.6e}] um\n", box.get_y_min(), box.get_y_max());
        fmt::print("    z: [{:.6e}, {:.6e}] um\n", box.get_z_min(), box.get_z_max());
    }
    fmt::print("  effective depth: {:.6e} um\n", config.self_consistent_options_2d.m_effective_depth_um);
    fmt::print("  particle z period: {:.6e} um\n", config.self_consistent_options_2d.m_particle_z_period_um);
    fmt::print("  initial electrons: {}\n", config.number_electrons_start);
    fmt::print("  initial holes: {}\n", config.number_holes_start);
    fmt::print("  initial position: ({:.6e}, {:.6e}, {:.6e})\n",
               config.starting_position_um.x(),
               config.starting_position_um.y(),
               config.starting_position_um.z());
    if (device_options.m_enable_scheduled_particle_injection) {
        fmt::print("  scheduled injection: enabled\n");
        fmt::print("    time: {:.6e} s\n", device_options.m_scheduled_injection_time_s);
        fmt::print("    position: ({:.6e}, {:.6e}, {:.6e}) um\n",
                   device_options.m_scheduled_injection_position_um.x(),
                   device_options.m_scheduled_injection_position_um.y(),
                   device_options.m_scheduled_injection_position_um.z());
        fmt::print("    type: {}\n",
                   device_options.m_scheduled_injection_type == carrier_type::electron ? "electron" : "hole");
        fmt::print("    weight: {:.6e}\n", device_options.m_scheduled_injection_weight);
    }

    const auto                                start = std::chrono::high_resolution_clock::now();
    self_consistent_device_admc_simulation_2d simulation(simulation_device,
                                                         device_options,
                                                         config.self_consistent_options_2d,
                                                         material_database,
                                                         config.starting_position_um,
                                                         config.number_electrons_start,
                                                         config.number_holes_start,
                                                         config.random_seed);
    simulation.set_prefix_export_trajectory_filename((std::filesystem::path(output_dir) / "trajectory").string());
    simulation.run_self_consistent_transport_simulation();
    const auto                          stop    = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> elapsed = stop - start;

    const auto history_file = std::filesystem::path(output_dir) / "device_history.csv";
    simulation.export_history_to_csv(history_file.string());
    fmt::print("Wrote {}\n", history_file.string());
    fmt::print("Remaining electrons: {}\n", simulation.get_number_electrons());
    fmt::print("Remaining holes: {}\n", simulation.get_number_holes());
    fmt::print("Simulation completed in {:.3f} s\n", elapsed.count());

    const auto manifest_file = std::filesystem::path(output_dir) / "simulation_manifest.txt";
    write_manifest(manifest_file,
                   config,
                   output_dir,
                   material_repository,
                   mesh->get_dimension(),
                   static_cast<std::size_t>(mesh->get_nb_vertices()),
                   static_cast<std::size_t>(mesh->get_nb_regions()),
                   elapsed.count(),
                   simulation);
    fmt::print("Wrote {}\n", manifest_file.string());
}

}  // namespace uepm::ADMC
