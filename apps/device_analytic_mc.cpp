/**
 * @file device_analytic_mc.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-05-26
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#include <fmt/chrono.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <tclap/CmdLine.h>

#include <chrono>
#include <filesystem>
#include <stdexcept>
#include <string>

#include "device.hpp"
#include "device_amc_simulation.hpp"
#include "mesh.hpp"
#include "msh_file.hpp"
#include "vector.hpp"

namespace {

uepm::mesh::vector3 make_vector3(double x, double y, double z) { return uepm::mesh::vector3{x, y, z}; }

void validate_material(const std::string& material_symbol) {
    if (material_symbol != "Si") {
        throw std::invalid_argument("Only Si is currently supported by the analytical AMC transport model.");
    }
}

void validate_options(const uepm::amc::options_device_amc& options) {
    if (options.m_time_step <= 0.0) {
        throw std::invalid_argument("--dt must be positive.");
    }

    if (options.m_t_max <= 0.0) {
        throw std::invalid_argument("--time must be positive.");
    }

    if (options.m_lattice_temperature < 0.0) {
        throw std::invalid_argument("--temperature must be non-negative.");
    }

    if (options.m_max_energy_eV <= 0.0) {
        throw std::invalid_argument("--max-energy must be positive.");
    }

    if (options.m_self_scattering_safety_factor <= 0.0) {
        throw std::invalid_argument("--gamma-safety must be positive.");
    }

    if (options.m_gamma_max_energy_samples < 2) {
        throw std::invalid_argument("--gamma-samples must be at least 2.");
    }

    if (options.m_frequency_export_trajectory <= 0) {
        throw std::invalid_argument("--export-frequency must be positive.");
    }

    if (options.m_nb_threads <= 0) {
        throw std::invalid_argument("--nthreads must be positive.");
    }

    if (options.m_avalanche_threshold == 0) {
        throw std::invalid_argument("--avalanche-threshold must be positive.");
    }

    if (options.m_max_number_particle == 0) {
        throw std::invalid_argument("--max-particles must be positive.");
    }

    if (options.m_avalanche_threshold > options.m_max_number_particle) {
        throw std::invalid_argument("--avalanche-threshold cannot be larger than --max-particles.");
    }
}

std::string make_default_output_directory(const std::string& mesh_filename) {
    const std::string mesh_stem = std::filesystem::path(mesh_filename).stem().string();

    return fmt::format("device_amc_{}", mesh_stem);
}

}  // namespace

int main(int argc, const char** argv) {
    try {
        TCLAP::CmdLine cmd("Analytical Monte Carlo device simulation.", ' ', "1.0");

        TCLAP::ValueArg<std::string> arg_device_mesh("", "device-mesh", "Path to the device mesh .msh file.", true, "", "path");

        TCLAP::ValueArg<std::string> arg_material("m",
                                                  "material",
                                                  "Material symbol. Currently only Si is supported.",
                                                  false,
                                                  "Si",
                                                  "string");

        TCLAP::ValueArg<std::string> arg_output_dir("d", "outdir", "Output directory.", false, "", "path");
        TCLAP::ValueArg<std::string> arg_simulation_name("", "name", "Simulation name.", false, "device_amc", "string");
        TCLAP::ValueArg<double> arg_time("t", "time", "Final simulation time in seconds.", false, 1.0e-12, "seconds");
        TCLAP::ValueArg<double> arg_dt("", "dt", "Global device time step in seconds.", false, 1.0e-15, "seconds");
        TCLAP::ValueArg<double> arg_temperature("T", "temperature", "Lattice temperature in K.", false, 300.0, "K");
        TCLAP::ValueArg<double> arg_max_energy("e",
                                               "max-energy",
                                               "Maximum energy in eV used for gamma_max computation.",
                                               false,
                                               10.0,
                                               "eV");

        TCLAP::ValueArg<double> arg_gamma_safety("", "gamma-safety", "Safety factor applied to gamma_max.", false, 1.2, "double");
        TCLAP::ValueArg<int> arg_gamma_samples("",
                                               "gamma-samples",
                                               "Number of energy samples used for gamma_max computation.",
                                               false,
                                               1000,
                                               "int");
        TCLAP::ValueArg<int> arg_nb_threads("j", "nthreads", "Number of OpenMP threads.", false, 1, "int");
        TCLAP::ValueArg<int> arg_seed("", "seed", "Random seed.", false, 0, "int");
        TCLAP::ValueArg<std::size_t> arg_number_electrons("", "nelectrons", "Initial number of electrons.", false, 1, "integer");
        TCLAP::ValueArg<std::size_t> arg_number_holes("", "nholes", "Initial number of holes.", false, 0, "integer");
        TCLAP::ValueArg<double> arg_start_x("", "x0", "Initial particle x position in meters.", true, 0.0, "m");
        TCLAP::ValueArg<double> arg_start_y("", "y0", "Initial particle y position in meters.", false, 0.0, "m");
        TCLAP::ValueArg<double> arg_start_z("", "z0", "Initial particle z position in meters.", false, 0.0, "m");
        TCLAP::ValueArg<std::size_t> arg_max_particles("",
                                                       "max-particles",
                                                       "Hard maximum number of active particles.",
                                                       false,
                                                       1000000,
                                                       "integer");

        TCLAP::ValueArg<std::size_t> arg_avalanche_threshold("",
                                                             "avalanche-threshold",
                                                             "Particle count threshold used to stop the simulation.",
                                                             false,
                                                             1000,
                                                             "integer");
        TCLAP::SwitchArg arg_disable_impact_ionization("", "disable-impact-ionization", "Disable impact-ionization computation.", false);
        TCLAP::SwitchArg arg_disable_particle_creation("",
                                                       "disable-particle-creation",
                                                       "Compute impact ionization but do not create electron-hole pairs.",
                                                       false);
        TCLAP::SwitchArg arg_keep_history("H", "keep-particle-history", "Store full particle trajectories in memory.", false);
        TCLAP::SwitchArg arg_export_time_steps("E",
                                               "export-time-steps",
                                               "Export particle state periodically during the simulation.",
                                               false);

        TCLAP::ValueArg<int> arg_export_frequency("", "export-frequency", "Export one time step every N iterations.", false, 10, "integer");
        TCLAP::SwitchArg arg_do_not_stop_without_electrons("",
                                                           "keep-going-without-electrons",
                                                           "Do not stop when no electrons remain in the device.",
                                                           false);

        cmd.add(arg_device_mesh);
        cmd.add(arg_material);
        cmd.add(arg_output_dir);
        cmd.add(arg_simulation_name);
        cmd.add(arg_time);
        cmd.add(arg_dt);
        cmd.add(arg_temperature);
        cmd.add(arg_max_energy);
        cmd.add(arg_gamma_safety);
        cmd.add(arg_gamma_samples);
        cmd.add(arg_nb_threads);
        cmd.add(arg_seed);
        cmd.add(arg_number_electrons);
        cmd.add(arg_number_holes);
        cmd.add(arg_start_x);
        cmd.add(arg_start_y);
        cmd.add(arg_start_z);
        cmd.add(arg_max_particles);
        cmd.add(arg_avalanche_threshold);
        cmd.add(arg_disable_impact_ionization);
        cmd.add(arg_disable_particle_creation);
        cmd.add(arg_keep_history);
        cmd.add(arg_export_time_steps);
        cmd.add(arg_export_frequency);
        cmd.add(arg_do_not_stop_without_electrons);

        cmd.parse(argc, argv);

        const std::string device_mesh_file = arg_device_mesh.getValue();
        const std::string material_symbol  = arg_material.getValue();

        validate_material(material_symbol);

        uepm::amc::options_device_amc options;
        options.m_t_max                                = arg_time.getValue();
        options.m_time_step                            = arg_dt.getValue();
        options.m_lattice_temperature                  = arg_temperature.getValue();
        options.m_max_energy_eV                        = arg_max_energy.getValue();
        options.m_self_scattering_safety_factor        = arg_gamma_safety.getValue();
        options.m_gamma_max_energy_samples             = static_cast<std::size_t>(arg_gamma_samples.getValue());
        options.m_nb_threads                           = arg_nb_threads.getValue();
        options.m_max_number_particle                  = arg_max_particles.getValue();
        options.m_avalanche_threshold                  = arg_avalanche_threshold.getValue();
        options.m_activate_impact_ionization           = !arg_disable_impact_ionization.getValue();
        options.m_particle_creation_activated          = !arg_disable_particle_creation.getValue();
        options.m_stop_simu_when_no_electron_remaining = !arg_do_not_stop_without_electrons.getValue();
        options.m_keep_particles_history               = arg_keep_history.getValue();
        options.m_export_time_step                     = arg_export_time_steps.getValue();
        options.m_frequency_export_trajectory          = arg_export_frequency.getValue();

        validate_options(options);

        const std::string output_dir =
            arg_output_dir.getValue().empty() ? make_default_output_directory(device_mesh_file) : arg_output_dir.getValue();

        std::filesystem::create_directories(output_dir);

        const std::string trajectory_dir = fmt::format("{}/trajectory", output_dir);

        if (options.m_export_time_step || options.m_keep_particles_history) {
            std::filesystem::create_directories(trajectory_dir);
        }

        uepm::file::msh_file mesh_file(device_mesh_file);
        mesh_file.read_mesh();
        mesh_file.read_states();

        uepm::mesh::mesh* mesh = mesh_file.get_p_mesh();

        if (mesh == nullptr) {
            throw std::runtime_error("Failed to load mesh.");
        }

        uepm::device::device simulation_device(mesh);

        const auto starting_position = make_vector3(arg_start_x.getValue(), arg_start_y.getValue(), arg_start_z.getValue());

        fmt::print("Device AMC simulation\n");
        fmt::print("  mesh: {}\n", device_mesh_file);
        fmt::print("  material: {}\n", material_symbol);
        fmt::print("  output directory: {}\n", output_dir);
        fmt::print("  mesh vertices: {}\n", mesh->get_nb_vertices());
        fmt::print("  final time: {:.6e} s\n", options.m_t_max);
        fmt::print("  time step: {:.6e} s\n", options.m_time_step);
        fmt::print("  temperature: {:.3f} K\n", options.m_lattice_temperature);
        fmt::print("  initial electrons: {}\n", arg_number_electrons.getValue());
        fmt::print("  initial holes: {}\n", arg_number_holes.getValue());
        fmt::print("  starting position: ({:.6e}, {:.6e}, {:.6e}) m\n",
                   starting_position.x(),
                   starting_position.y(),
                   starting_position.z());

        uepm::amc::device_amc_simulation simulation(simulation_device,
                                                    options,
                                                    arg_simulation_name.getValue(),
                                                    starting_position,
                                                    arg_number_electrons.getValue(),
                                                    arg_number_holes.getValue(),
                                                    arg_seed.getValue());

        simulation.set_prefix_export_trajectory_filename(fmt::format("{}/time_step", trajectory_dir));

        const auto start = std::chrono::high_resolution_clock::now();

        simulation.run();

        const auto                          end     = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> elapsed = end - start;

        fmt::print("Simulation completed in {:.3f} seconds.\n", elapsed.count());
        fmt::print("Final simulation time: {:.6e} s\n", simulation.get_current_time().value_or(0.0));
        fmt::print("Remaining electrons: {}\n", simulation.get_number_electrons());
        fmt::print("Remaining holes: {}\n", simulation.get_number_holes());

        const std::string history_file = fmt::format("{}/device_history.csv", output_dir);

        simulation.export_history_to_csv(history_file);

        fmt::print("Wrote {}\n", history_file);

        if (options.m_keep_particles_history) {
            const std::string trajectory_prefix = fmt::format("{}/particle_", trajectory_dir);

            simulation.export_all_trajectories_as_csv(trajectory_prefix);

            fmt::print("Wrote particle trajectories to {}\n", trajectory_dir);
        }

        return 0;
    } catch (const TCLAP::ArgException& error) {
        fmt::print(stderr, "Argument error: {} for argument {}\n", error.error(), error.argId());
        return 1;
    } catch (const std::exception& error) {
        fmt::print(stderr, "Error: {}\n", error.what());
        return 1;
    }
}