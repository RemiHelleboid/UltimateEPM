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

#include "amc_self_consistent_device_simulation_3d.hpp"
#include "device.hpp"
#include "materials.hpp"
#include "mesh.hpp"
#include "msh_file.hpp"
#include "vector.hpp"

namespace {

uepm::mesh::vector3 make_vector3(double x, double y, double z) { return uepm::mesh::vector3{x, y, z}; }

void validate_material_symbol(const std::string& material_symbol) {
    if (material_symbol != "Si") {
        throw std::invalid_argument("Only Si is currently supported by the analytical AMC transport model.");
    }
}

void validate_device_options(const uepm::amc::options_device_amc& options) {
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

    if (options.m_max_number_particle == 0) {
        throw std::invalid_argument("--max-particles must be positive.");
    }

    if (options.m_avalanche_threshold == 0) {
        throw std::invalid_argument("--avalanche-threshold must be positive.");
    }

    if (options.m_avalanche_threshold > options.m_max_number_particle) {
        throw std::invalid_argument("--avalanche-threshold cannot be larger than --max-particles.");
    }

    if (options.m_frequency_export_trajectory <= 0) {
        throw std::invalid_argument("--export-frequency must be positive.");
    }

    if (options.m_nb_threads <= 0) {
        throw std::invalid_argument("--nthreads must be positive.");
    }
}

void validate_self_consistent_options(const uepm::amc::options_self_consistent_device_amc_3d& options) {
    if (options.m_poisson_frequency == 0) {
        throw std::invalid_argument("--poisson-frequency must be positive.");
    }
}

std::string make_default_output_directory(const std::string& mesh_file) {
    return fmt::format("self_consistent_amc_3d_{}", std::filesystem::path(mesh_file).stem().string());
}

void add_default_contacts(uepm::device::device& simulation_device, uepm::mesh::mesh& mesh) {
    constexpr double eps_electrode    = 0.001;
    constexpr double ohmic_resistance = 0.0;

    const double min_x = mesh.get_bounding_box().get_x_min();
    const double max_x = mesh.get_bounding_box().get_x_max();
    const double min_y = mesh.get_bounding_box().get_y_min();
    const double max_y = mesh.get_bounding_box().get_y_max();
    const double min_z = mesh.get_bounding_box().get_z_min();
    const double max_z = mesh.get_bounding_box().get_z_max();

    const uepm::mesh::vector3 anode_corner1(min_x - eps_electrode, min_y - eps_electrode, min_z - eps_electrode);
    const uepm::mesh::vector3 anode_corner2(min_x + eps_electrode, max_y + eps_electrode, max_z + eps_electrode);
    simulation_device.add_contact("anode", anode_corner1, anode_corner2, ohmic_resistance);

    const uepm::mesh::vector3 cathode_corner1(max_x - eps_electrode, min_y - eps_electrode, min_z - eps_electrode);
    const uepm::mesh::vector3 cathode_corner2(max_x + eps_electrode, max_y + eps_electrode, max_z + eps_electrode);
    simulation_device.add_contact("cathode", cathode_corner1, cathode_corner2, ohmic_resistance);

    fmt::print("Added device contacts:\n");
    fmt::print("  anode   x in [{:.6e}, {:.6e}]\n", anode_corner1.x(), anode_corner2.x());
    fmt::print("  cathode x in [{:.6e}, {:.6e}]\n", cathode_corner1.x(), cathode_corner2.x());
}

}  // namespace

int main(int argc, const char** argv) {
    try {
        TCLAP::CmdLine cmd("3D self-consistent analytical Monte Carlo device simulation.", ' ', "1.0");

        TCLAP::ValueArg<std::string> arg_device_mesh("", "device-mesh", "Path to the 3D device mesh file.", true, "", "path");

        TCLAP::ValueArg<std::string> arg_material_file("",
                                                       "material-file",
                                                       "Path to the material parameter file used by the Poisson solver.",
                                                       false,
                                                       std::string(PROJECT_SRC_DIR) + "/parameter_files/materials-chel.yaml",
                                                       "path");

        TCLAP::ValueArg<std::string> arg_material("m",
                                                  "material",
                                                  "AMC transport material. Currently only Si is supported.",
                                                  false,
                                                  "Si",
                                                  "string");

        TCLAP::ValueArg<std::string> arg_output_dir("d", "outdir", "Output directory.", false, "", "path");

        TCLAP::ValueArg<std::string> arg_simulation_name("", "name", "Simulation name.", false, "self_consistent_amc_3d", "string");

        TCLAP::ValueArg<double> arg_time("t", "time", "Final simulation time in seconds.", false, 1.0e-12, "s");

        TCLAP::ValueArg<double> arg_dt("", "dt", "Synchronized device Monte Carlo time step in seconds.", false, 1.0e-15, "s");

        TCLAP::ValueArg<double> arg_temperature("T", "temperature", "Lattice temperature in K.", false, 300.0, "K");

        TCLAP::ValueArg<double> arg_max_energy("e",
                                               "max-energy",
                                               "Maximum carrier energy in eV used to precompute gamma_max.",
                                               false,
                                               10.0,
                                               "eV");

        TCLAP::ValueArg<double> arg_gamma_safety("", "gamma-safety", "Safety factor applied to gamma_max.", false, 1.2, "double");

        TCLAP::ValueArg<std::size_t> arg_gamma_samples("",
                                                       "gamma-samples",
                                                       "Number of energy samples used to precompute gamma_max.",
                                                       false,
                                                       1000,
                                                       "integer");

        TCLAP::ValueArg<std::size_t> arg_poisson_frequency("",
                                                           "poisson-frequency",
                                                           "Number of transport steps between two Poisson updates.",
                                                           false,
                                                           10,
                                                           "integer");

        TCLAP::ValueArg<double> arg_anode_voltage("", "anode-voltage", "Dirichlet voltage applied to the anode.", false, 0.0, "V");

        TCLAP::ValueArg<double> arg_cathode_voltage("", "cathode-voltage", "Dirichlet voltage applied to the cathode.", false, 0.0, "V");

        TCLAP::ValueArg<double> arg_start_x("", "x0", "Initial particle x position in mesh units.", true, 0.0, "double");

        TCLAP::ValueArg<double> arg_start_y("", "y0", "Initial particle y position in mesh units.", false, 0.0, "double");

        TCLAP::ValueArg<double> arg_start_z("", "z0", "Initial particle z position in mesh units.", false, 0.0, "double");

        TCLAP::ValueArg<std::size_t> arg_number_electrons("", "nelectrons", "Initial number of electrons.", false, 1, "integer");

        TCLAP::ValueArg<std::size_t> arg_number_holes("", "nholes", "Initial number of holes.", false, 0, "integer");

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

        TCLAP::ValueArg<int> arg_nb_threads("j", "nthreads", "Number of threads requested by the simulation.", false, 1, "integer");

        TCLAP::ValueArg<int> arg_seed("", "seed", "Random seed.", false, 0, "integer");

        TCLAP::SwitchArg arg_disable_impact_ionization("", "disable-impact-ionization", "Disable impact-ionization computation.", false);

        TCLAP::SwitchArg arg_disable_particle_creation("",
                                                       "disable-particle-creation",
                                                       "Compute impact ionization but do not create electron-hole pairs.",
                                                       false);

        TCLAP::SwitchArg arg_keep_particle_history("H", "keep-particle-history", "Store full particle trajectories.", false);

        TCLAP::SwitchArg arg_export_time_steps("E", "export-time-steps", "Export particle state periodically.", false);

        TCLAP::ValueArg<int> arg_export_frequency("",
                                                  "export-frequency",
                                                  "Export one time step every N iterations.",
                                                  false,
                                                  100,
                                                  "integer");

        TCLAP::SwitchArg arg_keep_going_without_electrons("",
                                                          "keep-going-without-electrons",
                                                          "Do not stop when no electrons remain in the device.",
                                                          false);

        cmd.add(arg_device_mesh);
        cmd.add(arg_material_file);
        cmd.add(arg_material);
        cmd.add(arg_output_dir);
        cmd.add(arg_simulation_name);
        cmd.add(arg_time);
        cmd.add(arg_dt);
        cmd.add(arg_temperature);
        cmd.add(arg_max_energy);
        cmd.add(arg_gamma_safety);
        cmd.add(arg_gamma_samples);
        cmd.add(arg_poisson_frequency);
        cmd.add(arg_anode_voltage);
        cmd.add(arg_cathode_voltage);
        cmd.add(arg_start_x);
        cmd.add(arg_start_y);
        cmd.add(arg_start_z);
        cmd.add(arg_number_electrons);
        cmd.add(arg_number_holes);
        cmd.add(arg_max_particles);
        cmd.add(arg_avalanche_threshold);
        cmd.add(arg_nb_threads);
        cmd.add(arg_seed);
        cmd.add(arg_disable_impact_ionization);
        cmd.add(arg_disable_particle_creation);
        cmd.add(arg_keep_particle_history);
        cmd.add(arg_export_time_steps);
        cmd.add(arg_export_frequency);
        cmd.add(arg_keep_going_without_electrons);

        cmd.parse(argc, argv);

        const std::string mesh_file       = arg_device_mesh.getValue();
        const std::string material_symbol = arg_material.getValue();

        validate_material_symbol(material_symbol);

        uepm::amc::options_device_amc device_options;
        device_options.m_t_max                                = arg_time.getValue();
        device_options.m_time_step                            = arg_dt.getValue();
        device_options.m_lattice_temperature                  = arg_temperature.getValue();
        device_options.m_max_energy_eV                        = arg_max_energy.getValue();
        device_options.m_self_scattering_safety_factor        = arg_gamma_safety.getValue();
        device_options.m_gamma_max_energy_samples             = arg_gamma_samples.getValue();
        device_options.m_max_number_particle                  = arg_max_particles.getValue();
        device_options.m_avalanche_threshold                  = arg_avalanche_threshold.getValue();
        device_options.m_nb_threads                           = arg_nb_threads.getValue();
        device_options.m_activate_impact_ionization           = !arg_disable_impact_ionization.getValue();
        device_options.m_particle_creation_activated          = !arg_disable_particle_creation.getValue();
        device_options.m_stop_simu_when_no_electron_remaining = !arg_keep_going_without_electrons.getValue();
        device_options.m_keep_particles_history               = arg_keep_particle_history.getValue();
        device_options.m_export_time_step                     = arg_export_time_steps.getValue();
        device_options.m_frequency_export_trajectory          = arg_export_frequency.getValue();

        uepm::amc::options_self_consistent_device_amc_3d self_consistent_options;
        self_consistent_options.m_poisson_frequency = arg_poisson_frequency.getValue();
        self_consistent_options.m_anode_voltage     = arg_anode_voltage.getValue();
        self_consistent_options.m_cathode_voltage   = arg_cathode_voltage.getValue();

        validate_device_options(device_options);
        validate_self_consistent_options(self_consistent_options);

        const std::string output_dir =
            arg_output_dir.getValue().empty() ? make_default_output_directory(mesh_file) : arg_output_dir.getValue();

        const std::string trajectory_dir = fmt::format("{}/trajectory", output_dir);

        std::filesystem::create_directories(output_dir);

        if (device_options.m_export_time_step || device_options.m_keep_particles_history) {
            std::filesystem::create_directories(trajectory_dir);
        }

        fmt::print("Loading mesh: {}\n", mesh_file);

        uepm::file::msh_file msh_file(mesh_file);
        msh_file.read_mesh();
        msh_file.read_states();

        uepm::mesh::mesh* mesh = msh_file.get_p_mesh();

        if (mesh == nullptr) {
            throw std::runtime_error("Mesh loading failed.");
        }

        if (mesh->get_dimension() != 3) {
            throw std::runtime_error("This app requires a 3D mesh.");
        }

        fmt::print("Loading materials: {}\n", arg_material_file.getValue());
        std::string material_file = PROJECT_SRC_DIR + std::string("/examples/materials/materials.yaml");
        uepm::physic::material::list_materials list_of_materials;
        list_of_materials.load_materials_from_file(material_file);
        uepm::device::device simulation_device(mesh);
        add_default_contacts(simulation_device, *mesh);

        const auto starting_position = make_vector3(arg_start_x.getValue(), arg_start_y.getValue(), arg_start_z.getValue());
        fmt::print("Self-consistent AMC 3D simulation\n");
        fmt::print("  mesh vertices: {}\n", mesh->get_nb_vertices());
        fmt::print("  output directory: {}\n", output_dir);
        fmt::print("  material: {}\n", material_symbol);
        fmt::print("  final time: {:.6e} s\n", device_options.m_t_max);
        fmt::print("  time step: {:.6e} s\n", device_options.m_time_step);
        fmt::print("  Poisson frequency: {}\n", self_consistent_options.m_poisson_frequency);
        fmt::print("  anode voltage: {:.6e} V\n", self_consistent_options.m_anode_voltage);
        fmt::print("  cathode voltage: {:.6e} V\n", self_consistent_options.m_cathode_voltage);
        fmt::print("  initial electrons: {}\n", arg_number_electrons.getValue());
        fmt::print("  initial holes: {}\n", arg_number_holes.getValue());
        fmt::print("  initial position: ({:.6e}, {:.6e}, {:.6e})\n", starting_position.x(), starting_position.y(), starting_position.z());

        uepm::amc::self_consistent_device_amc_simulation_3d simulation(simulation_device,
                                                                       device_options,
                                                                       self_consistent_options,
                                                                       list_of_materials,
                                                                       arg_simulation_name.getValue(),
                                                                       starting_position,
                                                                       arg_number_electrons.getValue(),
                                                                       arg_number_holes.getValue(),
                                                                       arg_seed.getValue());

        simulation.set_prefix_export_trajectory_filename(fmt::format("{}/time_step", trajectory_dir));

        const auto start = std::chrono::high_resolution_clock::now();

        simulation.run_self_consistent_transport_simulation();

        const auto                          stop    = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> elapsed = stop - start;

        fmt::print("Simulation completed in {:.3f} s\n", elapsed.count());
        fmt::print("Remaining electrons: {}\n", simulation.get_number_electrons());
        fmt::print("Remaining holes: {}\n", simulation.get_number_holes());

        const std::string history_file = fmt::format("{}/device_history.csv", output_dir);

        simulation.export_history_to_csv(history_file);

        fmt::print("Wrote {}\n", history_file);

        if (device_options.m_keep_particles_history) {
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
