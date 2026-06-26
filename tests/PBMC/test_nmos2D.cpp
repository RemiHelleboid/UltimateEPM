/**
 * @file test_poisson_2d.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-06
 *
 * @copyright Copyright (c) 2021
 *
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <fmt/core.h>

#include "device.hpp"
#include "doctest/doctest.h"
#include "materials.hpp"
#include "msh_file.hpp"
#include "pbmc_device_config.hpp"
#include "pbmc_device_runner.hpp"
#include "pbmc_device_setup.hpp"
#include "pbmc_material_model.hpp"
#include "pbmc_run_manifest.hpp"
#include "pbmc_scattering_model.hpp"
#include "pbmc_transport_kernel.hpp"
#include "poisson_solver_2d.hpp"
#include "vtkWriter.hpp"

namespace {

TEST_CASE("Testing Poisson 2d on a NMOS.") {
    std::string mesh_file_path = "./nmos_test.msh";
    fmt::print("Loading mesh: {}\n", mesh_file_path);

    uepm::file::msh_file msh_file(mesh_file_path);
    msh_file.read_mesh();
    msh_file.read_states();

    uepm::mesh::mesh* p_mesh = msh_file.get_p_mesh();
    if (p_mesh == nullptr) {
        throw std::runtime_error("Mesh loading failed.");
    }

    const int mesh_dimension = p_mesh->get_dimension();
    if (mesh_dimension != 2 && mesh_dimension != 3) {
        throw std::runtime_error("Only 2D and 3D meshes are supported.");
    }

    uepm::PBMC::pbmc_transport_config cfg;
    cfg.m_carrier_type                  = uepm::PBMC::particle_type::electron;
    cfg.m_lattice_temperature           = 300.0;
    cfg.m_max_energy_eV                 = 1.0;
    cfg.m_self_scattering_safety_factor = 1.2;
    cfg.m_gamma_max_energy_samples      = 100;
    cfg.m_enable_impact_ionization      = false;
    cfg.m_enable_impurity_scattering    = false;
    fmt::print("Mesh dimension: {}D\n", mesh_dimension);
    const uepm::physics::material_repository material_repository = uepm::physics::material_repository{};
    fmt::print("Loading materials from repository: {}\n", material_repository.root().string());

    uepm::physics::material_database material_database = material_repository.load_all_materials();
    uepm::fem::poisson_solver_2d     MyPoissonSolver(p_mesh, p_mesh->get_nb_vertices(), material_database);

    fmt::print("Computing stiffness matrix ...\n");
    MyPoissonSolver.compute_stiffness_matrix();
    fmt::print("Computing second member ...\n");
    MyPoissonSolver.compute_second_member(0.0);

    fmt::print("Applying Dirichlet conditions ...\n");
    MyPoissonSolver.apply_dirichlet_condition("source", 0.0);
    MyPoissonSolver.apply_dirichlet_condition("body", 0.0);
    MyPoissonSolver.apply_dirichlet_condition("drain", 0.1);
    MyPoissonSolver.apply_dirichlet_condition("gate", 0.0);

    fmt::print("Solving Poisson system ...\n");
    MyPoissonSolver.decompose_matrix();
    MyPoissonSolver.solve_system();
    fmt::print("Exporting solution to mesh ...\n");
    MyPoissonSolver.add_solution_to_mesh_functions("PoissonSolution", true);
    fmt::print("Test completed successfully.\n");
    MyPoissonSolver.print_infos();
    MyPoissonSolver.export_solution_csv("nmos_test_poisson_solution.csv");

    const std::string FileName = "nmos_test_poisson_solution.vtk";
    uepm::file::export_as_vtk(*(p_mesh), FileName);
    const std::string FileNameMSH = "nmos_test_poisson_solution.msh";
    msh_file.export_as_msh(FileNameMSH, {}, 1);
}

TEST_CASE("Testing PBMC Simulation 2d on a NMOS.") {
    std::string mesh_file_path = "./nmos_test.msh";
    fmt::print("Loading mesh: {}\n", mesh_file_path);

    // static const std::string file_input_test_msh = PROJECT_SRC_DIR +
    // std::string("/tests/test_data/config_nmos.yaml");
    static const std::string config_file_input_test_msh = std::string("./config_nmos.yaml");

    auto config = uepm::PBMC::load_device_pbmc_config(config_file_input_test_msh);
    // uepm::PBMC::run_self_consistent_device_pbmc_simulation(run_config);

    uepm::PBMC::validate_material_symbol(config.material_symbol);
    auto device_options               = config.device_options;
    device_options.m_simulation_name  = config.simulation_name;
    device_options.m_output_directory = config.output_dir;

    auto self_consistent_options_2d = config.self_consistent_options_2d;

    device_options.validate();
    self_consistent_options_2d.validate();

    const std::string output_dir =
        config.output_dir.empty() ? uepm::PBMC::make_default_output_directory(config.mesh_file) : config.output_dir;
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
    device_options.m_material_model = uepm::PBMC::load_pbmc_material_model(material_repository, common_material);

    uepm::device::device simulation_device(mesh);
    uepm::PBMC::add_collecting_contacts(simulation_device, *mesh, {"source", "drain"});

    config.self_consistent_options_2d.m_common.m_initialize_particles_from_doping = true;

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
    fmt::print("  built-in potential: {}\n",
               config.self_consistent_options_2d.m_common.m_enable_built_in_potential ? "enabled" : "disabled");
    if (config.self_consistent_options_2d.m_common.m_enable_built_in_potential) {
        fmt::print("  intrinsic concentration at {:.3f} K: {:.6e} cm^-3\n",
                   config.device_options.m_lattice_temperature,
                   uepm::PBMC::silicon_intrinsic_concentration_cm_3(config.device_options.m_lattice_temperature));
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

    const auto&                     common_options = self_consistent_options_2d.m_common;
    uepm::PBMC::simulation_manifest manifest;

    uepm::PBMC::self_consistent_device_pbmc_simulation_2d simulation(simulation_device,
                                                                     device_options,
                                                                     config.self_consistent_options_2d,
                                                                     material_database,
                                                                     config.starting_position,
                                                                     config.number_electrons_start,
                                                                     config.number_holes_start,
                                                                     config.seed_random_generator);

    simulation.set_prefix_export_trajectory_filename(trajectory_dir);

    simulation.run_self_consistent_transport_simulation();
}
}  // namespace
