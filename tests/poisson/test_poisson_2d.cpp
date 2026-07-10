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

#include "doctest/doctest.h"
#include "msh_file.hpp"
#include "poisson_solver_2d.hpp"
#include "vtkWriter.hpp"

TEST_CASE("Testing Poisson 2d on a NMOS.") {
    //     std::string mesh_file_path = "./nmos_test.msh";
    //     fmt::print("Loading mesh: {}\n", mesh_file_path);

    //     uepm::file::msh_file msh_file(mesh_file_path);
    //     msh_file.read_mesh();
    //     msh_file.read_states();

    //     uepm::mesh::mesh* mesh = msh_file.get_p_mesh();
    //     if (mesh == nullptr) {
    //         throw std::runtime_error("Mesh loading failed.");
    //     }

    //     const int mesh_dimension = mesh->get_dimension();
    //     if (mesh_dimension != 2 && mesh_dimension != 3) {
    //         throw std::runtime_error("Only 2D and 3D meshes are supported.");
    //     }

    // }
    //     uepm::PBMC::pbmc_transport_config cfg;
    //     cfg.m_carrier_type                     = carrier_type;
    //     cfg.m_lattice_temperature              = options.m_lattice_temperature;
    //     cfg.m_max_energy_eV                    = options.m_max_energy_eV;
    //     cfg.m_self_scattering_safety_factor    = options.m_self_scattering_safety_factor;
    //     cfg.m_gamma_max_energy_samples         = options.m_gamma_max_energy_samples;
    //     cfg.m_enable_impact_ionization         = options.m_activate_impact_ionization;
    //     cfg.m_enable_impurity_scattering       = options.m_enable_impurity_scattering;
    //     cfg.m_impurity_density_source          = impurity_density_source::particle_local;
    //     cfg.m_background_impurity_density_cm_3 = 0.0;
    //     if (options.m_enable_impurity_scattering) {
    //         cfg.m_impurity_scattering_model = options.m_impurity_scattering_model;
    //         cfg.m_impurity_screening_model  = options.m_impurity_screening_model;
    //     }

    //     fmt::print("Mesh dimension: {}D\n", mesh_dimension);
    //     const uepm::physics::material_repository material_repository =
    //         config.material_root.empty() ? uepm::physics::material_repository{}
    //                                      : uepm::physics::material_repository{config.material_root};
    //     fmt::print("Loading materials from repository: {}\n", material_repository.root().string());

    //     uepm::physics::material_database material_database = material_repository.load_all_materials();
    //     const auto&                      common_material   = material_database.require(config.material_symbol);
    //     device_options.m_material_model                    = load_pbmc_material_model(material_repository,
    //     common_material);

    //     uepm::device::device simulation_device(mesh);
    //     // add_default_PN_contacts(simulation_device, *mesh);
}