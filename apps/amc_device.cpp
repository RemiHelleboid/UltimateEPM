/**
 * @file amc_device.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-26
 *
 * @copyright Copyright (c) 2026
 *
 */

#include <tclap/CmdLine.h>

#include <string>

#include "amc_device_runner.hpp"
#include "amc_run_manifest.hpp"
#include "amc_device_setup.hpp"

int main(int argc, const char** argv) {
    try {
        TCLAP::CmdLine               cmd("Self-consistent analytical Monte Carlo device simulation.", ' ', "1.0");
        TCLAP::ValueArg<std::string> arg_device_mesh("",
                                                     "device-mesh",
                                                     "Path to the device mesh file.",
                                                     true,
                                                     "",
                                                     "path");
        TCLAP::ValueArg<std::string> arg_material_file(
            "",
            "material-file",
            "Path to the material parameter file used by the Poisson solver.",
            false,
            std::string(PROJECT_SRC_DIR) + "/examples/materials/materials.yaml",
            "path");
        TCLAP::ValueArg<std::string> arg_material("m",
                                                  "material",
                                                  "AMC transport material. Currently only Si is supported.",
                                                  false,
                                                  "Si",
                                                  "string");

        TCLAP::ValueArg<std::string> arg_output_dir("d", "outdir", "Output directory.", false, "", "path");
        TCLAP::ValueArg<std::string> arg_simulation_name("",
                                                         "name",
                                                         "Simulation name.",
                                                         false,
                                                         "self_consistent_amc",
                                                         "string");
        TCLAP::ValueArg<double>      arg_time("t", "time", "Final simulation time in seconds.", false, 1.0e-12, "s");
        TCLAP::ValueArg<double>      arg_dt("",
                                            "dt",
                                            "Synchronized device Monte Carlo time step in seconds.",
                                            false,
                                            1.0e-15,
                                            "s");
        TCLAP::ValueArg<double> arg_temperature("T", "temperature", "Lattice temperature in K.", false, 300.0, "K");
        TCLAP::ValueArg<double> arg_max_energy("e",
                                               "max-energy",
                                               "Maximum carrier energy in eV used to precompute gamma_max.",
                                               false,
                                               10.0,
                                               "eV");

        TCLAP::ValueArg<double> arg_gamma_safety("",
                                                 "gamma-safety",
                                                 "Safety factor applied to gamma_max.",
                                                 false,
                                                 1.2,
                                                 "double");

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

        TCLAP::ValueArg<double>      arg_anode_voltage("",
                                                       "anode-voltage",
                                                       "Dirichlet voltage applied to the anode.",
                                                       false,
                                                       0.0,
                                                       "V");
        TCLAP::ValueArg<double>      arg_cathode_voltage("",
                                                         "cathode-voltage",
                                                         "Dirichlet voltage applied to the cathode.",
                                                         false,
                                                         0.0,
                                                         "V");
        TCLAP::ValueArg<double>      arg_start_x("",
                                                 "x0",
                                                 "Initial particle x position in mesh units.",
                                                 false,
                                                 0.0,
                                                 "double");
        TCLAP::ValueArg<double>      arg_start_y("",
                                                 "y0",
                                                 "Initial particle y position in mesh units.",
                                                 false,
                                                 0.0,
                                                 "double");
        TCLAP::ValueArg<double>      arg_start_z("",
                                                 "z0",
                                                 "Initial particle z position in mesh units.",
                                                 false,
                                                 0.0,
                                                 "double");
        TCLAP::ValueArg<std::size_t> arg_number_electrons("",
                                                          "nelectrons",
                                                          "Initial number of electrons.",
                                                          false,
                                                          1,
                                                          "integer");
        TCLAP::ValueArg<std::size_t> arg_number_holes("", "nholes", "Initial number of holes.", false, 0, "integer");
        TCLAP::ValueArg<std::size_t> arg_max_particles("",
                                                       "max-particles",
                                                       "Hard maximum number of active particles.",
                                                       false,
                                                       1000000000,
                                                       "integer");

        TCLAP::ValueArg<double> arg_avalanche_voltage_drop(
            "",
            "avalanche-voltage-drop",
            "Absolute quench-circuit voltage drop used to detect avalanche, in volts.",
            false,
            1.0,
            "volts");
        TCLAP::ValueArg<double> arg_quench_high_field(
            "",
            "quench-high-field",
            "Particle electric-field threshold that resets the successful-quench quiet window, in V/cm.",
            false,
            1.0e5,
            "V/cm");
        TCLAP::ValueArg<double> arg_quench_quiet_time(
            "",
            "quench-quiet-time",
            "Required time without high-field particles or impact ionization after avalanche, in seconds.",
            false,
            1.0e-11,
            "seconds");
        TCLAP::ValueArg<int>         arg_nb_threads("j",
                                                    "nthreads",
                                                    "Number of threads requested by the simulation.",
                                                    false,
                                                    1,
                                                    "integer");

        TCLAP::ValueArg<int> arg_seed("", "seed", "Random seed.", false, 0, "integer");
        TCLAP::SwitchArg     arg_disable_impact_ionization("",
                                                           "disable-impact-ionization",
                                                           "Disable impact-ionization computation.",
                                                           false);
        TCLAP::SwitchArg     arg_disable_particle_creation(
            "",
            "disable-particle-creation",
            "Compute impact ionization but do not create electron-hole pairs.",
            false);
        TCLAP::SwitchArg arg_keep_particle_history("H",
                                                   "keep-particle-history",
                                                   "Store full particle trajectories.",
                                                   false);

        TCLAP::SwitchArg arg_export_time_steps("E", "export-time-steps", "Export particle state periodically.", false);
        TCLAP::ValueArg<int>    arg_export_frequency("",
                                                     "export-frequency",
                                                     "Export one time step every N iterations.",
                                                     false,
                                                     100,
                                                     "integer");
        TCLAP::SwitchArg        arg_keep_going_without_electrons("",
                                                                 "keep-going-without-electrons",
                                                                 "Do not stop when no electrons remain in the device.",
                                                                 false);
        TCLAP::ValueArg<double> arg_effective_depth(
            "",
            "effective-depth",
            "Effective physical depth represented by a 2D simulation, in microns.",
            false,
            1.0,
            "um");
        TCLAP::ValueArg<double> arg_particle_z_period("",
                                                      "particle-z-period",
                                                      "Numerical periodic z length used by 2D particles, in microns.",
                                                      false,
                                                      1.0,
                                                      "um");

        TCLAP::SwitchArg arg_inject_particle("",
                                             "inject-particle",
                                             "Inject one scheduled particle during the simulation.",
                                             false);

        TCLAP::ValueArg<double> arg_inject_time("",
                                                "inject-time",
                                                "Scheduled injection time, in seconds.",
                                                false,
                                                0.0,
                                                "s");

        TCLAP::ValueArg<double> arg_inject_x("",
                                             "inject-x",
                                             "Scheduled injection x position, in microns.",
                                             false,
                                             0.0,
                                             "um");

        TCLAP::ValueArg<double> arg_inject_y("",
                                             "inject-y",
                                             "Scheduled injection y position, in microns.",
                                             false,
                                             0.0,
                                             "um");

        TCLAP::ValueArg<double> arg_inject_z("",
                                             "inject-z",
                                             "Scheduled injection z position, in microns.",
                                             false,
                                             0.0,
                                             "um");

        TCLAP::ValueArg<std::string> arg_inject_type("",
                                                     "inject-type",
                                                     "Scheduled injected particle type: electron, e, hole, or h.",
                                                     false,
                                                     "electron",
                                                     "type");

        TCLAP::ValueArg<double> arg_inject_weight("",
                                                  "inject-weight",
                                                  "Scheduled injected particle numerical weight.",
                                                  false,
                                                  1.0,
                                                  "weight");

        TCLAP::SwitchArg arg_disable_doping_init_particles("",
                                                           "disable-doping-init-particles",
                                                           "Disable initialization of particles from doping.",
                                                           false);

        TCLAP::SwitchArg arg_enable_impurity_scattering("",
                                                        "enable-impurity-scattering",
                                                        "Enable impurity scattering.",
                                                        false);

        TCLAP::ValueArg<std::string> arg_impurity_model("",
                                                        "impurity-model",
                                                        "Impurity scattering model: mobility or screened-coulomb.",
                                                        false,
                                                        "mobility",
                                                        "model");

        TCLAP::ValueArg<double> arg_quench_resistance("R", "resistance", "Quench resistance.", false, 1.0, "weight");

        TCLAP::ValueArg<double> arg_quench_capacitance("C", "capacitance", "Quench capacitance.", false, 1.0, "weight");

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
        cmd.add(arg_effective_depth);
        cmd.add(arg_particle_z_period);
        cmd.add(arg_start_x);
        cmd.add(arg_start_y);
        cmd.add(arg_start_z);
        cmd.add(arg_number_electrons);
        cmd.add(arg_number_holes);
        cmd.add(arg_max_particles);
        cmd.add(arg_avalanche_voltage_drop);
        cmd.add(arg_quench_high_field);
        cmd.add(arg_quench_quiet_time);
        cmd.add(arg_nb_threads);
        cmd.add(arg_seed);
        cmd.add(arg_disable_impact_ionization);
        cmd.add(arg_disable_particle_creation);
        cmd.add(arg_keep_particle_history);
        cmd.add(arg_export_time_steps);
        cmd.add(arg_export_frequency);
        cmd.add(arg_keep_going_without_electrons);
        cmd.add(arg_inject_particle);
        cmd.add(arg_inject_time);
        cmd.add(arg_inject_x);
        cmd.add(arg_inject_y);
        cmd.add(arg_inject_z);
        cmd.add(arg_inject_type);
        cmd.add(arg_inject_weight);
        cmd.add(arg_disable_doping_init_particles);
        cmd.add(arg_enable_impurity_scattering);
        cmd.add(arg_impurity_model);
        cmd.add(arg_quench_resistance);
        cmd.add(arg_quench_capacitance);

        cmd.parse(argc, argv);

        const std::string mesh_file       = arg_device_mesh.getValue();
        const std::string material_symbol = arg_material.getValue();

        uepm::amc::options_device_amc device_options;
        device_options.m_t_max                                = arg_time.getValue();
        device_options.m_time_step                            = arg_dt.getValue();
        device_options.m_lattice_temperature                  = arg_temperature.getValue();
        device_options.m_max_energy_eV                        = arg_max_energy.getValue();
        device_options.m_self_scattering_safety_factor        = arg_gamma_safety.getValue();
        device_options.m_gamma_max_energy_samples             = arg_gamma_samples.getValue();
        device_options.m_max_number_particle                  = arg_max_particles.getValue();
        device_options.m_nb_threads                           = arg_nb_threads.getValue();
        device_options.m_activate_impact_ionization           = !arg_disable_impact_ionization.getValue();
        device_options.m_particle_creation_activated          = !arg_disable_particle_creation.getValue();
        device_options.m_stop_simu_when_no_electron_remaining = !arg_keep_going_without_electrons.getValue();
        device_options.m_keep_particles_history               = arg_keep_particle_history.getValue();
        device_options.m_export_time_step                     = arg_export_time_steps.getValue();
        device_options.m_frequency_export_trajectory          = arg_export_frequency.getValue();
        device_options.m_enable_scheduled_particle_injection  = arg_inject_particle.getValue();
        device_options.m_output_directory                     = arg_output_dir.getValue();
        if (device_options.m_enable_scheduled_particle_injection) {
            auto& injection    = device_options.m_scheduled_particle_injection;
            injection.m_time_s = arg_inject_time.getValue();
            injection.m_position_um =
                uepm::mesh::vector3{arg_inject_x.getValue(), arg_inject_y.getValue(), arg_inject_z.getValue()};
            injection.m_particle_type = uepm::amc::parse_particle_type(arg_inject_type.getValue());
            injection.m_weight        = arg_inject_weight.getValue();
        }
        device_options.m_enable_impurity_scattering = arg_enable_impurity_scattering.getValue();
        device_options.m_impurity_scattering_model  = uepm::amc::parse_impurity_model(arg_impurity_model.getValue());

        device_options.validate();

        uepm::amc::options_self_consistent_device_amc_common sc_options;
        sc_options.m_poisson_frequency                                 = arg_poisson_frequency.getValue();
        sc_options.m_anode_voltage                                     = arg_anode_voltage.getValue();
        sc_options.m_cathode_voltage                                   = arg_cathode_voltage.getValue();
        sc_options.m_initialize_particles_from_doping                  = !arg_disable_doping_init_particles.getValue();
        sc_options.m_passive_quench_circuit.m_enabled                  = true;
        sc_options.m_passive_quench_circuit.m_bias_voltage_V           = sc_options.m_cathode_voltage;
        sc_options.m_passive_quench_circuit.m_initial_device_voltage_V = sc_options.m_cathode_voltage;
        sc_options.m_passive_quench_circuit.m_resistance_ohm           = arg_quench_resistance.getValue();
        sc_options.m_passive_quench_circuit.m_capacitance_F            = arg_quench_capacitance.getValue();
        sc_options.m_quench_biased_contact                             = uepm::amc::quench_biased_contact::cathode;
        sc_options.m_ramo_current_to_quench_current_sign               = -1.0;
        sc_options.m_avalanche_voltage_drop_threshold_V                 = arg_avalanche_voltage_drop.getValue();
        sc_options.m_quench_high_field_threshold_V_per_cm               = arg_quench_high_field.getValue();
        sc_options.m_quench_quiet_time_s                                = arg_quench_quiet_time.getValue();
        sc_options.validate();

        uepm::amc::self_consistent_device_amc_run_config run_config;
        run_config.mesh_file       = mesh_file;
        run_config.material_file   = arg_material_file.getValue();
        run_config.material_symbol = material_symbol;
        run_config.output_dir      = arg_output_dir.getValue();
        run_config.simulation_name = arg_simulation_name.getValue();
        run_config.command_line    = uepm::amc::command_line_from_arguments(argc, argv);
        run_config.starting_position =
            uepm::mesh::vector3{arg_start_x.getValue(), arg_start_y.getValue(), arg_start_z.getValue()};
        run_config.number_electrons_start                            = arg_number_electrons.getValue();
        run_config.number_holes_start                                = arg_number_holes.getValue();
        run_config.seed_random_generator                             = arg_seed.getValue();
        run_config.device_options                                    = device_options;
        run_config.self_consistent_options_2d.m_common               = sc_options;
        run_config.self_consistent_options_2d.m_effective_depth_um   = arg_effective_depth.getValue();
        run_config.self_consistent_options_2d.m_particle_z_period_um = arg_particle_z_period.getValue();
        run_config.self_consistent_options_3d.m_common               = sc_options;

        run_self_consistent_device_amc_simulation(run_config);

        return 0;
    } catch (const TCLAP::ArgException& error) {
        fmt::print(stderr, "Argument error: {} for argument {}\n", error.error(), error.argId());
        return 1;
    } catch (const std::exception& error) {
        fmt::print(stderr, "Error: {}\n", error.what());
        return 1;
    }
}
