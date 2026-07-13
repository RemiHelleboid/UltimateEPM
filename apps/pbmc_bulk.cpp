/**
 * @file pbmc_bulk.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-08
 *
 * @copyright Copyright (c) 2026
 *
 */

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <tclap/CmdLine.h>

#include <chrono>
#include <ctime>
#include <filesystem>
#include <iostream>
#include <string>
#include <thread>

#include "bulk_pbmc_simulation.hpp"
#include "pbmc_device_setup.hpp"
#include "pbmc_run_manifest.hpp"

namespace {

uepm::PBMC::particle_type parse_particle_type(const std::string& value) {
    if (value == "electron") {
        return uepm::PBMC::particle_type::electron;
    }

    if (value == "hole") {
        return uepm::PBMC::particle_type::hole;
    }

    throw std::invalid_argument("invalid particle type");
}

void validate_material(const std::string& material_symbol) {
    if (material_symbol != "Si") {
        throw std::invalid_argument("only Si is currently supported by the analytical PBMC model");
    }
}

}  // namespace

int main(int argc, const char** argv) {
    try {
        TCLAP::CmdLine cmd("Analytical bulk Monte Carlo simulation.", ' ', "1.0");

        TCLAP::ValueArg<std::string> arg_material("m",
                                                  "material",
                                                  "epm_material symbol. Currently only Si is supported.",
                                                  false,
                                                  "Si",
                                                  "string");

        TCLAP::ValueArg<std::string> arg_output_dir("d",
                                                    "outdir",
                                                    "Output directory for results.",
                                                    false,
                                                    "",
                                                    "string");

        TCLAP::ValueArg<std::string> arg_particle_type("p",
                                                       "part-type",
                                                       "Carrier type: electron or hole.",
                                                       false,
                                                       "electron",
                                                       "string");

        TCLAP::ValueArg<std::string> arg_runner("",
                                                "runner",
                                                "Simulation runner: self-scattering or fixed-step.",
                                                false,
                                                "self-scattering",
                                                "string");

        TCLAP::ValueArg<int> arg_number_particles("N", "npart", "Number of particles to simulate.", false, 1, "int");

        TCLAP::ValueArg<int> arg_number_threads("j", "nthreads", "Number of threads to use.", false, 1, "int");

        TCLAP::ValueArg<double> arg_final_time("t",
                                               "time",
                                               "Simulation final time in seconds.",
                                               false,
                                               1.0e-12,
                                               "double");

        TCLAP::ValueArg<double> arg_time_step("",
                                              "dt",
                                              "Fixed time step used by the fixed-step runner in seconds.",
                                              false,
                                              5.0e-15,
                                              "double");

        TCLAP::ValueArg<double> arg_temperature("T",
                                                "temperature",
                                                "Lattice temperature in kelvin.",
                                                false,
                                                300.0,
                                                "double");

        TCLAP::ValueArg<double> arg_electric_field_x("",
                                                     "Ex",
                                                     "Electric field in x direction in V/cm.",
                                                     false,
                                                     0.0,
                                                     "double");

        TCLAP::ValueArg<double> arg_impurity_density(
            "",
            "impurity-density",
            "Background impurity density in cm^-3, used to compute impurity scattering rates.",
            false,
            1.0e10,
            "double");

        TCLAP::SwitchArg             arg_enable_impact_ionization("",
                                                      "enable-impact-ionization",
                                                      "Enable impact ionization scattering.",
                                                      false);
        TCLAP::SwitchArg             arg_enable_impurity_scattering("",
                                                        "enable-impurity-scattering",
                                                        "Enable impurity scattering.",
                                                        false);
        TCLAP::ValueArg<std::string> arg_impurity_model("",
                                                        "impurity-model",
                                                        "Impurity scattering model: mobility or screened-coulomb.",
                                                        false,
                                                        "mobility",
                                                        "string");
        TCLAP::ValueArg<std::string> arg_impurity_screening(
            "",
            "impurity-screening",
            "Screened-Coulomb screening: debye (analytic, default) or full (finite-temperature numerical).",
            false,
            "debye",
            "string");

        TCLAP::ValueArg<double> arg_max_energy("e",
                                               "max-energy",
                                               "Maximum carrier energy used to compute gamma_max in eV.",
                                               false,
                                               10.0,
                                               "double");

        TCLAP::ValueArg<double> arg_warmup_fraction("",
                                                    "warmup",
                                                    "Fraction of simulation time ignored for steady-state averages.",
                                                    false,
                                                    0.2,
                                                    "double");

        TCLAP::ValueArg<double> arg_gamma_safety("",
                                                 "gamma-safety",
                                                 "Safety factor applied to computed maximum self-scattering rate.",
                                                 false,
                                                 1.2,
                                                 "double");

        TCLAP::ValueArg<int> arg_gamma_samples("",
                                               "gamma-samples",
                                               "Number of energy samples used to compute gamma_max.",
                                               false,
                                               1000,
                                               "int");

        TCLAP::SwitchArg arg_plot_with_python("P", "plot", "Call a Python script after the Monte Carlo run.", false);

        TCLAP::SwitchArg arg_export_history("E", "export", "Export particle histories to CSV files.", false);

        cmd.add(arg_material);
        cmd.add(arg_output_dir);
        cmd.add(arg_particle_type);
        cmd.add(arg_runner);
        cmd.add(arg_number_particles);
        cmd.add(arg_number_threads);
        cmd.add(arg_final_time);
        cmd.add(arg_time_step);
        cmd.add(arg_temperature);
        cmd.add(arg_electric_field_x);
        cmd.add(arg_impurity_density);
        cmd.add(arg_enable_impact_ionization);
        cmd.add(arg_enable_impurity_scattering);
        cmd.add(arg_impurity_model);
        cmd.add(arg_impurity_screening);
        cmd.add(arg_max_energy);
        cmd.add(arg_warmup_fraction);
        cmd.add(arg_gamma_safety);
        cmd.add(arg_gamma_samples);
        cmd.add(arg_plot_with_python);
        cmd.add(arg_export_history);

        cmd.parse(argc, argv);

        const std::string material_symbol      = arg_material.getValue();
        const std::string particle_type_string = arg_particle_type.getValue();
        const std::string runner               = arg_runner.getValue();

        validate_material(material_symbol);

        const auto carrier_type = parse_particle_type(particle_type_string);

        const int number_particles = arg_number_particles.getValue();
        if (number_particles <= 0) {
            throw std::invalid_argument("number of particles must be positive");
        }
        const int gamma_samples = arg_gamma_samples.getValue();
        if (gamma_samples < 2) {
            throw std::invalid_argument("gamma-samples must be at least 2");
        }
        const double final_time = arg_final_time.getValue();
        if (final_time <= 0.0) {
            throw std::invalid_argument("simulation final time must be positive");
        }
        const double time_step = arg_time_step.getValue();
        if (time_step <= 0.0) {
            throw std::invalid_argument("time step must be positive");
        }
        const double temperature = arg_temperature.getValue();
        if (temperature < 0.0) {
            throw std::invalid_argument("temperature must be non-negative");
        }
        const double max_energy_eV = arg_max_energy.getValue();
        if (max_energy_eV <= 0.0) {
            throw std::invalid_argument("max energy must be positive");
        }
        const double warmup_fraction = arg_warmup_fraction.getValue();
        if (warmup_fraction < 0.0 || warmup_fraction >= 1.0) {
            throw std::invalid_argument("warmup fraction must be in [0, 1)");
        }
        const double gamma_safety = arg_gamma_safety.getValue();
        if (gamma_safety <= 0.0) {
            throw std::invalid_argument("gamma safety factor must be positive");
        }
        const std::size_t nb_threads = static_cast<std::size_t>(arg_number_threads.getValue());
        if (nb_threads == 0) {
            throw std::invalid_argument("number of threads must be positive");
        }
        const double impurity_density_cm_3 = arg_impurity_density.getValue();
        if (impurity_density_cm_3 < 0.0) {
            throw std::invalid_argument("impurity density must be non-negative");
        }
        const std::string                     impurity_model = arg_impurity_model.getValue();
        uepm::PBMC::impurity_scattering_model parsed_impurity_model;
        if (impurity_model == "mobility") {
            parsed_impurity_model = uepm::PBMC::impurity_scattering_model::mobility_empirical;
        } else if (impurity_model == "screened-coulomb") {
            parsed_impurity_model = uepm::PBMC::impurity_scattering_model::screened_coulomb;
        } else {
            throw std::invalid_argument("--impurity-model must be either 'mobility' or 'screened-coulomb'");
        }
        const auto parsed_impurity_screening =
            uepm::PBMC::parse_impurity_screening_model(arg_impurity_screening.getValue());
        const std::string output_dir = [&] {
            const std::string requested_output_dir = arg_output_dir.getValue();
            if (!requested_output_dir.empty()) {
                return requested_output_dir;
            }
            return std::string(fmt::format("bulk_pbmc_{}_{}", material_symbol, particle_type_string));
        }();
        std::filesystem::create_directories(output_dir);

        constexpr double V_per_cm_to_V_per_m = 100.0;

        uepm::PBMC::bulk_pbmc_simulation_config config;
        config.m_material_model                = uepm::PBMC::load_pbmc_material_model(material_symbol);
        config.m_carrier_type                  = carrier_type;
        config.m_record_history                = arg_export_history.getValue();
        config.m_lattice_temperature           = temperature;
        config.m_final_time                    = final_time;
        config.m_time_step                     = time_step;
        config.m_number_of_particles           = static_cast<std::size_t>(number_particles);
        config.m_electric_field                = {arg_electric_field_x.getValue() * V_per_cm_to_V_per_m, 0.0, 0.0};
        config.m_max_energy_eV                 = max_energy_eV;
        config.m_warmup_fraction               = warmup_fraction;
        config.m_self_scattering_safety_factor = gamma_safety;
        config.m_gamma_max_energy_samples      = static_cast<std::size_t>(gamma_samples);
        config.m_nb_threads                    = nb_threads;
        config.m_enable_impact_ionization      = arg_enable_impact_ionization.getValue();
        config.m_enable_impurity_scattering    = arg_enable_impurity_scattering.getValue();
        config.m_impurity_density_cm_3         = impurity_density_cm_3;
        config.m_impurity_scattering_model     = parsed_impurity_model;
        config.m_impurity_screening_model      = parsed_impurity_screening;

        fmt::print("Running bulk PBMC simulation\n");
        fmt::print("epm_material: {}\n", material_symbol);
        fmt::print("Carrier type: {}\n", particle_type_string);
        fmt::print("Runner: {}\n", runner);
        fmt::print("Number of particles: {}\n", config.m_number_of_particles);
        fmt::print("Temperature: {:.3f} K\n", config.m_lattice_temperature);
        fmt::print("Final time: {:.6e} s\n", config.m_final_time);
        fmt::print("Fixed time step: {:.6e} s\n", config.m_time_step);
        fmt::print("Electric field: ({:.6e}, {:.6e}, {:.6e}) V/m\n",
                   config.m_electric_field.x(),
                   config.m_electric_field.y(),
                   config.m_electric_field.z());
        fmt::print("  impurity scattering: {}\n", config.m_enable_impurity_scattering ? "enabled" : "disabled");
        fmt::print("  impurity density: {:.6e} cm^-3\n", config.m_impurity_density_cm_3);
        fmt::print(
            "  impurity model: {}\n",
            config.m_enable_impurity_scattering
                ? (config.m_impurity_scattering_model == uepm::PBMC::impurity_scattering_model::mobility_empirical
                       ? "mobility-empirical"
                       : "screened-coulomb")
                : "N/A");
        fmt::print("  impurity screening: {}\n",
                   config.m_enable_impurity_scattering &&
                           config.m_impurity_scattering_model == uepm::PBMC::impurity_scattering_model::screened_coulomb
                       ? uepm::PBMC::impurity_screening_model_name(config.m_impurity_screening_model)
                       : "N/A");
        fmt::print("Gamma max energy: {:.6f} eV\n", config.m_max_energy_eV);
        fmt::print("Gamma safety factor: {:.6f}\n", config.m_self_scattering_safety_factor);
        fmt::print("Gamma samples: {}\n", config.m_gamma_max_energy_samples);
        fmt::print("Warmup fraction: {:.6f}\n", config.m_warmup_fraction);
        fmt::print("Output directory: {}\n", output_dir);

        uepm::PBMC::bulk_pbmc_simulation simulation{config};
        simulation.initialize();

        const std::string started_at_utc = uepm::PBMC::current_utc_timestamp();
        const auto        start_time     = std::chrono::high_resolution_clock::now();

        if (runner == "self-scattering") {
            simulation.run_self_scattering_emc();
        } else if (runner == "fixed-step") {
            simulation.run();
        } else {
            throw std::invalid_argument("runner must be either self-scattering or fixed-step");
        }

        const auto                          end_time = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> elapsed  = end_time - start_time;

        fmt::print("Simulation completed in {:.3f} seconds.\n", elapsed.count());

        const std::string timestamp   = std::to_string(std::time(nullptr));
        const std::string file_prefix = fmt::format("{}/simulation_results_{}", output_dir, timestamp);

        if (arg_export_history.getValue()) {
            simulation.export_particles_history_to_csv(file_prefix);
        }

        simulation.export_observables_to_csv(fmt::format("{}/observables.csv", output_dir));

        if (arg_plot_with_python.getValue()) {
            fmt::print("Plot option requested, but no plotting hook is currently configured.\n");
        }

        const auto&  observables = simulation.observables();
        const auto&  ii_stats    = simulation.impact_ionization_statistics();
        const double mean_velocity_x_m_per_s =
            observables.accumulated_time_s > 0.0
                ? observables.weighted_velocity_x_m2_per_s2 / observables.accumulated_time_s
                : 0.0;
        const double mean_energy_eV = observables.accumulated_time_s > 0.0
                                          ? observables.weighted_kinetic_energy_eV_s / observables.accumulated_time_s
                                          : 0.0;

        uepm::PBMC::simulation_manifest manifest;
        manifest.add("run", "simulation_type", "bulk_PBMC");
        manifest.add("run", "status", "completed");
        manifest.add("run", "started_at_utc", started_at_utc);
        manifest.add("run", "finished_at_utc", uepm::PBMC::current_utc_timestamp());
        manifest.add("run", "elapsed_seconds", elapsed.count());
        manifest.add("run", "command_line", uepm::PBMC::command_line_from_arguments(argc, argv));
        manifest.add("run", "working_directory", std::filesystem::current_path().string());

        manifest.add("build", "project_version", uepm::PBMC::pbmc_project_version());
        manifest.add("build", "build_type", uepm::PBMC::pbmc_build_type());
        manifest.add("build", "compiler", uepm::PBMC::pbmc_compiler());
        manifest.add("build", "hardware_concurrency", static_cast<std::size_t>(std::thread::hardware_concurrency()));

        manifest.add("input", "material", material_symbol);
        manifest.add("input", "carrier_type", particle_type_string);
        manifest.add("input", "runner", runner);
        manifest.add("input", "output_directory", std::filesystem::absolute(output_dir).string());

        manifest.add("simulation", "number_of_particles", config.m_number_of_particles);
        manifest.add("simulation", "requested_threads", config.m_nb_threads);
        manifest.add("simulation",
                     "random_seed_policy",
                     runner == "self-scattering" ? "thread_seed_base_1234" : "random_device");
        manifest.add("simulation", "lattice_temperature_K", config.m_lattice_temperature);
        manifest.add("simulation", "final_time_s", config.m_final_time);
        manifest.add("simulation", "fixed_time_step_s", config.m_time_step);
        manifest.add("simulation", "warmup_fraction", config.m_warmup_fraction);
        manifest.add("simulation", "record_history", config.m_record_history);
        manifest.add("simulation", "electric_field_x_V_per_m", config.m_electric_field.x());
        manifest.add("simulation", "electric_field_y_V_per_m", config.m_electric_field.y());
        manifest.add("simulation", "electric_field_z_V_per_m", config.m_electric_field.z());

        manifest.add("scattering", "impact_ionization_enabled", config.m_enable_impact_ionization);
        manifest.add("scattering", "impurity_scattering_enabled", config.m_enable_impurity_scattering);
        manifest.add("scattering", "impurity_density_cm_3", config.m_impurity_density_cm_3);
        manifest.add("scattering",
                     "impurity_model",
                     config.m_impurity_scattering_model == uepm::PBMC::impurity_scattering_model::mobility_empirical
                         ? "mobility-empirical"
                         : "screened-coulomb");
        manifest.add("scattering",
                     "impurity_screening",
                     uepm::PBMC::impurity_screening_model_name(config.m_impurity_screening_model));
        manifest.add("scattering", "gamma_max_energy_eV", config.m_max_energy_eV);
        manifest.add("scattering", "gamma_safety_factor", config.m_self_scattering_safety_factor);
        manifest.add("scattering", "gamma_energy_samples", config.m_gamma_max_energy_samples);
        manifest.add("scattering", "computed_gamma_max_s_1", simulation.gamma_max());

        manifest.add("results", "mean_velocity_x_m_per_s", mean_velocity_x_m_per_s);
        manifest.add("results", "mean_kinetic_energy_eV", mean_energy_eV);
        manifest.add("results", "diffusion_x_m2_per_s", observables.diffusion.coefficient_m2_per_s[0]);
        manifest.add("results", "diffusion_y_m2_per_s", observables.diffusion.coefficient_m2_per_s[1]);
        manifest.add("results", "diffusion_z_m2_per_s", observables.diffusion.coefficient_m2_per_s[2]);
        manifest.add("results", "accumulated_carrier_time_s", observables.accumulated_time_s);
        manifest.add("results",
                     "acoustic_events",
                     simulation.count_scattering_events(uepm::PBMC::scattering_event::acoustic));
        manifest.add("results",
                     "intervalley_absorption_events",
                     simulation.count_scattering_events(uepm::PBMC::scattering_event::intervalley_absorption));
        manifest.add("results",
                     "intervalley_emission_events",
                     simulation.count_scattering_events(uepm::PBMC::scattering_event::intervalley_emission));
        manifest.add("results",
                     "impurity_events",
                     simulation.count_scattering_events(uepm::PBMC::scattering_event::impurity));
        manifest.add("results",
                     "impact_ionization_events",
                     simulation.count_scattering_events(uepm::PBMC::scattering_event::impact_ionization));
        manifest.add("results",
                     "self_scattering_events",
                     simulation.count_scattering_events(uepm::PBMC::scattering_event::self_scattering));
        manifest.add("results", "impact_ionization_rate_per_carrier_s_1", ii_stats.event_rate_per_carrier_s_1());
        manifest.add("results", "impact_ionization_coefficient_cm_1", ii_stats.ionization_coefficient_cm_1());

        manifest.add("outputs", "observables_csv", fmt::format("{}/observables.csv", output_dir));
        manifest.add("outputs", "particle_history_exported", config.m_record_history);
        if (config.m_record_history) {
            manifest.add("outputs", "particle_history_prefix", file_prefix);
        }

        const std::filesystem::path manifest_file = std::filesystem::path(output_dir) / "simulation_manifest.txt";
        manifest.add("outputs", "simulation_manifest", manifest_file.string());
        manifest.write(manifest_file);
        fmt::print("Wrote {}\n", manifest_file.string());

        return 0;
    } catch (const TCLAP::ArgException& error) {
        fmt::print(stderr, "Argument error: {} for argument {}\n", error.error(), error.argId());
        return 1;
    } catch (const std::exception& error) {
        fmt::print(stderr, "Error: {}\n", error.what());
        return 1;
    }
}
