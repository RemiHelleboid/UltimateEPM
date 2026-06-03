/**
 * @file electron_phonon.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-02-10
 *
 * @copyright Copyright (c) 2024
 *
 */
/**
 * @file electron_phonon.cpp
 * @brief Analytical bulk Monte Carlo simulation driver.
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

#include "bulk_amc_simulation.hpp"

namespace {

uepm::amc::particle_type parse_particle_type(const std::string& value) {
    if (value == "electron") {
        return uepm::amc::particle_type::electron;
    }

    if (value == "hole") {
        return uepm::amc::particle_type::hole;
    }

    throw std::invalid_argument("invalid particle type");
}

void validate_material(const std::string& material_symbol) {
    if (material_symbol != "Si") {
        throw std::invalid_argument("only Si is currently supported by the analytical AMC model");
    }
}

}  // namespace

int main(int argc, const char** argv) {
    try {
        TCLAP::CmdLine cmd("Analytical bulk Monte Carlo simulation.", ' ', "1.0");

        TCLAP::ValueArg<std::string> arg_material("m",
                                                  "material",
                                                  "Material symbol. Currently only Si is supported.",
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

        TCLAP::SwitchArg arg_enable_impact_ionization("",
                                                      "enable-impact-ionization",
                                                      "Enable impact ionization scattering.",
                                                      false);
        TCLAP::SwitchArg arg_enable_impurity_scattering("",
                                                        "enable-impurity-scattering",
                                                        "Enable impurity scattering.",
                                                        false);

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
        const std::string output_dir = [&] {
            const std::string requested_output_dir = arg_output_dir.getValue();
            if (!requested_output_dir.empty()) {
                return requested_output_dir;
            }
            return fmt::format("bulk_amc_{}_{}_T{:.1f}_Ex{:.3e}",
                               material_symbol,
                               particle_type_string,
                               temperature,
                               arg_electric_field_x.getValue());
        }();
        std::filesystem::create_directories(output_dir);

        constexpr double V_per_cm_to_V_per_m = 100.0;

        uepm::amc::bulk_amc_simulation_config config;
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

        fmt::print("Running bulk AMC simulation\n");
        fmt::print("Material: {}\n", material_symbol);
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
        fmt::print("Gamma max energy: {:.6f} eV\n", config.m_max_energy_eV);
        fmt::print("Gamma safety factor: {:.6f}\n", config.m_self_scattering_safety_factor);
        fmt::print("Gamma samples: {}\n", config.m_gamma_max_energy_samples);
        fmt::print("Warmup fraction: {:.6f}\n", config.m_warmup_fraction);
        fmt::print("Output directory: {}\n", output_dir);

        uepm::amc::bulk_amc_simulation simulation{config};
        simulation.initialize();

        const auto start_time = std::chrono::high_resolution_clock::now();

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

        return 0;
    } catch (const TCLAP::ArgException& error) {
        fmt::print(stderr, "Argument error: {} for argument {}\n", error.error(), error.argId());
        return 1;
    } catch (const std::exception& error) {
        fmt::print(stderr, "Error: {}\n", error.what());
        return 1;
    }
}
