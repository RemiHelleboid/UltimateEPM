/**
 * @file admc_bulk.cpp
 * @brief Command-line bulk ADMC simulation.
 */

#include <fmt/core.h>
#include <tclap/CmdLine.h>

#include <chrono>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "bulk_admc_simulation.hpp"

namespace {

uepm::ADMC::carrier_type parse_carrier_type(const std::string& value) {
    if (value == "electron") {
        return uepm::ADMC::carrier_type::electron;
    }
    if (value == "hole") {
        return uepm::ADMC::carrier_type::hole;
    }
    throw std::invalid_argument("carrier type must be 'electron' or 'hole'");
}

void export_observables(const std::filesystem::path&            filename,
                        const uepm::ADMC::bulk_admc_simulation& simulation,
                        uepm::ADMC::carrier_type                carrier_type) {
    std::ofstream stream(filename);
    if (!stream.is_open()) {
        throw std::runtime_error("cannot open output file: " + filename.string());
    }

    const auto&  config    = simulation.config();
    const auto&  diffusion = simulation.diffusion_observables();
    const auto&  particles = simulation.particles();
    const auto&  state     = particles.front().state();
    const double time_s    = simulation.current_time_s();
    const double drift_x   = diffusion.mean_displacement_m[0] / time_s;
    const double drift_y   = diffusion.mean_displacement_m[1] / time_s;
    const double drift_z   = diffusion.mean_displacement_m[2] / time_s;

    stream << "carrier_type,number_particles,time_s,time_step_s,temperature_K,doping_concentration_cm_3,"
              "electric_field_x_V_per_m,electric_field_y_V_per_m,electric_field_z_V_per_m,"
              "mobility_m2_per_V_s,einstein_diffusion_m2_per_s,"
              "measured_drift_x_m_per_s,measured_drift_y_m_per_s,measured_drift_z_m_per_s,"
              "measured_diffusion_x_m2_per_s,measured_diffusion_y_m2_per_s,measured_diffusion_z_m2_per_s\n";
    stream << uepm::ADMC::carrier_type_name(carrier_type) << ',' << particles.size() << ',' << time_s << ','
           << config.time_step_s << ',' << config.environment.lattice_temperature_K << ','
           << config.environment.doping_concentration_cm_3 << ',' << config.environment.electric_field_V_per_m.x()
           << ',' << config.environment.electric_field_V_per_m.y() << ','
           << config.environment.electric_field_V_per_m.z() << ',' << state.mobility_m2_per_V_s << ','
           << state.diffusion_m2_per_s << ',' << drift_x << ',' << drift_y << ',' << drift_z << ','
           << diffusion.coefficient_m2_per_s[0] << ',' << diffusion.coefficient_m2_per_s[1] << ','
           << diffusion.coefficient_m2_per_s[2] << '\n';
}

}  // namespace

int main(int argc, const char** argv) {
    try {
        TCLAP::CmdLine               cmd("Bulk advection-diffusion Monte Carlo simulation.", ' ', "1.0");
        TCLAP::ValueArg<std::string> arg_output("d", "outdir", "Output directory.", false, "bulk_admc", "string");
        TCLAP::ValueArg<std::string> arg_carrier("p",
                                                 "part-type",
                                                 "Carrier type: electron or hole.",
                                                 false,
                                                 "electron",
                                                 "string");
        TCLAP::ValueArg<int>         arg_particles("N", "npart", "Number of particles.", false, 10000, "int");
        TCLAP::ValueArg<double> arg_time("t", "time", "Final simulation time in seconds.", false, 1.0e-11, "double");
        TCLAP::ValueArg<double> arg_dt("", "dt", "Time step in seconds.", false, 1.0e-14, "double");
        TCLAP::ValueArg<double> arg_temperature("T",
                                                "temperature",
                                                "Lattice temperature in K.",
                                                false,
                                                300.0,
                                                "double");
        TCLAP::ValueArg<double> arg_doping("",
                                           "doping",
                                           "Absolute doping concentration in cm^-3.",
                                           false,
                                           1.0e16,
                                           "double");
        TCLAP::ValueArg<double> arg_ex("", "Ex", "Electric field x component in V/cm.", false, 0.0, "double");
        TCLAP::ValueArg<double> arg_ey("", "Ey", "Electric field y component in V/cm.", false, 0.0, "double");
        TCLAP::ValueArg<double> arg_ez("", "Ez", "Electric field z component in V/cm.", false, 0.0, "double");
        TCLAP::ValueArg<std::uint64_t> arg_seed("", "seed", "Random-number seed.", false, 5489u, "uint64");
        TCLAP::SwitchArg arg_export_history("E", "export", "Export particle histories to CSV files.", false);

        cmd.add(arg_output);
        cmd.add(arg_carrier);
        cmd.add(arg_particles);
        cmd.add(arg_time);
        cmd.add(arg_dt);
        cmd.add(arg_temperature);
        cmd.add(arg_doping);
        cmd.add(arg_ex);
        cmd.add(arg_ey);
        cmd.add(arg_ez);
        cmd.add(arg_seed);
        cmd.add(arg_export_history);
        cmd.parse(argc, argv);

        if (arg_particles.getValue() < 2) {
            throw std::invalid_argument("directional diffusion requires at least two particles");
        }

        const auto                              carrier_type = parse_carrier_type(arg_carrier.getValue());
        uepm::ADMC::bulk_admc_simulation_config config;
        config.number_electrons =
            carrier_type == uepm::ADMC::carrier_type::electron ? static_cast<std::size_t>(arg_particles.getValue()) : 0;
        config.number_holes =
            carrier_type == uepm::ADMC::carrier_type::hole ? static_cast<std::size_t>(arg_particles.getValue()) : 0;
        config.time_step_s                           = arg_dt.getValue();
        config.final_time_s                          = arg_time.getValue();
        config.random_seed                           = arg_seed.getValue();
        config.record_history                        = arg_export_history.getValue();
        config.environment.lattice_temperature_K     = arg_temperature.getValue();
        config.environment.doping_concentration_cm_3 = arg_doping.getValue();
        constexpr double V_per_cm_to_V_per_m         = 100.0;
        config.environment.electric_field_V_per_m    = {arg_ex.getValue() * V_per_cm_to_V_per_m,
                                                        arg_ey.getValue() * V_per_cm_to_V_per_m,
                                                        arg_ez.getValue() * V_per_cm_to_V_per_m};

        const std::filesystem::path output_directory(arg_output.getValue());
        std::filesystem::create_directories(output_directory);

        uepm::ADMC::bulk_admc_simulation simulation(config);
        const auto                       start = std::chrono::steady_clock::now();
        simulation.run();
        const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

        export_observables(output_directory / "observables.csv", simulation, carrier_type);
        if (arg_export_history.getValue()) {
            const std::string file_prefix =
                (output_directory / ("simulation_results_" + std::to_string(std::time(nullptr)))).string();
            simulation.export_particles_history_to_csv(file_prefix);
        }

        const auto& diffusion = simulation.diffusion_observables();
        const auto& state     = simulation.particles().front().state();
        fmt::print("Bulk ADMC completed in {:.3f} s\n", elapsed);
        fmt::print("Carrier: {}, particles: {}, simulated time: {:.6e} s\n",
                   uepm::ADMC::carrier_type_name(carrier_type),
                   simulation.particles().size(),
                   simulation.current_time_s());
        fmt::print("Mobility: {:.6e} m^2/(V.s)\n", state.mobility_m2_per_V_s);
        fmt::print("Einstein diffusion: {:.6e} m^2/s\n", state.diffusion_m2_per_s);
        fmt::print("Measured diffusion: Dx={:.6e}, Dy={:.6e}, Dz={:.6e} m^2/s\n",
                   diffusion.coefficient_m2_per_s[0],
                   diffusion.coefficient_m2_per_s[1],
                   diffusion.coefficient_m2_per_s[2]);
        fmt::print("Results: {}\n", (output_directory / "observables.csv").string());
        return 0;
    } catch (const TCLAP::ArgException& error) {
        fmt::print(stderr, "Argument error: {} ({})\n", error.error(), error.argId());
    } catch (const std::exception& error) {
        fmt::print(stderr, "Error: {}\n", error.what());
    }
    return 1;
}
