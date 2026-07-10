/**
 * @file pbmc_rates.cpp
 * @brief Export analytical PBMC electron/hole phonon scattering rates vs energy.
 */

#include <fmt/format.h>
#include <tclap/CmdLine.h>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "pbmc_material_model.hpp"
#include "pbmc_particle.hpp"
#include "pbmc_scattering_model.hpp"

namespace {

uepm::PBMC::particle_type parse_particle_type(const std::string& value) {
    if (value == "electron") {
        return uepm::PBMC::particle_type::electron;
    }
    if (value == "hole") {
        return uepm::PBMC::particle_type::hole;
    }
    throw std::invalid_argument("part-type must be electron or hole");
}

void ensure_parent_directory(const std::filesystem::path& filename) {
    const auto parent = filename.parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
}

std::vector<double> make_energy_grid(double max_energy_eV, double energy_step_eV) {
    if (max_energy_eV < 0.0) {
        throw std::invalid_argument("max-energy must be non-negative");
    }
    if (energy_step_eV <= 0.0) {
        throw std::invalid_argument("energy-step must be positive");
    }

    std::vector<double> energies;
    for (double energy = 0.0; energy <= max_energy_eV + 0.5 * energy_step_eV; energy += energy_step_eV) {
        energies.push_back(std::min(energy, max_energy_eV));
        if (energies.back() == max_energy_eV) {
            break;
        }
    }
    return energies;
}

void write_electron_rates(const uepm::PBMC::pbmc_material_model& material,
                          std::size_t                            valley_index,
                          double                                 temperature_K,
                          const std::vector<double>&             energies_eV,
                          const std::filesystem::path&           output_file) {
    if (valley_index >= material.m_electron_valleys.size()) {
        throw std::out_of_range("electron valley index is outside the PBMC material model");
    }

    const auto& valley = material.m_electron_valleys[valley_index];

    ensure_parent_directory(output_file);
    std::ofstream out(output_file);
    if (!out) {
        throw std::runtime_error(fmt::format("failed to open {}", output_file.string()));
    }

    out << "energy_eV,acoustic_s_1,intervalley_absorption_s_1,intervalley_emission_s_1,total_s_1";
    for (const auto& branch : material.m_electron_intervalley_transitions) {
        out << ',' << branch.m_name << "_abs_s_1," << branch.m_name << "_em_s_1";
    }
    out << '\n';

    for (const double energy_eV : energies_eV) {
        const double acoustic =
            uepm::PBMC::acoustic_scattering_rate(valley, material.m_electron_acoustic, energy_eV, temperature_K);

        std::vector<std::pair<double, double>> branch_rates;
        branch_rates.reserve(material.m_electron_intervalley_transitions.size());

        double intervalley_absorption = 0.0;
        double intervalley_emission   = 0.0;
        for (const auto& branch : material.m_electron_intervalley_transitions) {
            const double absorption =
                uepm::PBMC::intervalley_scattering_rate(valley,
                                                        branch,
                                                        material.m_electron_acoustic.mass_density_kg_per_m3,
                                                        energy_eV,
                                                        true,
                                                        temperature_K);
            const double emission =
                uepm::PBMC::intervalley_scattering_rate(valley,
                                                        branch,
                                                        material.m_electron_acoustic.mass_density_kg_per_m3,
                                                        energy_eV,
                                                        false,
                                                        temperature_K);
            intervalley_absorption += absorption;
            intervalley_emission += emission;
            branch_rates.emplace_back(absorption, emission);
        }

        const double total = acoustic + intervalley_absorption + intervalley_emission;
        out << fmt::format("{:.12g},{:.12e},{:.12e},{:.12e},{:.12e}",
                           energy_eV,
                           acoustic,
                           intervalley_absorption,
                           intervalley_emission,
                           total);
        for (const auto& [absorption, emission] : branch_rates) {
            out << fmt::format(",{:.12e},{:.12e}", absorption, emission);
        }
        out << '\n';
    }
}

void write_hole_rates(const uepm::PBMC::pbmc_material_model& material,
                      std::size_t                            band_index,
                      double                                 temperature_K,
                      const std::vector<double>&             energies_eV,
                      const std::filesystem::path&           output_file) {
    if (band_index >= material.m_hole_bands.size()) {
        throw std::out_of_range("hole band index is outside the PBMC material model");
    }

    const auto& band = material.m_hole_bands[band_index];

    std::vector<uepm::PBMC::hole_optical_transition> transitions;
    for (const auto& transition : material.m_hole_optical_transitions) {
        if (transition.initial_band == band_index) {
            transitions.push_back(transition);
        }
    }

    ensure_parent_directory(output_file);
    std::ofstream out(output_file);
    if (!out) {
        throw std::runtime_error(fmt::format("failed to open {}", output_file.string()));
    }

    out << "energy_eV,acoustic_s_1,optical_absorption_s_1,optical_emission_s_1,total_s_1";
    for (const auto& transition : transitions) {
        out << ',' << transition.name << "_abs_s_1," << transition.name << "_em_s_1";
    }
    out << '\n';

    for (const double energy_eV : energies_eV) {
        const double acoustic =
            uepm::PBMC::acoustic_scattering_rate(band, material.m_hole_acoustic, energy_eV, temperature_K);

        std::vector<std::pair<double, double>> transition_rates;
        transition_rates.reserve(transitions.size());

        double optical_absorption = 0.0;
        double optical_emission   = 0.0;
        for (const auto& transition : transitions) {
            const auto&  final_band = material.m_hole_bands.at(transition.final_band);
            const double absorption =
                uepm::PBMC::optical_scattering_rate_holes(final_band,
                                                          transition,
                                                          material.m_hole_acoustic.mass_density_kg_per_m3,
                                                          energy_eV,
                                                          true,
                                                          temperature_K);
            const double emission =
                uepm::PBMC::optical_scattering_rate_holes(final_band,
                                                          transition,
                                                          material.m_hole_acoustic.mass_density_kg_per_m3,
                                                          energy_eV,
                                                          false,
                                                          temperature_K);
            optical_absorption += absorption;
            optical_emission += emission;
            transition_rates.emplace_back(absorption, emission);
        }

        const double total = acoustic + optical_absorption + optical_emission;
        out << fmt::format("{:.12g},{:.12e},{:.12e},{:.12e},{:.12e}",
                           energy_eV,
                           acoustic,
                           optical_absorption,
                           optical_emission,
                           total);
        for (const auto& [absorption, emission] : transition_rates) {
            out << fmt::format(",{:.12e},{:.12e}", absorption, emission);
        }
        out << '\n';
    }
}

}  // namespace

int main(int argc, const char** argv) {
    try {
        TCLAP::CmdLine cmd("Export analytical PBMC phonon scattering rates vs energy.", ' ', "1.0");

        TCLAP::ValueArg<std::string> arg_material("m", "material", "Material symbol.", false, "Si", "string");
        TCLAP::ValueArg<std::string> arg_pbmc_set("", "pbmc-set", "PBMC parameter set.", false, "default", "string");
        TCLAP::ValueArg<std::string> arg_particle_type("p",
                                                       "part-type",
                                                       "Carrier type: electron or hole.",
                                                       false,
                                                       "electron",
                                                       "string");
        TCLAP::ValueArg<int> arg_valley_or_band("v", "valley", "Electron valley or hole band index.", false, 0, "int");
        TCLAP::ValueArg<double> arg_temperature("T",
                                                "temperature",
                                                "Lattice temperature in kelvin.",
                                                false,
                                                300.0,
                                                "double");
        TCLAP::ValueArg<double> arg_max_energy("e",
                                               "max-energy",
                                               "Maximum carrier kinetic energy in eV.",
                                               false,
                                               1.0,
                                               "double");
        TCLAP::ValueArg<double> arg_energy_step("", "energy-step", "Energy grid step in eV.", false, 0.005, "double");
        TCLAP::ValueArg<std::string> arg_output("o", "out", "Output CSV filename.", false, "pbmc_rates.csv", "string");

        cmd.add(arg_material);
        cmd.add(arg_pbmc_set);
        cmd.add(arg_particle_type);
        cmd.add(arg_valley_or_band);
        cmd.add(arg_temperature);
        cmd.add(arg_max_energy);
        cmd.add(arg_energy_step);
        cmd.add(arg_output);

        cmd.parse(argc, argv);

        const auto carrier_type = parse_particle_type(arg_particle_type.getValue());
        const int  index        = arg_valley_or_band.getValue();
        if (index < 0) {
            throw std::invalid_argument("valley index must be non-negative");
        }

        const auto material = uepm::PBMC::load_pbmc_material_model(arg_material.getValue(), arg_pbmc_set.getValue());
        const auto energies = make_energy_grid(arg_max_energy.getValue(), arg_energy_step.getValue());
        const auto output   = std::filesystem::path(arg_output.getValue());

        if (carrier_type == uepm::PBMC::particle_type::electron) {
            write_electron_rates(material,
                                 static_cast<std::size_t>(index),
                                 arg_temperature.getValue(),
                                 energies,
                                 output);
        } else {
            write_hole_rates(material, static_cast<std::size_t>(index), arg_temperature.getValue(), energies, output);
        }

        fmt::print("Wrote PBMC rates to {}\n", output.string());
    } catch (const std::exception& e) {
        fmt::print(stderr, "Error: {}\n", e.what());
        return 1;
    }

    return 0;
}
