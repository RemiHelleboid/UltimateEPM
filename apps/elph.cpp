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

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <tclap/CmdLine.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include "BandStructure.h"
#include "Material.h"
#include "Options.h"
#include "electron_phonon.hpp"
#include "fermi_level.hpp"

template <typename Derived>
struct fmt::formatter<Eigen::DenseBase<Derived>> : fmt::ostream_formatter {};

int export_result_mobility(const std::string     &filename,
                           const Eigen::Matrix3d &mu_tensor,
                           double                 mu_iso,
                           double                 Ef,
                           const Options         &my_options,
                           double                 max_energy,
                           double                 temperature,
                           std::size_t            nb_vtx,
                           std::size_t            nb_conduction_bands,
                           std::size_t            nb_valence_bands) {
    std::ofstream    file(filename);
    constexpr double mu_to_cm2Vs = 1e4;  // m^2/(V·s) to cm^2/(V·s)

    if (file) {
        file << "# Mobility tensor computed with EPP\n";
        file << "# Material : " << my_options.materialName << "\n";
        file << "# Number of vertices : " << nb_vtx << "\n";
        file << "# Number of conduction bands : " << nb_conduction_bands << "\n";
        file << "# Number of valence bands : " << nb_valence_bands << "\n";
        file << "# Energy range in eV : " << max_energy << "\n";
        file << "# Temperature in Kelvin : " << temperature << "\n";
        file << "# Fermi level in eV : " << Ef << "\n";
        file << "# Mobility tensor in cm^2/(V·s)\n";
        file << mu_tensor * mu_to_cm2Vs << "\n";
        file << "# Isotropic mobility in cm^2/(V·s)\n";
        file << mu_iso * mu_to_cm2Vs << "\n";
        file.close();
        fmt::print("Mobility tensor written to {}\n", filename);
        return 0;
    } else {
        fmt::print(std::cerr, "Error: could not write to file {}\n", filename);
        return 1;
    }
    return 0;
}

int main(int argc, char const *argv[]) {
    fmt::print("Starting UltimateEPM Electron-Phonon Calculations ... \n\n");

    TCLAP::CmdLine               cmd("EPP PROGRAM. COMPUTE BAND STRUCTURE ON A BZ MESH.", ' ', "1.0");
    TCLAP::ValueArg<std::string> arg_mesh_file("f", "meshbandfile", "File with BZ mesh and bands energy.", true, "bz.msh", "string");
    TCLAP::ValueArg<std::string> arg_phonon_rates("P", "phononrates", "File to load BZ phonon rates.", false, "bz_phonon.csv", "string");
    TCLAP::ValueArg<std::string> arg_material("m", "material", "Symbol of the material to use (Si, Ge, GaAs, ...)", true, "Si", "string");
    TCLAP::ValueArg<int>         arg_nb_energies("e", "nenergy", "Number of energies to compute", false, 250, "int");
    TCLAP::ValueArg<int>         arg_nb_conduction_bands("c", "ncbands", "Number of conduction bands to consider", false, -1, "int");
    TCLAP::ValueArg<int>         arg_nb_valence_bands("v", "nvbands", "Number of valence bands to consider", false, -1, "int");
    TCLAP::ValueArg<int>         arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    TCLAP::ValueArg<double>      arg_temperature("T", "temperature", "Temperature in Kelvin.", false, 300.0, "double");
    TCLAP::ValueArg<double>      arg_band_gap("g", "bandgap", "Band gap energy in eV.", false, 1.12, "double");
    TCLAP::ValueArg<double>      arg_energy_range("E",
                                             "energy_window",
                                             "Energy window around the band gap to consider (in eV).",
                                             false,
                                             0.3,
                                             "double");
    TCLAP::SwitchArg             arg_export_rates("X", "export-rates", "Export electron-phonon rates in k.", false);
    TCLAP::SwitchArg plot_with_python("p", "plot", "Call a python script after the computation to plot the band structure.", false);
    TCLAP::SwitchArg use_irr_wedge("w", "wedge", "Consider only the irreducible wedge of the BZ.", false);
    TCLAP::SwitchArg plot_with_knkpnp("K", "knkpnp", "Compute and store the full (n,k) -> (n',k') transition rate matrices.", false);
    TCLAP::SwitchArg use_unit_defpot("U", "unitdefpot", "Keep deformation potential to 1.0.", false);
    cmd.add(plot_with_python);
    cmd.add(arg_mesh_file);
    cmd.add(arg_material);
    cmd.add(arg_nb_conduction_bands);
    cmd.add(arg_nb_valence_bands);
    cmd.add(arg_nb_energies);
    cmd.add(arg_nb_threads);
    cmd.add(use_irr_wedge);
    cmd.add(plot_with_knkpnp);
    cmd.add(arg_temperature);
    cmd.add(arg_energy_range);
    cmd.add(arg_export_rates);
    cmd.add(arg_phonon_rates);
    cmd.add(use_unit_defpot);
    cmd.add(arg_band_gap);
    cmd.parse(argc, argv);

    auto start = std::chrono::high_resolution_clock::now();

    uepm::pseudopotential::Materials materials;
    const std::string                file_material_parameters = std::string(PROJECT_SRC_DIR) + "/parameter_files/materials-chel.yaml";
    materials.load_material_parameters(file_material_parameters);

    Options my_options;
    my_options.materialName                          = arg_material.getValue();
    my_options.nrLevels                              = arg_nb_conduction_bands.getValue() + arg_nb_valence_bands.getValue();
    my_options.nrThreads                             = arg_nb_threads.getValue();
    const int         number_energies                = arg_nb_energies.getValue();
    const int         nb_conduction_bands            = arg_nb_conduction_bands.getValue();
    const int         nb_valence_bands               = arg_nb_valence_bands.getValue();
    const double      max_energy                     = arg_energy_range.getValue();  // eV
    const double      temperature                    = arg_temperature.getValue();
    bool              irreducible_wedge_only         = use_irr_wedge.getValue();
    const std::string mesh_band_input_file           = arg_mesh_file.getValue();
    const std::string phonon_file                    = std::string(PROJECT_SRC_DIR) + "/parameter_files/phonon_kamakura.yaml";
    const bool        shift_conduction_band          = true;
    const bool        set_positive_valence_band      = false;
    const bool        export_rates                   = arg_export_rates.getValue();
    bool              use_unit_deformation_potential = use_unit_defpot.getValue();
    bool              phonon_rates_provided          = arg_phonon_rates.isSet();
    std::string       phonon_rates_file              = "";
    if (phonon_rates_provided) {
        phonon_rates_file = arg_phonon_rates.getValue();
    }
    double band_gap = arg_band_gap.getValue();

    uepm::pseudopotential::Material current_material = materials.materials.at(arg_material.getValue());

    uepm::mesh_bz::ElectronPhonon ElectronPhonon{current_material};
    ElectronPhonon.set_temperature(temperature);
    ElectronPhonon.set_number_threads_mesh_ops(my_options.nrThreads);
    ElectronPhonon.set_max_energy_global(max_energy);

    ElectronPhonon.read_mesh_geometry_from_msh_file(mesh_band_input_file);
    ElectronPhonon.build_search_tree();
    ElectronPhonon.read_mesh_bands_from_msh_file(mesh_band_input_file,
                                                 nb_conduction_bands,
                                                 nb_valence_bands,
                                                 shift_conduction_band,
                                                 set_positive_valence_band);
    ElectronPhonon.load_phonon_parameters(phonon_file);
    ElectronPhonon.set_nb_bands_elph(nb_conduction_bands);

    std::size_t           nb_vtx = ElectronPhonon.get_number_vertices();
    std::filesystem::path name_path(mesh_band_input_file);
    std::string           name_stem     = name_path.stem().string();
    auto                  stamp_params  = fmt::format("_T{}K_C{}V{}_N{}", temperature, nb_conduction_bands, nb_valence_bands, nb_vtx);
    std::string           prefix_export = name_stem + stamp_params;

    const double energy_windows_guard = 10.0 * uepm::constants::k_b_eV * temperature;
    if (max_energy < energy_windows_guard) {
        fmt::print("Warning: energy window {:.3f} eV is small compared to thermal energy scale {:.3f} eV at T = {:.1f} K.\n",
                   max_energy,
                   energy_windows_guard,
                   temperature);
    }
    if (phonon_rates_provided) {
        ElectronPhonon.read_phonon_scattering_rates_from_file(phonon_rates_file);
    } else {
        ElectronPhonon.compute_electron_phonon_rates_over_mesh(max_energy, irreducible_wedge_only);
    }

    if (export_rates && !phonon_rates_provided) {
        const std::string rates_file = prefix_export + "_eph_rates.msh";
        ElectronPhonon.export_rate_values(rates_file);
    }
    ElectronPhonon.test_elph();

    ElectronPhonon.apply_scissor(band_gap);  // eV
    // Solve for Fermi level and export CSV
    uepm::mesh_bz::fermi::Options fermi_options;
    fermi_options.nE                = 250;  // number of energy points for DOS interpolation
    fermi_options.threads           = my_options.nrThreads;
    fermi_options.use_interp        = false;        // use interpolation when computing DOS at given energy
    fermi_options.T_K               = temperature;  // temperature for Fermi-Dirac
    const bool use_iw               = true;         // use only irreducible wedge for DOS and Fermi level
    fermi_options.abs_max_energy_eV = 1.0;          // absolute max energy to consider (both conduction and valence)

    auto result = uepm::mesh_bz::fermi::solve_fermi(ElectronPhonon, fermi_options, use_iw);
    if (result.success) {
        fmt::print("Fermi level found: EF = {:.6f} eV\n", result.EF_eV);
        fmt::print("  p = {:.6e} cm^-3\n", result.p_m3 * 1e-6);
        fmt::print("  n = {:.6e} cm^-3\n", result.n_m3 * 1e-6);
    } else {
        fmt::print(std::cerr, "Fermi level not found.\n");
        return 1;
    }
    const double Ef = result.EF_eV;
    const double T  = temperature;

    const auto       mu_tensor   = ElectronPhonon.compute_electron_MRTA_mobility_tensor(Ef, T);
    const double     mu_iso      = ElectronPhonon.compute_electron_MRTA_mobility_isotropic(Ef, T);
    constexpr double mu_to_cm2Vs = 1e4;  // m^2/(V·s) to cm^2/(V·s)
    Eigen::Matrix3d  M           = mu_tensor * mu_to_cm2Vs;
    fmt::print("\n\nAt T = {:.1f} K and EF = {:.6f} eV:\n\n", temperature, Ef);
    fmt::print("μ_iso = {:.2f} cm^2/(V·s)\n\n", mu_iso * mu_to_cm2Vs);
    fmt::print("tensor = \n{} cm^2/(V*s)\n\n\n", fmt::streamed(M));

    double mean_energy = ElectronPhonon.mean_electron_energy_equilibrium(Ef, T, true);
    fmt::print("Mean electron energy above CBM at equilibrium: {:.6f} eV\n", mean_energy);

    std::string output_mobility     = name_stem + "_mobility.txt";
    auto        out                 = name_stem + stamp_params + "_mobility.txt";
    std::string out_rates_vs_energy = prefix_export + "_elph_results.csv";

    export_result_mobility(out_rates_vs_energy,
                           mu_tensor,
                           mu_iso,
                           Ef,
                           my_options,
                           max_energy,
                           temperature,
                           nb_vtx,
                           nb_conduction_bands,
                           nb_valence_bands);

    auto stop     = std::chrono::high_resolution_clock::now();
    auto duration = stop - start;
    fmt::print("\nTotal time : {:.2f} seconds\n\n\n", std::chrono::duration<double>(duration).count());

    if (plot_with_python.getValue()) {
        std::string  rates_vs_Energy_file = prefix_export + "_rates_vs_energy.csv";
        const double energy_step          = 0.001;  // energy step in eV
        ElectronPhonon.compute_plot_electron_phonon_rates_vs_energy_over_mesh(max_energy, energy_step, rates_vs_Energy_file);
        std::string command = "python3 " + std::string(PROJECT_SRC_DIR) + "/python/plots/plot_phonon_rate.py -f " + prefix_export;
        fmt::print("Running command: {}\n", command);
        int pyRes = std::system(command.c_str());
        if (pyRes != 0) {
            fmt::print(std::cerr, "Error: Python script returned non-zero exit code {}\n", pyRes);
        }
    }

    return 0;
}
