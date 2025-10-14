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

#include <tclap/CmdLine.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "BandStructure.h"
#include "Material.h"
#include "Options.h"
#include "bz_mesh.hpp"
#include "bz_meshfile.hpp"
#include "bz_states.hpp"
#include "electron_phonon.hpp"
#include "fermi_level.hpp"

int main(int argc, char const *argv[]) {
    TCLAP::CmdLine               cmd("EPP PROGRAM. COMPUTE BAND STRUCTURE ON A BZ MESH.", ' ', "1.0");
    TCLAP::ValueArg<std::string> arg_mesh_file("f", "meshbandfile", "File with BZ mesh and bands energy.", true, "bz.msh", "string");
    TCLAP::ValueArg<std::string> arg_material("m", "material", "Symbol of the material to use (Si, Ge, GaAs, ...)", true, "Si", "string");
    TCLAP::ValueArg<int>         arg_nb_energies("e", "nenergy", "Number of energies to compute", false, 250, "int");
    TCLAP::ValueArg<int>         arg_nb_conduction_bands("c", "ncbands", "Number of conduction bands to consider", false, -1, "int");
    TCLAP::ValueArg<int>         arg_nb_valence_bands("v", "nvbands", "Number of valence bands to consider", false, -1, "int");
    TCLAP::ValueArg<int>         arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    TCLAP::ValueArg<double>      arg_temperature("T", "temperature", "Temperature in Kelvin.", false, 300.0, "double");
    TCLAP::ValueArg<double>      arg_energy_range("E",
                                             "energy_window",
                                             "Energy window around the band gap to consider (in eV).",
                                             false,
                                             0.3,
                                             "double");
    TCLAP::SwitchArg plot_with_python("P", "plot", "Call a python script after the computation to plot the band structure.", false);
    TCLAP::SwitchArg plot_with_wedge("w", "wedge", "Consider only the irreducible wedge of the BZ.", false);
    TCLAP::SwitchArg plot_with_knkpnp("K", "knkpnp", "Compute and store the full (n,k) -> (n',k') transition rate matrices.", false);
    cmd.add(plot_with_python);
    cmd.add(arg_mesh_file);
    cmd.add(arg_material);
    cmd.add(arg_nb_conduction_bands);
    cmd.add(arg_nb_valence_bands);
    cmd.add(arg_nb_energies);
    cmd.add(arg_nb_threads);
    cmd.add(plot_with_wedge);
    cmd.add(plot_with_knkpnp);
    cmd.add(arg_temperature);
    cmd.add(arg_energy_range);

    cmd.parse(argc, argv);

    auto start = std::chrono::high_resolution_clock::now();

    uepm::pseudopotential::Materials materials;
    const std::string                file_material_parameters = std::string(PROJECT_SRC_DIR) + "/parameter_files/materials-local.yaml";
    materials.load_material_parameters(file_material_parameters);

    Options my_options;
    my_options.materialName          = arg_material.getValue();
    my_options.nrLevels              = arg_nb_conduction_bands.getValue() + arg_nb_valence_bands.getValue();
    my_options.nrThreads             = arg_nb_threads.getValue();
    const int    number_energies     = arg_nb_energies.getValue();
    const int    nb_conduction_bands = arg_nb_conduction_bands.getValue();
    const int    nb_valence_bands    = arg_nb_valence_bands.getValue();
    const double max_energy          = arg_energy_range.getValue();  // eV
    const double temperature         = arg_temperature.getValue();

    uepm::pseudopotential::Material current_material = materials.materials.at(arg_material.getValue());

    const std::string             mesh_band_input_file = arg_mesh_file.getValue();
    uepm::mesh_bz::ElectronPhonon ElectronPhonon{current_material};
    ElectronPhonon.set_nb_threads(my_options.nrThreads);

    // const std::string phonon_file = std::string(PROJECT_SRC_DIR) + "/parameter_files/phonon_michaillat.yaml";
    const std::string phonon_file = std::string(PROJECT_SRC_DIR) + "/parameter_files/phonon_kamakura.yaml";

    ElectronPhonon.read_mesh_geometry_from_msh_file(mesh_band_input_file);
    std::cout << "Keeping " << nb_conduction_bands << " conduction bands and " << nb_valence_bands << " valence bands." << std::endl;
    const bool shift_conduction_band     = true;
    const bool set_positive_valence_band = false;
    ElectronPhonon.read_mesh_bands_from_msh_file(mesh_band_input_file,
                                                 nb_conduction_bands,
                                                 nb_valence_bands,
                                                 shift_conduction_band,
                                                 set_positive_valence_band);

    ElectronPhonon.set_nb_bands_elph(nb_conduction_bands);

    ElectronPhonon.apply_scissor(1.12);  // eV
    std::cout << std::scientific;
    // Solve for Fermi level and export CSV
    uepm::mesh_bz::fermi::Options fermi_options;
    fermi_options.nE         = 1000;  // number of energy points for DOS interpolation
    fermi_options.threads    = my_options.nrThreads;
    fermi_options.use_interp = false;   // use interpolation when computing DOS at given energy
    fermi_options.T_K        = 300.0;  // temperature for Fermi-Dirac

    auto result = uepm::mesh_bz::fermi::solve_fermi(ElectronPhonon, fermi_options);
    if (result.success) {
        std::cout << "Fermi level found: EF = " << result.EF_eV << " eV\n";
        std::cout << "  p = " << result.p_m3 * 1e-6 << " cm^-3\n";
        std::cout << "  n = " << result.n_m3 * 1e-6 << " cm^-3\n";
    } else {
        std::cout << "Fermi level not found.\n";
    }
    const double Ef        = result.EF_eV;
    const double T         = 300.0;

    ElectronPhonon.apply_scissor(-1.12);  // eV

    constexpr double energy_step_dos = 0.002;  // eV
    const double     max_energy_dos  = max_energy + 0.1;
    ElectronPhonon.precompute_dos_tetra(energy_step_dos, max_energy_dos);

    ElectronPhonon.load_phonon_parameters(phonon_file);
    bool irreducible_wedge_only = plot_with_wedge.getValue();
    bool populate_nk_npkp       = plot_with_knkpnp.getValue();

    std::cout << "Max energy: " << max_energy << " eV" << std::endl;
    const double energy_step = 0.001;  // eV
    ElectronPhonon.compute_electron_phonon_rates_over_mesh(max_energy, irreducible_wedge_only, populate_nk_npkp);
    ElectronPhonon.export_rate_values("rates_all.csv");

    ElectronPhonon.compute_plot_electron_phonon_rates_vs_energy_over_mesh(my_options.nrLevels,
                                                                          max_energy,
                                                                          energy_step,
                                                                          "rates_vs_energy.csv");


    ElectronPhonon.apply_scissor(1.12);  // eV
    const auto   mu_tensor = ElectronPhonon.compute_electron_MRTA_mobility_tensor(Ef, T);
    const double mu_iso    = ElectronPhonon.compute_electron_MRTA_mobility_isotropic(Ef, T);
    std::cout << "μ_iso = " << mu_iso << " m^2/(V·s)\n";
    std::cout << "μ tensor:\n" << mu_tensor << std::endl;

    // // ElectronPhonon.add_electron_phonon_rates_to_mesh(mesh_band_input_file, "rates.msh");
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n\n\n" << std::endl;

    if (plot_with_python.getValue()) {
        std::string command = "python3 " + std::string(PROJECT_SRC_DIR) + "/python/plots/plot_phonon_rate.py -f rates_vs_energy.csv";
        std::cout << "Running command: " << command << std::endl;
        std::system(command.c_str());
    }

    return 0;
}
