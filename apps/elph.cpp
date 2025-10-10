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


int main(int argc, char const *argv[])
{
    TCLAP::CmdLine               cmd("EPP PROGRAM. COMPUTE BAND STRUCTURE ON A BZ MESH.", ' ', "1.0");
    TCLAP::ValueArg<std::string> arg_mesh_file("f", "meshbandfile", "File with BZ mesh and bands energy.", true, "bz.msh", "string");
    TCLAP::ValueArg<std::string> arg_material("m", "material", "Symbol of the material to use (Si, Ge, GaAs, ...)", true, "Si", "string");
    TCLAP::ValueArg<int>         arg_nb_energies("e", "nenergy", "Number of energies to compute", false, 250, "int");
    TCLAP::ValueArg<int>         arg_nb_conduction_bands("c", "ncbands", "Number of conduction bands to consider", false, -1, "int");
    TCLAP::ValueArg<int>         arg_nb_valence_bands("v", "nvbands", "Number of valence bands to consider", false, -1, "int");
    TCLAP::ValueArg<int>         arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    TCLAP::ValueArg<double>      arg_temperature("T", "temperature", "Temperature in Kelvin.", false, 300.0, "double");
    TCLAP::ValueArg<double>      arg_energy_range("E", "energy_window", "Energy window around the band gap to consider (in eV).", false, 0.3, "double");
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

    EmpiricalPseudopotential::Materials materials;
    const std::string                   file_material_parameters = std::string(PROJECT_SRC_DIR) + "/parameter_files/materials-local.yaml";
    materials.load_material_parameters(file_material_parameters);

    Options my_options;
    my_options.materialName = arg_material.getValue();
    my_options.nrLevels     = arg_nb_conduction_bands.getValue() + arg_nb_valence_bands.getValue();
    my_options.nrThreads    = arg_nb_threads.getValue();
    int number_energies     = arg_nb_energies.getValue();

    EmpiricalPseudopotential::Material current_material = materials.materials.at(arg_material.getValue());

    const std::string mesh_band_input_file = arg_mesh_file.getValue();
    bz_mesh::ElectronPhonon   ElectronPhonon{current_material};
    ElectronPhonon.set_nb_threads(my_options.nrThreads);

    // const std::string phonon_file = std::string(PROJECT_SRC_DIR) + "/parameter_files/phonon_michaillat.yaml";
    const std::string phonon_file = std::string(PROJECT_SRC_DIR) + "/parameter_files/phonon_kamakura.yaml";
    
    
    ElectronPhonon.read_mesh_geometry_from_msh_file(mesh_band_input_file);
    ElectronPhonon.read_mesh_bands_from_msh_file(mesh_band_input_file, my_options.nrLevels);
    constexpr double energy_step_dos = 0.005; // eV
    ElectronPhonon.precompute_dos_tetra(energy_step_dos);

    unsigned int nb_bands = ElectronPhonon.get_number_bands();
    std::cout << "Number of bands: " << nb_bands << std::endl;
    if (my_options.nrLevels > nb_bands) {
        std::cout << "Number of bands requested is greater than the number of bands in the mesh file. Resetting to " << nb_bands << std::endl;
        my_options.nrLevels = nb_bands;
    }
    
    
    ElectronPhonon.load_phonon_parameters(phonon_file);
    bool irreducible_wedge_only = plot_with_wedge.getValue();
    bool populate_nk_npkp = plot_with_knkpnp.getValue();

    // const double max_energy  = 6.0;  // eV
    const double max_energy  = arg_energy_range.getValue();  // eV
    const double temperature = arg_temperature.getValue();
    ElectronPhonon.set_temperature(temperature);
    std::cout << "Max energy: " << max_energy << " eV" << std::endl;
    const double energy_step = 0.05; // eV
    ElectronPhonon.compute_electron_phonon_rates_over_mesh(max_energy, irreducible_wedge_only, populate_nk_npkp);
    ElectronPhonon.export_rate_values("rates_all.csv");

    ElectronPhonon.compute_plot_electron_phonon_rates_vs_energy_over_mesh(my_options.nrLevels, max_energy, energy_step, "rates_vs_energy.csv", irreducible_wedge_only);

    // ElectronPhonon.add_electron_phonon_rates_to_mesh(mesh_band_input_file, "rates.msh");


    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds\n\n\n" << std::endl;

    if (plot_with_python.getValue()) {
        std::string command = "python3 " + std::string(PROJECT_SRC_DIR) + "/python/plots/plot_phonon_rate.py -f rates_vs_energy.csv";
        std::cout << "Running command: " << command << std::endl;
        std::system(command.c_str());
    }

    return 0;

}
