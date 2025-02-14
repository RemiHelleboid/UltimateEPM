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
    TCLAP::ValueArg<int>         arg_nb_bands("b", "nbands", "Number of bands to consider", false, 12, "int");
    TCLAP::ValueArg<int>         arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    TCLAP::SwitchArg plot_with_python("P", "plot", "Call a python script after the computation to plot the band structure.", false);
    cmd.add(plot_with_python);
    cmd.add(arg_mesh_file);
    cmd.add(arg_material);
    cmd.add(arg_nb_bands);
    cmd.add(arg_nb_energies);
    cmd.add(arg_nb_threads);

    cmd.parse(argc, argv);

    auto start = std::chrono::high_resolution_clock::now();

    EmpiricalPseudopotential::Materials materials;
    const std::string                   file_material_parameters = std::string(CMAKE_SOURCE_DIR) + "/parameter_files/materials.yaml";
    materials.load_material_parameters(file_material_parameters);

    Options my_options;
    my_options.materialName = arg_material.getValue();
    my_options.nrLevels     = arg_nb_bands.getValue();
    my_options.nrThreads    = arg_nb_threads.getValue();
    int number_energies     = arg_nb_energies.getValue();

    EmpiricalPseudopotential::Material current_material = materials.materials.at(arg_material.getValue());

    const std::string mesh_band_input_file = arg_mesh_file.getValue();
    bz_mesh::ElectronPhonon   ElectronPhonon{current_material};
    const std::string phonon_file = std::string(CMAKE_SOURCE_DIR) + "/parameter_files/phonon_michaillat.yaml";
    ElectronPhonon.load_phonon_parameters(phonon_file);


    ElectronPhonon.read_mesh_geometry_from_msh_file(mesh_band_input_file);
    ElectronPhonon.read_mesh_bands_from_msh_file(mesh_band_input_file);

    unsigned int nb_bands = ElectronPhonon.get_number_bands();
    std::cout << "Number of bands: " << nb_bands << std::endl;
    if (my_options.nrLevels > nb_bands) {
        std::cout << "Number of bands requested is greater than the number of bands in the mesh file. Resetting to " << nb_bands << std::endl;
        my_options.nrLevels = nb_bands;
    }


    ElectronPhonon.compute_electron_phonon_rates_over_mesh();

    ElectronPhonon.export_rate_values("rates_all.csv");

    ElectronPhonon.compute_plot_electron_phonon_rates_vs_energy_over_mesh(my_options.nrLevels, 10.0, 0.01, "rates_vs_energy.csv");

    ElectronPhonon.add_electron_phonon_rates_to_mesh(mesh_band_input_file, "rates.msh");


    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " seconds" << std::endl;

    return 0;

}
