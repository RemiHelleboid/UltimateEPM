/**
 * @file fullstates.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-12-20
 *
 * @copyright Copyright (c) 2022
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

int main(int argc, char *argv[]) {
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

    bool                                nonlocal_epm = false;
    bool                                enable_soc   = false;
    EmpiricalPseudopotential::Materials materials;
    std::string                         file_material_parameters = std::string(PROJECT_SRC_DIR) + "/parameter_files/materials-cohen.yaml";
    if (nonlocal_epm) {
        file_material_parameters = std::string(PROJECT_SRC_DIR) + "/parameter_files/materials.yaml";
    }
    std::cout << "Loading material parameters from " << file_material_parameters << std::endl;
    materials.load_material_parameters(file_material_parameters);

    Options my_options;
    my_options.materialName = arg_material.getValue();
    my_options.nrLevels     = arg_nb_bands.getValue();
    my_options.nrThreads    = arg_nb_threads.getValue();
    my_options.print_options();
    int  nb_bands_to_use = arg_nb_bands.getValue();
    auto start           = std::chrono::high_resolution_clock::now();

    EmpiricalPseudopotential::Material current_material = materials.materials.at(arg_material.getValue());

    bz_mesh::BZ_States my_bz_mesh(current_material);
    my_bz_mesh.set_nb_bands(nb_bands_to_use);
    EmpiricalPseudopotential::BandStructure band_structure{};

    int nb_nearest_neighbors = 10;
    band_structure.Initialize(current_material, nb_bands_to_use, {}, nb_nearest_neighbors, nonlocal_epm, enable_soc);
    auto basis = band_structure.get_basis_vectors();
    my_bz_mesh.set_basis_vectors(basis);

    my_bz_mesh.read_mesh_geometry_from_msh_file(arg_mesh_file.getValue());
    std::cout << "Mesh volume: " << my_bz_mesh.compute_mesh_volume() << std::endl;

    my_bz_mesh.compute_eigenstates(my_options.nrThreads);

    Vector3D<double> q_shift = Vector3D<double>{3e-14, 0.0, 0.0};
    my_bz_mesh.compute_shifted_eigenstates(q_shift, my_options.nrThreads);
    std::cout << "\n\n" << std::endl;
    // my_bz_mesh.export_full_eigenstates();

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s" << std::endl;

    std::vector<double> list_energy;
    double              min_energy  = 0.0;
    double              max_energy  = 20.0;
    double              energy_step = 0.01;
    for (double energy = min_energy; energy <= max_energy + energy_step; energy += energy_step) {
        list_energy.push_back(energy);
    }
    double eta_smearing = 0.05;
    std::cout << "Number of energies to compute: " << list_energy.size() << std::endl;
    auto start2 = std::chrono::high_resolution_clock::now();
    my_bz_mesh.compute_dielectric_function(list_energy, eta_smearing, my_options.nrThreads);
    my_bz_mesh.export_dielectric_function("./TEST_DIELECTRIC_FUNCTION_");
    auto                          end2             = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds2 = end2 - start2;
    std::cout << "Time Dielectric Function : " << elapsed_seconds2.count() << "s" << std::endl;

    std::string python_script =
        "import matplotlib.pyplot as plt\n"
        "import numpy as np\n"
        "energy, eps = np.loadtxt('./TEST_DIELECTRIC_FUNCTION__dielectric_function.csv', delimiter=',', unpack=True, skiprows=1)\n"
        "fig, ax = plt.subplots(figsize=(8, 6))\n"
        "ax.plot(energy, eps, label='Dielectric Function')\n"
        "ax.set_xlabel('Energy (eV)')\n"
        "ax.set_ylabel('Dielectric Function (ε)')\n"
        "ax.set_title('Dielectric Function vs Energy')\n"
        "ax.legend()\n"
        "ax.grid()\n"
        "plt.savefig('dielectric_function.png')\n"
        "plt.show()\n";
    std::cout << "Run python script to plot dielectric function." << std::endl;
    std::string command = "python -c \"" + python_script + "\"";
    system(command.c_str());
    return 0;
}