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
#include "dielectric_mesh.hpp"
#include "impact_ionization.hpp"

int main(int argc, char *argv[]) {
    TCLAP::CmdLine               cmd("EPP PROGRAM. COMPUTE BAND STRUCTURE ON A BZ MESH.", ' ', "1.0");
    TCLAP::ValueArg<std::string> arg_mesh_file("f", "meshbandfile", "File with BZ mesh and bands energy.", true, "bz.msh", "string");
    TCLAP::ValueArg<std::string> arg_dielectric_file("y", "dielectric", "File (.msh) with the dielectric function.", true, "", "string");
    TCLAP::ValueArg<std::string> arg_material("m", "material", "Symbol of the material to use (Si, Ge, GaAs, ...)", true, "Si", "string");
    TCLAP::ValueArg<int>         arg_nb_energies("e", "nenergy", "Number of energies to compute", false, 250, "int");
    TCLAP::ValueArg<int>         arg_nb_bands("b", "nbands", "Number of bands to consider", false, 16, "int");
    TCLAP::ValueArg<int>         arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    TCLAP::ValueArg<double>      arg_radius_BZ("R",
                                          "radiusBZ",
                                          "Max norm of G vector for which the BZ center in G is taken into account",
                                          false,
                                          0,
                                          "double");
    TCLAP::SwitchArg plot_with_python("P", "plot", "Call a python script after the computation to plot the band structure.", false);
    cmd.add(plot_with_python);
    cmd.add(arg_mesh_file);
    cmd.add(arg_dielectric_file);
    cmd.add(arg_material);
    cmd.add(arg_nb_bands);
    cmd.add(arg_nb_energies);
    cmd.add(arg_nb_threads);
    cmd.add(arg_radius_BZ);

    cmd.parse(argc, argv);

    bool                                nonlocal_epm = false;
    bool                                enable_soc   = false;
    uepm::pseudopotential::Materials materials;
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
    int  nb_threads      = arg_nb_threads.getValue();
    auto start           = std::chrono::high_resolution_clock::now();

    uepm::pseudopotential::Material current_material = materials.materials.at(arg_material.getValue());

    uepm::mesh_bz::ImpactIonization my_impact_ionization(current_material, arg_mesh_file.getValue());
    my_impact_ionization.read_dielectric_file(arg_dielectric_file.getValue());
    my_impact_ionization.interp_test_dielectric_function("test_dielectric_function.csv");

    my_impact_ionization.set_max_radius_G0_BZ(arg_radius_BZ.getValue());
    my_impact_ionization.compute_eigenstates(nb_threads);
    int idxn1 = 0;
    int idxk1 = 0;
    my_impact_ionization.compute_impact_ionization_rate(idxn1, idxk1);

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s" << std::endl;

    return 0;
}