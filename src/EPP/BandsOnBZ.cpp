/**
 * @file compute_bands_on_bz_mesh.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-07-07
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <tclap/CmdLine.h>

#include <filesystem>

#include "BandStructure.h"
#include "Material.h"
#include "Options.h"
#include "bz_meshfile.hpp"

int main(int argc, char* argv[]) {
    TCLAP::CmdLine               cmd("EPP PROGRAM. COMPUTE BAND STRUCTURE ON A BZ MESH.", ' ', "1.0");
    TCLAP::ValueArg<std::string> arg_mesh_file("f", "meshfile", "Name to print", true, "bz.msh", "string");
    TCLAP::ValueArg<std::string> arg_material("m", "material", "Symbol of the material to use (Si, Ge, GaAs, ...)", true, "Si", "string");
    TCLAP::ValueArg<int>         arg_nb_bands("b", "nbands", "Number of bands to compute", false, 12, "int");
    TCLAP::ValueArg<int>         arg_nearest_neighbors("n",
                                               "nearestNeighbors",
                                               "number of nearest neiborgs to consider for the EPP calculation.",
                                               false,
                                               10,
                                               "int");
    TCLAP::ValueArg<int>         arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    cmd.add(arg_mesh_file);
    cmd.add(arg_material);
    cmd.add(arg_nb_bands);
    cmd.add(arg_nearest_neighbors);
    cmd.add(arg_nb_threads);

    cmd.parse(argc, argv);

    EmpiricalPseudopotential::Materials materials;

    Options my_options;
    my_options.materialName     = arg_material.getValue();
    my_options.nrLevels         = arg_nb_bands.getValue();
    my_options.nearestNeighbors = arg_nearest_neighbors.getValue();
    my_options.nrThreads        = arg_nb_threads.getValue();
    my_options.print_options();

    EmpiricalPseudopotential::Material mat = materials.materials.at(my_options.materialName);

    // bz_mesh my_mesh("mesh.msh");
    const std::string mesh_filename = arg_mesh_file.getValue();
    bz_mesh           my_mesh(mesh_filename);
    my_mesh.read_mesh();
    std::vector<Vector3D<double>>& mesh_kpoints = my_mesh.get_kpoints();

    EmpiricalPseudopotential::BandStructure my_bandstructure;
    my_bandstructure.Initialize(mat, my_options.nrLevels, mesh_kpoints, my_options.nearestNeighbors);

    if (my_options.nrThreads > 1) {
        my_bandstructure.Compute_parralel(my_options.nrThreads);
    } else {
        my_bandstructure.Compute();
    }
    my_bandstructure.AdjustValues();

    unsigned int        band_index = 0;
    std::vector<double> band_0     = my_bandstructure.get_band(band_index);
    // export band in file
    // my_bandstructure.export_result_in_file_with_kpoints("BZ_BANDS_SI.csv");

    std::string out_file_bands = my_bandstructure.path_band_filename();

    my_mesh.add_all_bands_on_mesh(out_file_bands + "_all_bands.msh", my_bandstructure);

    return 0;
}