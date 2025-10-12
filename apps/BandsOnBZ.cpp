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
#include <fstream>
#include <iostream>

#include "BandStructure.h"
#include "Material.h"
#include "Options.h"
#include "bz_meshfile.hpp"

int main(int argc, char* argv[]) {
    TCLAP::CmdLine               cmd("EPP PROGRAM. COMPUTE BAND STRUCTURE ON A BZ MESH.", ' ', "1.0");
    TCLAP::ValueArg<std::string> arg_mesh_file("f", "meshfile", "Name to print", true, "bz.msh", "string");
    TCLAP::ValueArg<std::string> arg_material("m", "material", "Symbol of the material to use (Si, Ge, GaAs, ...)", true, "Si", "string");
    TCLAP::ValueArg<std::string> arg_outfile("o", "outfile", "Name of the output file", false, "", "string");
    TCLAP::ValueArg<int>         arg_nb_valence_bands("v", "nvbands", "Number of valence bands to export", false, 4, "int");
    TCLAP::ValueArg<int>         arg_nb_conduction_bands("c", "ncbands", "Number of conduction bands to export", false, 12, "int");
    TCLAP::ValueArg<int>         arg_nearest_neighbors("n",
                                               "nearestNeighbors",
                                               "number of nearest neiborgs to consider for the EPP calculation.",
                                               false,
                                               10,
                                               "int");
    TCLAP::SwitchArg arg_enable_nonlocal_correction("C", "nonlocal-correction", "Enable the non-local-correction for the EPM model", false);
    TCLAP::SwitchArg arg_enable_soc("s", "soc", "Enable the spin-orbit coupling for the EPM model", false);
    TCLAP::SwitchArg arg_cond_band_zero("z", "MinCondZero", "Shift the conduction band minimum to 0 eV", false);
    TCLAP::ValueArg<int> arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    cmd.add(arg_mesh_file);
    cmd.add(arg_material);
    cmd.add(arg_nb_valence_bands);
    cmd.add(arg_nb_conduction_bands);
    cmd.add(arg_outfile);
    cmd.add(arg_nearest_neighbors);
    cmd.add(arg_nb_threads);
    cmd.add(arg_enable_nonlocal_correction);
    cmd.add(arg_enable_soc);
    cmd.add(arg_cond_band_zero);

    cmd.parse(argc, argv);

    uepm::pseudopotential::Materials materials;
    const std::string                   file_material_parameters = std::string(PROJECT_SRC_DIR) + "/parameter_files/materials-local.yaml";
    materials.load_material_parameters(file_material_parameters);
    bool enable_nonlocal_correction = arg_enable_nonlocal_correction.isSet();
    bool enable_soc                 = arg_enable_soc.isSet();

    Options my_options;
    my_options.materialName = arg_material.getValue();

    constexpr int max_valence_bands = 8;
    my_options.nrLevels             = arg_nb_valence_bands.getValue() + arg_nb_conduction_bands.getValue() +
                          max_valence_bands;  // add extra bands for the calculations (won't be exported)
    my_options.nearestNeighbors = arg_nearest_neighbors.getValue();
    my_options.nrThreads        = arg_nb_threads.getValue();
    my_options.print_options();

    uepm::pseudopotential::Material mat = materials.materials.at(my_options.materialName);

    // bz_mesh my_mesh("mesh.msh");
    const std::string mesh_filename = arg_mesh_file.getValue();
    bz_mesh_points    my_mesh(mesh_filename);
    my_mesh.read_mesh();
    std::vector<Vector3D<double>>& mesh_kpoints = my_mesh.get_kpoints();

    auto start = std::chrono::high_resolution_clock::now();

    uepm::pseudopotential::BandStructure my_bandstructure;
    my_bandstructure
        .Initialize(mat, my_options.nrLevels, mesh_kpoints, my_options.nearestNeighbors, enable_nonlocal_correction, arg_enable_soc);

    if (my_options.nrThreads > 1) {
        my_bandstructure.Compute_parallel(my_options.nrThreads);
    } else {
        my_bandstructure.Compute();
    }
    my_bandstructure.AdjustValues(arg_cond_band_zero.getValue());

    auto end              = std::chrono::high_resolution_clock::now();
    auto total_time_count = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    // std::ofstream f_time("OpenMP_Bands_times.csv", std::ios::app);
    // f_time << my_options.nrThreads << "," << total_time_count / double(1000) << std::endl;
    // f_time.close();

    std::cout << "Total computation time: " << total_time_count / double(1000) << std::endl;

    std::filesystem::path in_path(mesh_filename);
    std::string           out_file_bands = in_path.stem().replace_extension("").string() + "_" + my_bandstructure.path_band_filename();

    if (arg_outfile.isSet()) {
        out_file_bands = arg_outfile.getValue();
    }

    my_mesh.add_all_bands_on_mesh(out_file_bands, my_bandstructure, arg_nb_valence_bands.getValue(), arg_nb_conduction_bands.getValue());

    return 0;
}