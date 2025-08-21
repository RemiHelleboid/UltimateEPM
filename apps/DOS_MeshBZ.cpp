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
#include "bz_mesh.hpp"
#include "bz_meshfile.hpp"

inline void export_multiple_vector_to_csv(const std::string                      &filename,
                                          const std::vector<std::string>         &header_columns,
                                          const std::vector<std::vector<double>> &value_vector_of_vector) {
    if (value_vector_of_vector.empty()) {
        return;
    }
    const std::size_t reference_vector_size = value_vector_of_vector[0].size();
    for (auto &&vector : value_vector_of_vector) {
        if (vector.size() != reference_vector_size) {
            std::cout << "ERROR WHEN EXPORTING VECTORS TO : " << filename << ", mismatch between vector sizes : " << reference_vector_size
                      << " != " << vector.size() << std::endl;
            return;
        }
    }
    const std::string DumbColumnName = "DumbColumn\n";
    const double      dumb_value     = 0.0;
    std::ofstream     csv_file(filename);
    for (auto &&col_name : header_columns) {
        csv_file << col_name << ",";
    }
    csv_file << DumbColumnName;
    for (std::size_t index_value = 0; index_value < value_vector_of_vector[0].size(); ++index_value) {
        for (std::size_t index_vector = 0; index_vector < value_vector_of_vector.size(); ++index_vector) {
            csv_file << value_vector_of_vector[index_vector][index_value] << ",";
        }
        csv_file << dumb_value << "\n";
    }
    csv_file.close();
}

int main(int argc, char *argv[]) {
    TCLAP::CmdLine               cmd("EPP PROGRAM. COMPUTE BAND STRUCTURE ON A BZ MESH.", ' ', "1.0");
    TCLAP::ValueArg<std::string> arg_mesh_file("f", "meshbandfile", "File with BZ mesh and bands energy.", true, "bz.msh", "string");
    TCLAP::ValueArg<std::string> arg_material("m", "material", "Symbol of the material to use (Si, Ge, GaAs, ...)", true, "Si", "string");
    TCLAP::ValueArg<int>         arg_nb_energies("e", "nenergy", "Number of energies to compute", false, 250, "int");
    TCLAP::ValueArg<int>         arg_nb_bands("b", "nbands", "Number of bands to consider", false, -1, "int");
    TCLAP::ValueArg<int>         arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    TCLAP::SwitchArg plot_with_python("P", "plot", "Call a python script after the computation to plot the band structure.", false);
    cmd.add(plot_with_python);
    cmd.add(arg_mesh_file);
    cmd.add(arg_material);
    cmd.add(arg_nb_bands);
    cmd.add(arg_nb_energies);
    cmd.add(arg_nb_threads);

    cmd.parse(argc, argv);

    EmpiricalPseudopotential::Materials materials;
    const std::string                   file_material_parameters = std::string(CMAKE_SOURCE_DIR) + "/parameter_files/materials-local.yaml";
    materials.load_material_parameters(file_material_parameters);

    Options my_options;
    my_options.materialName = arg_material.getValue();
    my_options.nrLevels     = arg_nb_bands.getValue();
    my_options.nrThreads    = arg_nb_threads.getValue();
    my_options.print_options();
    int  nb_bands_to_use = arg_nb_bands.getValue();
    auto start           = std::chrono::high_resolution_clock::now();

    EmpiricalPseudopotential::Material current_material = materials.materials.at(arg_material.getValue());

    const std::string mesh_band_input_file = arg_mesh_file.getValue();

    bz_mesh::MeshBZ my_bz_mesh{current_material};
    my_bz_mesh.read_mesh_geometry_from_msh_file(mesh_band_input_file);
    my_bz_mesh.read_mesh_bands_from_msh_file(mesh_band_input_file, nb_bands_to_use);

    std::cout << "Mesh volume: " << my_bz_mesh.compute_mesh_volume() << std::endl;
    double Vcell = std::pow(current_material.get_lattice_constant_meter(), 3) / 4.0;
    std::cout << "Vcell: " << Vcell << std::endl;
    const double VBZ_theory = std::pow(2.0 * M_PI, 3) / Vcell;   // m^-3
    const double VBZ_mesh   = my_bz_mesh.compute_mesh_volume();  // sum |tetra.volume| in your k units
    std::cout << "VBZ_mesh / VBZ_theory = " << (VBZ_mesh / VBZ_theory) << "\n";
    my_bz_mesh.set_reduce_bz_factor(VBZ_theory/VBZ_mesh);

    std::vector<std::vector<double>> list_list_dos{};
    std::vector<std::string>         list_header = {};

    std::cout << "Compute DOS on " << nb_bands_to_use << " bands.\n";

    for (int band_index = 0; band_index < nb_bands_to_use; ++band_index) {
        std::vector<std::vector<double>> lists_energies_dos =
            my_bz_mesh.compute_dos_band_at_band_auto(band_index, arg_nb_energies.getValue(), my_options.nrThreads);
        list_list_dos.push_back(lists_energies_dos[0]);
        list_list_dos.push_back(lists_energies_dos[1]);
        list_header.push_back("energy_band_" + std::to_string(band_index));
        list_header.push_back("dos_band_" + std::to_string(band_index));
        bz_mesh::Tetra::print_stat_iso_computing();
        bz_mesh::Tetra::reset_stat_iso_computing();
    }

    auto end              = std::chrono::high_resolution_clock::now();
    auto total_time_count = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::ofstream f_time("OpenMP_DOS_times.csv", std::ios::app);
    f_time << my_options.nrThreads << "," << total_time_count / double(1000) << std::endl;
    f_time.close();

    std::cout << "Total computation time: " << total_time_count / double(1000) << std::endl;

    std::filesystem::path in_path(mesh_band_input_file);
    std::string           out_file_bands = "DOS_" + in_path.stem().replace_extension("").string();

    export_multiple_vector_to_csv(out_file_bands + ".csv", list_header, list_list_dos);

    const std::string python_plot_dos  = std::string(CMAKE_SOURCE_DIR) + "/python/plots/plot_density_of_states.py";
    bool              call_python_plot = plot_with_python.isSet();
    if (call_python_plot) {
        std::string python_call = "python3 " + python_plot_dos + " --file " + out_file_bands + ".csv";
        std::cout << "Executing: " << python_call << std::endl;
        int succes_plot = system(python_call.c_str());
        if (succes_plot != 0) {
            std::cout << "Error while calling python script to plot DOS.\n";
        }
    }

    return 0;
}