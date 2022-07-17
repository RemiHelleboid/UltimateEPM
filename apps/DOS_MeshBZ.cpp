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

inline void export_multiple_vector_to_csv(const std::string &                     filename,
                                          const std::vector<std::string> &        header_columns,
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
    TCLAP::ValueArg<int>         arg_nb_bands("b", "nbands", "Number of bands to consider", false, 12, "int");
    TCLAP::ValueArg<int>         arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    cmd.add(arg_mesh_file);
    cmd.add(arg_material);
    cmd.add(arg_nb_bands);
    cmd.add(arg_nb_threads);

    cmd.parse(argc, argv);

    EmpiricalPseudopotential::Materials materials;

    Options my_options;
    my_options.materialName = arg_material.getValue();
    my_options.nrLevels     = arg_nb_bands.getValue();
    my_options.nrThreads    = arg_nb_threads.getValue();
    my_options.print_options();

    bz_mesh::MeshBZ my_bz_mesh{};
    my_bz_mesh.read_mesh_geometry_from_msh_file(arg_mesh_file.getValue());
    my_bz_mesh.read_mesh_bands_from_msh_file(arg_mesh_file.getValue());
    my_bz_mesh.compute_iso_surface(1.10, 0);

    std::cout << "Mesh volume: " << my_bz_mesh.compute_mesh_volume() << std::endl;

    std::vector<std::vector<double>> list_list_dos{};
    std::vector<std::string>         list_header = {};

    std::size_t number_bands = my_bz_mesh.get_number_bands();

    std::cout << "Compute DOS on " << number_bands << " bands.\n";

    for (int band_index = 0; band_index < number_bands; ++band_index) {
        std::vector<std::vector<double>> lists_energies_dos = my_bz_mesh.compute_dos_band_at_band_auto(band_index, 100);
        list_list_dos.push_back(lists_energies_dos[0]);
        list_list_dos.push_back(lists_energies_dos[1]);
        list_header.push_back("energy_band_" + std::to_string(band_index));
        list_header.push_back("dos_band_" + std::to_string(band_index));
    }

    export_multiple_vector_to_csv("DOS_BANDS.csv", list_header, list_list_dos);


    return 0;
}