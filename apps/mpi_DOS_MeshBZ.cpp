/**
 * @file mpi_DOS_MeshBZ.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-08-03
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <limits.h>
#include <mpi.h>  // must have a system with an MPI library
#include <stdint.h>
#include <stdio.h>   //printf
#include <stdlib.h>  //malloc
#include <tclap/CmdLine.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

#include "BandStructure.h"
#include "Material.h"
#include "Options.h"
#include "bz_mesh.hpp"
#include "bz_meshfile.hpp"

#if SIZE_MAX == UCHAR_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
#error "what is happening here?"
#endif

#define MASTER 0

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
    std::cout << "EXPORTING VECTORS TO : " << filename << std::endl;
    std::cout << "Number column: " << value_vector_of_vector.size() << std::endl;
    std::cout << "Number rows: " << value_vector_of_vector[0].size() << std::endl;
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

template <typename T>
std::vector<T> linspace(T x_min, T x_max, std::size_t number_points) {
    std::vector<T> list_x;
    list_x.resize(number_points);
    if (number_points == 0) {
        return list_x;
    }
    if (number_points == 1) {
        list_x[0] = x_min;
        return list_x;
    }
    double dx = (x_max - x_min) / (number_points - 1);
    for (std::size_t index_value = 0; index_value < number_points; ++index_value) {
        list_x[index_value] = x_min + dx * index_value;
    }
    return list_x;
}

int main(int argc, char *argv[]) {
    // Initialize the MPI environment.
    // MPI_Status status;

    int number_processes;
    int process_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &number_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

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

    uepm::pseudopotential::Materials materials;
    const std::string                   file_material_parameters = std::string(PROJECT_SRC_DIR) + "/parameter_files/materials.yaml";
    materials.load_material_parameters(file_material_parameters);

    Options my_options;
    my_options.materialName = arg_material.getValue();
    my_options.nrLevels     = arg_nb_bands.getValue();
    my_options.nrThreads    = arg_nb_threads.getValue();
    int number_energies     = arg_nb_energies.getValue();

    uepm::pseudopotential::Material current_material = materials.materials.at(arg_material.getValue());

    const std::string mesh_band_input_file = arg_mesh_file.getValue();
    uepm::mesh_bz::MeshBZ   my_bz_mesh{current_material};
    my_bz_mesh.read_mesh_geometry_from_msh_file(mesh_band_input_file);
    my_bz_mesh.read_mesh_bands_from_msh_file(mesh_band_input_file);

    std::size_t         number_bands = my_bz_mesh.get_number_bands();
    std::vector<int>    list_bands;
    std::vector<double> list_energies;
    std::size_t         total_number_dos = 0;

    if (process_rank == MASTER) {
        std::cout << "Number of bands : " << number_bands << std::endl;
        for (std::size_t index_band = 0; index_band < number_bands; index_band++) {
            double              margin_eV      = 0.01;
            auto                min_max_energy = my_bz_mesh.get_min_max_energy_at_band(index_band);
            std::vector<double> energies =
                linspace<double>(min_max_energy.first - margin_eV, min_max_energy.second + margin_eV, number_energies);
            std::cout << "Band " << index_band << " : " << min_max_energy.first << " " << min_max_energy.second << std::endl;
            std::vector<int> bands = std::vector<int>(number_energies, index_band);
            list_bands.insert(list_bands.end(), bands.begin(), bands.end());
            list_energies.insert(list_energies.end(), energies.begin(), energies.end());
        }
        total_number_dos = list_energies.size();

        std::cout << "Size bands : " << list_bands.size() << std::endl;
        std::cout << "Size energies : " << list_energies.size() << std::endl;
    }
    MPI_Bcast(&total_number_dos, 1, my_MPI_SIZE_T, MASTER, MPI_COMM_WORLD);

    double t_start = MPI_Wtime();

    if (process_rank == MASTER) {
        std::cout << "Total number of energies to compute : " << total_number_dos << std::endl;
    }

    // Define the number of elements each process will handle.
    int count     = (total_number_dos / number_processes);
    int remainder = (total_number_dos % number_processes);

    std::vector<int> counts_dos_per_process(number_processes);
    std::vector<int> displacements_dos_per_process(number_processes);
    for (int i = 0; i < number_processes - 1; i++) {
        counts_dos_per_process[i]        = count;
        displacements_dos_per_process[i] = i * count;
    }
    counts_dos_per_process.back()        = (count + remainder);
    displacements_dos_per_process.back() = ((number_processes - 1) * count);

    std::vector<int>    chunk_band_indices(counts_dos_per_process[process_rank]);
    std::vector<double> chunk_energies(counts_dos_per_process[process_rank]);

    MPI_Scatterv(list_bands.data(),
                 counts_dos_per_process.data(),
                 displacements_dos_per_process.data(),
                 MPI_INT,
                 chunk_band_indices.data(),
                 counts_dos_per_process[process_rank],
                 MPI_INT,
                 MASTER,
                 MPI_COMM_WORLD);
    MPI_Scatterv(list_energies.data(),
                 counts_dos_per_process.data(),
                 displacements_dos_per_process.data(),
                 MPI_DOUBLE,
                 chunk_energies.data(),
                 counts_dos_per_process[process_rank],
                 MPI_DOUBLE,
                 MASTER,
                 MPI_COMM_WORLD);

    std::vector<double> dos_values(counts_dos_per_process[process_rank]);
    std::cout << "Process " << process_rank << " will compute " << counts_dos_per_process[process_rank] << " DOS values." << std::endl;

    for (int index_dos = 0; index_dos < counts_dos_per_process[process_rank]; index_dos++) {
        dos_values[index_dos] = my_bz_mesh.compute_dos_at_energy_and_band(chunk_energies[index_dos], chunk_band_indices[index_dos]);
    }

    std::vector<double> all_dos_values;
    if (process_rank == MASTER) {
        all_dos_values.resize(total_number_dos);
    }

    MPI_Gatherv(dos_values.data(),
                counts_dos_per_process[process_rank],
                MPI_DOUBLE,
                all_dos_values.data(),
                counts_dos_per_process.data(),
                displacements_dos_per_process.data(),
                MPI_DOUBLE,
                MASTER,
                MPI_COMM_WORLD);

    double t_end = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    if (process_rank == MASTER) {
        std::ofstream time_dos("time_dos.csv", std::ios::app);
        time_dos << number_processes << "," << t_end - t_start << std::endl;
        time_dos.close();
        std::cout << "------------------------------------------------------------------\n";
        std::cout << "Total time : " << t_end - t_start << std::endl;
        std::cout << "------------------------------------------------------------------\n";

        std::cout << "Exporting DOS values to file." << std::endl;
        std::vector<std::vector<double>> dos_at_bands(number_bands * 2);
        for (std::size_t index_value_dos = 0; index_value_dos < total_number_dos; index_value_dos++) {
            int    band_index = list_bands[index_value_dos];
            double energy     = list_energies[index_value_dos];
            double dos        = all_dos_values[index_value_dos];
            dos_at_bands[2 * band_index].push_back(energy);
            dos_at_bands[2 * band_index + 1].push_back(dos);
        }
        std::vector<std::string> list_header;
        for (std::size_t index_band = 0; index_band < number_bands; index_band++) {
            list_header.push_back("energy_band_" + std::to_string(index_band));
            list_header.push_back("dos_band_" + std::to_string(index_band));
        }
        std::filesystem::path in_path(mesh_band_input_file);
        std::string           out_file_bands = "DOS_MPI_" + in_path.stem().replace_extension("").string();
        std::cout << "Exporting DOS values to file: " << out_file_bands << std::endl;
        export_multiple_vector_to_csv(out_file_bands + ".csv", list_header, dos_at_bands);

        const std::string python_plot_dos  = std::string(PROJECT_SRC_DIR) + "/python/plots/plot_density_of_states.py";
        bool              call_python_plot = plot_with_python.isSet();
        if (call_python_plot) {
            std::string python_call = "nohup python3 " + python_plot_dos + " --file " + out_file_bands + ".csv &";
            std::cout << "Executing: " << python_call << std::endl;
            int succes_plot = system(python_call.c_str());
            if (succes_plot != 0) {
                std::cout << "Error while executing python script to plot DOS." << std::endl;
            }
        }
    }

    MPI_Finalize();
    return 0;
}
