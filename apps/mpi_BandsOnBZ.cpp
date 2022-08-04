/**
 * @file mpi_BandsOnBZ.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief MPI version of the BandsOnBZ program.
 * @version 0.1
 * @date 2022-08-03
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <mpi.h>     // must have a system with an MPI library
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
#include "bz_meshfile.hpp"

#define MASTER 0

// Define the structure to be serialized.
// It is a simple 3D vector of double.
// We define a unique method, which computes the norm of the vector.
typedef struct vector_k {
    double m_kx;
    double m_ky;
    double m_kz;

    vector_k() : m_kx(0), m_ky(0), m_kz(0) {}
    vector_k(double kx, double ky, double kz) : m_kx(kx), m_ky(ky), m_kz(kz) {}
    void set_k(double kx, double ky, double kz) {
        m_kx = kx;
        m_ky = ky;
        m_kz = kz;
    }
    double           norm() const { return std::sqrt(m_kx * m_kx + m_ky * m_ky + m_kz * m_kz); }
    Vector3D<double> to_Vector3D() const { return Vector3D<double>(m_kx, m_ky, m_kz); }
} vector_k;

int main(int argc, char* argv[]) {
    // Initialize the MPI environment.
    MPI_Status status;

    int number_processes;
    int process_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &number_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);


    TCLAP::CmdLine               cmd("EPP PROGRAM. COMPUTE BAND STRUCTURE ON A BZ MESH.", ' ', "1.0");
    TCLAP::ValueArg<std::string> arg_mesh_file("f", "meshfile", "Name to print", true, "bz.msh", "string");
    TCLAP::ValueArg<std::string> arg_material("m", "material", "Symbol of the material to use (Si, Ge, GaAs, ...)", true, "Si", "string");
    TCLAP::ValueArg<std::string> arg_outfile("o", "outfile", "Name of the output file", false, "", "string");
    TCLAP::ValueArg<int>         arg_nb_bands("b", "nbands", "Number of bands to compute", false, 12, "int");
    TCLAP::ValueArg<int>         arg_nearest_neighbors("n",
                                               "nearestNeighbors",
                                               "number of nearest neiborgs to consider for the EPP calculation.",
                                               false,
                                               10,
                                               "int");
    TCLAP::SwitchArg arg_enable_nonlocal_correction("C", "nonlocal-correction", "Enable the non-local-correction for the EPM model", false);
    TCLAP::ValueArg<int> arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    cmd.add(arg_mesh_file);
    cmd.add(arg_material);
    cmd.add(arg_nb_bands);
    cmd.add(arg_outfile);
    cmd.add(arg_nearest_neighbors);
    cmd.add(arg_nb_threads);
    cmd.add(arg_enable_nonlocal_correction);

    cmd.parse(argc, argv);

    EmpiricalPseudopotential::Materials materials;
    const std::string                   file_material_parameters = std::string(CMAKE_SOURCE_DIR) + "/parameter_files/materials.yaml";
    materials.load_material_parameters(file_material_parameters);

    Options my_options;
    my_options.materialName     = arg_material.getValue();
    my_options.nrLevels         = arg_nb_bands.getValue();
    my_options.nearestNeighbors = arg_nearest_neighbors.getValue();
    my_options.nrThreads        = arg_nb_threads.getValue();
    // my_options.print_options();

    EmpiricalPseudopotential::Material mat = materials.materials.at(my_options.materialName);

    // Create a new MPI type for the struct k_vector.
    MPI_Datatype k_vector_type;
    const int    number_item_k_vector = 3;
    MPI_Datatype type[3]              = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    int          block_lengths[3]     = {1, 1, 1};

    MPI_Aint offsets[3];
    offsets[0] = offsetof(vector_k, m_kx);
    offsets[1] = offsetof(vector_k, m_ky);
    offsets[2] = offsetof(vector_k, m_kz);

    MPI_Type_create_struct(number_item_k_vector, block_lengths, offsets, type, &k_vector_type);
    MPI_Type_commit(&k_vector_type);
    // END Creating a new MPI type for the struct k_vector.

    std::vector<vector_k> all_k_vectors;
    long                  number_k_vectors = 0;
    const std::string mesh_filename = arg_mesh_file.getValue();
    bz_mesh_points    my_mesh(mesh_filename);
    if (process_rank == MASTER) {
        my_mesh.read_mesh();
        std::vector<Vector3D<double>>& mesh_k_points = my_mesh.get_kpoints();
        number_k_vectors                             = mesh_k_points.size();
        all_k_vectors.resize(mesh_k_points.size());
        for (size_t i = 0; i < mesh_k_points.size(); i++) {
            all_k_vectors[i].set_k(mesh_k_points[i].X, mesh_k_points[i].Y, mesh_k_points[i].Z);
        }
    }
    MPI_Bcast(&number_k_vectors, 1, MPI_LONG, MASTER, MPI_COMM_WORLD);
    int number_bands = my_options.nrLevels;

    // Scatter the k_vectors to all processes.
    double t_start = MPI_Wtime();

    // Define the number of elements each process will handle.
    int count     = (number_k_vectors / number_processes);
    int remainder = (number_k_vectors % number_processes);

    std::vector<int> counts_element_per_process(number_processes);
    std::vector<int> displacements_element_per_process(number_processes);
    for (int i = 0; i < number_processes - 1; i++) {
        counts_element_per_process[i]        = count;
        displacements_element_per_process[i] = i * count;
    }
    counts_element_per_process.back()        = (count + remainder);
    displacements_element_per_process.back() = ((number_processes - 1) * count);

    std::vector<vector_k> chunk_vector_of_k;
    chunk_vector_of_k.resize(counts_element_per_process[process_rank]);

    std::cout << "Process " << process_rank << " will handle " << counts_element_per_process[process_rank] << " k-points" << std::endl;

    // Scatter the vector of structs.
    MPI_Scatterv(all_k_vectors.data(),
                 counts_element_per_process.data(),
                 displacements_element_per_process.data(),
                 k_vector_type,
                 chunk_vector_of_k.data(),
                 counts_element_per_process[process_rank],
                 k_vector_type,
                 MASTER,
                 MPI_COMM_WORLD);

    std::vector<Vector3D<double>> Chunk_list_k_points;
    Chunk_list_k_points.resize(counts_element_per_process[process_rank]);
    for (int i = 0; i < counts_element_per_process[process_rank]; ++i) {
        Chunk_list_k_points[i] = chunk_vector_of_k[i].to_Vector3D();
    }

    bool                                    enable_nonlocal_correction = false;
    EmpiricalPseudopotential::BandStructure my_bandstructure;
    my_bandstructure.Initialize(mat, my_options.nrLevels, Chunk_list_k_points, my_options.nearestNeighbors, enable_nonlocal_correction);
    my_bandstructure.Compute();
    my_bandstructure.AdjustValues();

    std::vector<double> chunk_list_energies(counts_element_per_process[process_rank] * my_options.nrLevels);
    for (int i = 0; i < counts_element_per_process[process_rank]; ++i) {
        for (int j = 0; j < number_bands; ++j) {
            chunk_list_energies[i * number_bands + j] = my_bandstructure.get_energy_at_k_band(j, i);
        }
    }

    std::vector<int> gather_counts_element_per_process(number_processes);
    std::vector<int> gather_displacements_element_per_process(number_processes);
    for (int i = 0; i < number_processes; i++) {
        gather_counts_element_per_process[i]        = counts_element_per_process[i] * number_bands;
        gather_displacements_element_per_process[i] = displacements_element_per_process[i] * number_bands;
    }

    // Gather the results
    std::vector<double> all_energies_all_bands;
    if (process_rank == MASTER) {
        all_energies_all_bands.resize(number_k_vectors * number_bands);
    }

    MPI_Gatherv(chunk_list_energies.data(),
                counts_element_per_process[process_rank] * number_bands,
                MPI_DOUBLE,
                all_energies_all_bands.data(),
                gather_counts_element_per_process.data(),
                gather_displacements_element_per_process.data(),
                MPI_DOUBLE,
                MASTER,
                MPI_COMM_WORLD);

    double t_end = MPI_Wtime();

    if (process_rank == MASTER) {
        double total_time = t_end - t_start;
        // std::ofstream f_time("mpiBands_times.csv", std::ios::app);
        // f_time << number_processes << "," << total_time << std::endl;
        // f_time.close();
        std::cout << "------------------------------------------------------------\n";
        std::cout << "TOTAL COMPUTATION TIME : " << total_time << std::endl;
        std::cout << "------------------------------------------------------------\n";
    }

    if (process_rank == MASTER) {
        std::filesystem::path in_path(mesh_filename);
        std::string out_file_bands = in_path.stem().replace_extension("").string() + "_MPI_" + my_bandstructure.path_band_filename();
        my_mesh.add_all_bands_on_mesh(out_file_bands + "_all_bands.msh", all_energies_all_bands, number_bands);
    }

    MPI_Finalize();
}
