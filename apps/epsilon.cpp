/**
 * @file epsilon.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-11-26
 *
 * @copyright Copyright (c) 2022
 *
 */

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

#include <mpi.h>
#include <tclap/CmdLine.h>
#include <yaml-cpp/yaml.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <thread>

#include "BandStructure.h"
#include "DielectricFunction.hpp"
#include "Options.h"

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

void export_eps_result(const std::string& filename, const std::vector<double> energies, const std::vector<double>& eps, bool python_plot) {
    std::ofstream file(filename);
    file << "Energy,Epsilon" << std::endl;
    for (std::size_t i = 0; i < energies.size(); ++i) {
        file << energies[i] << "," << eps[i] << std::endl;
    }
    file.close();
    const std::string python_plot_band_structure_script = std::string(CMAKE_SOURCE_DIR) + "/python/plots/plot_eps_vs_energy.py";
    std::string       python_call                       = "python3 " + python_plot_band_structure_script + " --filename " + filename;
    // bool              call_python_plot                  = false;
    // bool call_python_plot = true;
    if (python_plot) {
        std::cout << "Executing: " << python_call << std::endl;
        int succes_plot = system(python_call.c_str());
        std::cout << "Succes plot: " << succes_plot << std::endl;
    }
}

Vector3D<double> get_q(double qxyz, int crystalo_dir) {
    if (crystalo_dir == 100) {
        return Vector3D<double>(qxyz, 0.0, 0.0);
    } else if (crystalo_dir == 110) {
        return Vector3D<double>(qxyz, qxyz, 0.0);
    } else if (crystalo_dir == 111) {
        return Vector3D<double>(qxyz, qxyz, qxyz);
    } else {
        throw std::runtime_error("Invalid crystal direction: " + std::to_string(crystalo_dir));
    }
}

int main(int argc, char** argv) {
    int number_processes;
    int process_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &number_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

    std::cout << "Process rank: " << process_rank << std::endl;
    std::cout << "Number of processes: " << number_processes << std::endl;

    std::cout << "EPSILON PROGRAM" << std::endl;
    TCLAP::CmdLine               cmd("Epsilon", ' ', "0.1");
    TCLAP::ValueArg<std::string> arg_yaml_config("c", "config", "YAML config file", true, "", "string");
    TCLAP::ValueArg<int>         arg_crystal_dir("d", "dir", "Crystalographic direction (100, 110, 111)", false, 100, "int");

    cmd.add(arg_yaml_config);
    cmd.add(arg_crystal_dir);
    cmd.parse(argc, argv);

    const std::string file_yaml_config = arg_yaml_config.getValue();
    const int         crystal_dir      = arg_crystal_dir.getValue();

    YAML::Node config = YAML::LoadFile(file_yaml_config);

    if (!config["material"]) {
        std::cout << "No material section in the config file" << std::endl;
        exit(0);
    }
    std::string material_name = config["material"].as<std::string>();
    std::cout << "Material: " << material_name << std::endl;
    int nb_nearest_neighbors = config["nearest-neigbors"].as<int>();
    int nb_bands             = config["nb-bands"].as<int>();

    std::cout << "Number of nearest neighbors: " << nb_nearest_neighbors << std::endl;
    std::cout << "Number of bands: " << nb_bands << std::endl;
    bool nonlocal_corrections = config["nonlocal"].as<bool>();
    std::cout << "Nonlocal corrections: " << nonlocal_corrections << std::endl;

    double min_energy   = config["min-energy"].as<double>();
    double max_energy   = config["max-energy"].as<double>();
    double energy_step  = config["step-energy"].as<double>();
    double eta_smearing = config["eta-smearing"].as<double>();

    std::cout << "Min energy: " << min_energy << std::endl;
    std::cout << "Max energy: " << max_energy << std::endl;
    std::cout << "Energy step: " << energy_step << std::endl;
    std::cout << "Eta smearing: " << eta_smearing << std::endl;

    int Nkx = config["Nkx"].as<int>();
    int Nky = config["Nky"].as<int>();
    int Nkz = config["Nkz"].as<int>();
    std::cout << "Nkx: " << Nkx << std::endl;
    std::cout << "Nky: " << Nky << std::endl;
    std::cout << "Nkz: " << Nkz << std::endl;

    int bz_sampling = config["bz-sampling"].as<int>();

    bool use_irreducible_wedge = (bz_sampling == 48) ? true : false;

    EmpiricalPseudopotential::Materials materials;
    const std::string                   file_material_parameters = std::string(CMAKE_SOURCE_DIR) + "/parameter_files/materials-local.yaml";
    materials.load_material_parameters(file_material_parameters);
    EmpiricalPseudopotential::Material      current_material = materials.materials.at("Si");
    EmpiricalPseudopotential::BandStructure band_structure{};

    band_structure.Initialize(current_material, nb_bands, {}, nb_nearest_neighbors, nonlocal_corrections);
    EmpiricalPseudopotential::DielectricFunction MyDielectricFunc(current_material, band_structure.get_basis_vectors(), nb_bands);

    double shift = 1.0e-2;
    MyDielectricFunc.generate_k_points_grid(Nkx, Nky, Nkz, shift, use_irreducible_wedge);
    // MyDielectricFunc.generate_k_points_random(5000);
    // MyDielectricFunc.export_kpoints("kpoints.dat");

    bool                irreducible_wedge = (bz_sampling == 1) ? true : false;
    std::vector<double> list_energy;
    for (double energy = min_energy; energy <= max_energy; energy += energy_step) {
        list_energy.push_back(energy);
    }

    // Create a new MPI type for the struct k_vector.
    MPI_Datatype k_vector_type;
    const int    number_item_k_vector = 3;
    MPI_Datatype type[3]              = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    int          block_lengths[3]     = {1, 1, 1};
    MPI_Aint     offsets[3];
    offsets[0] = offsetof(vector_k, m_kx);
    offsets[1] = offsetof(vector_k, m_ky);
    offsets[2] = offsetof(vector_k, m_kz);
    MPI_Type_create_struct(number_item_k_vector, block_lengths, offsets, type, &k_vector_type);
    MPI_Type_commit(&k_vector_type);
    // END Creating a new MPI type for the struct k_vector.

    std::size_t                   nb_qpoints;
    std::vector<Vector3D<double>> list_q;
    if (process_rank == 0) {
        std::cout << "Number of energies: " << list_energy.size() << std::endl;
        double min_q  = 1.0e-12;
        double max_q  = 1.0;
        double step_q = 0.01;
        for (double qx = min_q; qx <= max_q; qx += step_q) {
            Vector3D<double> q = get_q(qx, crystal_dir);
            list_q.push_back(q);
        }

        nb_qpoints = list_q.size();
    }
    std::cout << "Number of energies: " << list_energy.size() << std::endl;
    std::cout << "Crystalo dir: " << arg_crystal_dir.getValue() << std::endl;
    std::cout << "Number of q points: " << nb_qpoints << std::endl;

    MPI_Bcast(&nb_qpoints, 1, MPI_LONG, MASTER, MPI_COMM_WORLD);

    // Define the number of elements each process will handle.
    int count     = (nb_qpoints / number_processes);
    int remainder = (nb_qpoints % number_processes);
    std::cout << "Count: " << count << std::endl;

    std::vector<int> counts_element_per_process(number_processes);
    std::vector<int> displacements_element_per_process(number_processes);
    for (int i = 0; i < number_processes - 1; i++) {
        counts_element_per_process[i]        = count;
        displacements_element_per_process[i] = i * count;
    }
    counts_element_per_process.back()        = (count + remainder);
    displacements_element_per_process.back() = ((number_processes - 1) * count);

    std::vector<vector_k> chunk_vector_of_q;
    chunk_vector_of_q.resize(counts_element_per_process[process_rank]);
    std::cout << "Process " << process_rank << " will handle " << counts_element_per_process[process_rank] << " q-points" << std::endl;

    // Scatter the q-points to each process.
    MPI_Scatterv(list_q.data(),
                 counts_element_per_process.data(),
                 displacements_element_per_process.data(),
                 k_vector_type,
                 chunk_vector_of_q.data(),
                 counts_element_per_process[process_rank],
                 k_vector_type,
                 MASTER,
                 MPI_COMM_WORLD);

    std::vector<Vector3D<double>> chunk_list_q;
    for (auto& q : chunk_vector_of_q) {
        chunk_list_q.push_back(Vector3D<double>{q.m_kx, q.m_ky, q.m_kz});
    }
    for (auto& q_vect : chunk_list_q) {
        std::cout << "Process " << process_rank << " q: " << q_vect << std::endl;
        std::vector<double> list_epsilon    = MyDielectricFunc.compute_dielectric_function(q_vect, list_energy, eta_smearing);
        std::string         export_dir      = "Q" + std::to_string(crystal_dir) + "/";
        std::string         export_filename = export_dir + "epsilon_dir" + std::to_string(arg_crystal_dir.getValue()) + "_Qx" +
                                      std::to_string(q_vect.Length()) + "_Nxyz" + std::to_string(Nkx) + ".csv";
        bool run_python_script = false;
        export_eps_result(export_filename, list_energy, list_epsilon, run_python_script);
    }

    MPI_Finalize();
    return 0;
}