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

std::vector<Vector3D<double>> read_qpoint_dat_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    std::vector<Vector3D<double>> kpoints;
    std::string                   line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string       token;
        std::vector<double> kpoint;
        while (std::getline(ss, token, ' ')) {
            kpoint.push_back(std::stod(token));
        }
        kpoints.push_back(Vector3D<double>(kpoint[0], kpoint[1], kpoint[2]));
    }
    return kpoints;
}

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

    // Get processor name
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int  name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    std::cout << "Process " << process_rank << " of " << number_processes << " is on " << processor_name << std::endl;

    // std::cout << "EPSILON PROGRAM" << std::endl;
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
    std::string outdir = config["outdir"].as<std::string>();

    std::string material_name        = config["material"].as<std::string>();
    int         nb_nearest_neighbors = config["nearest-neigbors"].as<int>();
    int         nb_bands             = config["nb-bands"].as<int>();


    double min_energy   = config["min-energy"].as<double>();
    double max_energy   = config["max-energy"].as<double>();
    double energy_step  = config["step-energy"].as<double>();
    double eta_smearing = config["eta-smearing"].as<double>();

    int Nkx = config["Nkx"].as<int>();
    int Nky = config["Nky"].as<int>();
    int Nkz = config["Nkz"].as<int>();

    bool nonlocal_epm = false;
    if (config["nonlocal"]) {
        nonlocal_epm = config["nonlocal"].as<bool>();
    } 


    if (process_rank == 0) {
        std::cout << "Material: " << material_name << std::endl;
        std::cout << "Number of nearest neighbors: " << nb_nearest_neighbors << std::endl;
        std::cout << "Number of bands: " << nb_bands << std::endl;
        std::cout << "Nonlocal corrections: " << nonlocal_epm << std::endl;
        std::cout << "Min energy: " << min_energy << std::endl;
        std::cout << "Max energy: " << max_energy << std::endl;
        std::cout << "Energy step: " << energy_step << std::endl;
        std::cout << "Eta smearing: " << eta_smearing << std::endl;
        std::cout << "Nkx: " << Nkx << std::endl;
        std::cout << "Nky: " << Nky << std::endl;
        std::cout << "Nkz: " << Nkz << std::endl;
    }

    int bz_sampling = config["bz-sampling"].as<int>();

    bool use_irreducible_wedge = (bz_sampling == 48) ? true : false;

    EmpiricalPseudopotential::Materials materials;
    const std::string                   file_material_parameters = std::string(CMAKE_SOURCE_DIR) + "/parameter_files/materials-local.yaml";
    materials.load_material_parameters(file_material_parameters);
    EmpiricalPseudopotential::Material      current_material = materials.materials.at("Si");
    EmpiricalPseudopotential::BandStructure band_structure{};

    band_structure.Initialize(current_material, nb_bands, {}, nb_nearest_neighbors, nonlocal_epm);
    EmpiricalPseudopotential::DielectricFunction MyDielectricFunc(current_material, band_structure.get_basis_vectors(), nb_bands);

    double shift = 0.0;
    MyDielectricFunc.generate_k_points_grid(Nkx, Nky, Nkz, shift, use_irreducible_wedge);
    std::size_t nb_k_points = MyDielectricFunc.get_kpoints().size();
    if (process_rank == 0) {
        MyDielectricFunc.export_kpoints(outdir + "/kpoints.csv");
        std::cout << "Number of k-points: " << nb_k_points << std::endl;
    }

    bool                irreducible_wedge = (bz_sampling == 1) ? true : false;
    std::vector<double> list_energy;
    for (double energy = min_energy; energy <= max_energy + energy_step; energy += energy_step) {
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
    if (config["file-list-q"]) {
        std::string file_list_q = config["file-list-q"].as<std::string>();
        std::cout << "File list q: " << file_list_q << std::endl;
        list_q = read_qpoint_dat_file(file_list_q);
    } else {
        double                        min_q      = 1.0e-12;
        double                        max_q_norm = 4.0;
        double                        step_q     = 0.1e4;
        double                        qx         = min_q;
        Vector3D<double>              q          = get_q(qx, crystal_dir);
        while (q.Length() <= max_q_norm + step_q) {
            list_q.push_back(q);
            qx += step_q;
            q = get_q(qx, crystal_dir);
        }
    }
    nb_qpoints = list_q.size();
    if (process_rank == 0) {
        std::cout << "Number of energies: " << list_energy.size() << std::endl;
        std::cout << "Number of energies: " << list_energy.size() << std::endl;
        std::cout << "Crystalo dir: " << arg_crystal_dir.getValue() << std::endl;
        std::cout << "Number of q points: " << nb_qpoints << std::endl;
    }

    // Define the number of k-points each process will be responsible for.
    std::vector<int> counts_kpoints_per_process(number_processes);
    std::vector<int> displacements_kpoints_per_process(number_processes);
    int              nb_points = nb_k_points;
    while (nb_points > 0) {
        int displacement = 0;
        for (int i = 0; i < number_processes; i++) {
            counts_kpoints_per_process[i]++;
            displacements_kpoints_per_process[i] = displacement;
            displacement += counts_kpoints_per_process[i];
            nb_points--;
            if (nb_points <= 0) {
                break;
            }
        }
    }

    std::cout << "Process " << process_rank << " will handle " << counts_kpoints_per_process[process_rank] << " k-points" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    MyDielectricFunc.set_export_prefix(outdir + "/" + current_material.get_name() + "_");
    MyDielectricFunc.set_qpoints(list_q);
    MyDielectricFunc.set_energies(list_energy);
    MyDielectricFunc.set_offset_k_index(displacements_kpoints_per_process[process_rank]);
    MyDielectricFunc.set_nb_kpoints(counts_kpoints_per_process[process_rank]);
    MyDielectricFunc.set_non_local_epm(nonlocal_epm);

    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
    MyDielectricFunc.compute_dielectric_function(eta_smearing);
    MyDielectricFunc.clear_eigen_states();

    MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
    double end = MPI_Wtime();
    if (process_rank == 0) {
        std::cout << "Total computational time: " << end - start << std::endl;
    }

    // Gather the results from all the processes.
    std::vector<double> dielectric_function_per_process_flat;
    if (process_rank == 0) {
        dielectric_function_per_process_flat.resize(number_processes * nb_qpoints * list_energy.size());
    }
    std::vector<double> flattened_dielectric_function = MyDielectricFunc.get_flat_dielectric_function();
    MPI_Gather(flattened_dielectric_function.data(),
               flattened_dielectric_function.size(),
               MPI_DOUBLE,
               dielectric_function_per_process_flat.data(),
               flattened_dielectric_function.size(),
               MPI_DOUBLE,
               0,
               MPI_COMM_WORLD);
    // Reconstruct the dielectric function.
    if (process_rank == 0) {
        std::vector<std::vector<std::vector<double>>> dielectric_function_results(number_processes);
        for (int i = 0; i < number_processes; i++) {
            dielectric_function_results[i].resize(nb_qpoints);
            for (int j = 0; j < nb_qpoints; j++) {
                dielectric_function_results[i][j].resize(list_energy.size());
                for (int k = 0; k < list_energy.size(); k++) {
                    dielectric_function_results[i][j][k] =
                        dielectric_function_per_process_flat[i * nb_qpoints * list_energy.size() + j * list_energy.size() + k];
                }
            }
        }
        // Merge the results.
        EmpiricalPseudopotential::DielectricFunction dielectric_function =
            EmpiricalPseudopotential::DielectricFunction::merge_results(MyDielectricFunc,
                                                                        dielectric_function_results,
                                                                        counts_kpoints_per_process);
        dielectric_function.apply_kramers_kronig();
        std::cout << "END" << std::endl;

        std::filesystem::create_directories(outdir);
        // Write the results to a file.
        dielectric_function.export_dielectric_function("", true);
    }

    MPI_Finalize();
    return 0;
}