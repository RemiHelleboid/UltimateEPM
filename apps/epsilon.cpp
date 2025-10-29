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
        if (line.empty()) {
            continue;
        }
        std::stringstream   ss(line);
        std::string         token;
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
    const std::string python_plot_band_structure_script = std::string(PROJECT_SRC_DIR) + "/python/plots/plot_eps_vs_energy.py";
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

    if (process_rank == 0) {
        std::cout << "EPSILON PROGRAM" << std::endl;
        std::cout << "Number of processes: " << number_processes << std::endl;
    }
    std::cout << "Process " << process_rank << " of " << number_processes << " is on " << processor_name << std::endl;

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
    bool enable_soc   = false;
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

    uepm::pseudopotential::Materials materials;
    std::string                      file_material_parameters = std::string(PROJECT_SRC_DIR) + "/parameter_files/materials-chel.yaml";
    if (nonlocal_epm) {
        file_material_parameters = std::string(PROJECT_SRC_DIR) + "/parameter_files/materials.yaml";
    }
    std::cout << "Loading material parameters from " << file_material_parameters << std::endl;

    materials.load_material_parameters(file_material_parameters);
    uepm::pseudopotential::Material      current_material = materials.materials.at("Si");
    uepm::pseudopotential::BandStructure band_structure{};

    band_structure.Initialize(current_material, nb_bands, {}, nb_nearest_neighbors, nonlocal_epm, enable_soc);
    uepm::pseudopotential::DielectricFunction MyDielectricFunc(current_material, band_structure.get_basis_vectors(), nb_bands);

    double shift = 0.0;
    MyDielectricFunc.generate_k_points_grid(Nkx, Nky, Nkz, shift, use_irreducible_wedge);
    std::size_t nb_k_points = MyDielectricFunc.get_kpoints().size();
    if (process_rank == 0) {
        MyDielectricFunc.export_kpoints(outdir + "/kpoints.csv");
        std::cout << "Number of k-points: " << nb_k_points << std::endl;
    }

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
        double           min_q      = 5.0e-13;
        double           max_q_norm = 4.0;
        double           step_q     = 0.1e4;
        double           qx         = min_q;
        Vector3D<double> q          = get_q(qx, crystal_dir);
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
    std::vector<int> counts_kpoints_per_process(number_processes, 0);
    std::vector<int> displacements_kpoints_per_process(number_processes, 0);
    int              nb_points = nb_k_points;
    const int        Ntot      = static_cast<int>(nb_k_points);
    const int        base      = Ntot / number_processes;
    const int        rem       = Ntot % number_processes;

    for (int p = 0; p < number_processes; ++p) {
        counts_kpoints_per_process[p]        = base + (p < rem ? 1 : 0);
        displacements_kpoints_per_process[p] = (p == 0) ? 0 : displacements_kpoints_per_process[p - 1] + counts_kpoints_per_process[p - 1];
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
    MyDielectricFunc.compute_dielectric_function(eta_smearing, process_rank);
    MyDielectricFunc.clear_eigen_states();

    MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
    double end = MPI_Wtime();
    if (process_rank == 0) {
        std::cout << "Total computational time: " << end - start << std::endl;
    }

    // Gather the results from all the processes.
    // --- Gather results robustly (works for equal or variable payload sizes) ---
    const std::size_t Q  = nb_qpoints;
    const std::size_t E  = list_energy.size();
    const std::size_t QE = Q * E;

    std::vector<double> flat_local = MyDielectricFunc.get_flat_dielectric_function();
    const int           local_n    = static_cast<int>(flat_local.size());

    // 1) Gather per-rank sizes to rank 0
    std::vector<int> recvcounts, displs;
    if (process_rank == 0) {
        recvcounts.resize(number_processes, 0);
    }

    MPI_Gather(&local_n, 1, MPI_INT, process_rank == 0 ? recvcounts.data() : nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // 2) Allocate receive buffer on rank 0 and compute displacements
    std::vector<double> all_flat;
    if (process_rank == 0) {
        displs.resize(number_processes, 0);
        for (int p = 1; p < number_processes; ++p) {
            displs[p] = displs[p - 1] + recvcounts[p - 1];
        }
        const int total = (number_processes > 0) ? (displs.back() + recvcounts.back()) : 0;
        all_flat.resize(static_cast<std::size_t>(total));
    }

    // 3) Gatherv the payloads
    MPI_Gatherv(flat_local.data(),
                local_n,
                MPI_DOUBLE,
                process_rank == 0 ? all_flat.data() : nullptr,
                process_rank == 0 ? recvcounts.data() : nullptr,
                process_rank == 0 ? displs.data() : nullptr,
                MPI_DOUBLE,
                0,
                MPI_COMM_WORLD);

    // 4) Reconstruct a [process][q][e] cube on rank 0
    if (process_rank == 0) {
        std::vector<std::vector<std::vector<double>>> dielectric_function_results(
            number_processes,
            std::vector<std::vector<double>>(Q, std::vector<double>(E, 0.0)));

        for (int p = 0; p < number_processes; ++p) {
            const int         count = recvcounts[p];
            const std::size_t base  = static_cast<std::size_t>(displs[p]);

            if (count == static_cast<int>(QE)) {
                // Case A: rank p already reduced over its local k's (exactly Q*E values)
                for (std::size_t q = 0; q < Q; ++q) {
                    const double* src = &all_flat[base + q * E];
                    std::copy(src, src + E, dielectric_function_results[p][q].data());
                }
            } else if (count % static_cast<int>(QE) == 0) {
                // Case B: rank p sent local_k blocks of size Q*E -> sum over k
                const int k_loc = count / static_cast<int>(QE);
                for (int k = 0; k < k_loc; ++k) {
                    const std::size_t block = base + static_cast<std::size_t>(k) * QE;
                    for (std::size_t q = 0; q < Q; ++q) {
                        const double* src = &all_flat[block + q * E];
                        double*       dst = dielectric_function_results[p][q].data();
                        for (std::size_t e = 0; e < E; ++e) {
                            dst[e] += src[e];
                        }
                    }
                }
                // If your per-k values are averages instead of sums, divide here by k_loc.
                // for (std::size_t q = 0; q < Q; ++q)
                //     for (std::size_t e = 0; e < E; ++e) dielectric_function_results[p][q][e] /= k_loc;
            } else {
                throw std::runtime_error("Unexpected payload size from rank " + std::to_string(p));
            }
        }

        // 5) Merge and finish
        uepm::pseudopotential::DielectricFunction dielectric_function =
            uepm::pseudopotential::DielectricFunction::merge_results(MyDielectricFunc,
                                                                     dielectric_function_results,
                                                                     counts_kpoints_per_process);

        dielectric_function.apply_kramers_kronig();
        std::filesystem::create_directories(outdir);
        dielectric_function.export_dielectric_function("", true);
        std::cout << "END" << std::endl;
    }

    MPI_Finalize();
    return 0;
}