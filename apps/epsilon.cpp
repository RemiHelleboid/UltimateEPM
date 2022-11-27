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

#include <tclap/CmdLine.h>
#include <yaml-cpp/yaml.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <thread>

#include "BandStructure.h"
#include "DielectricFunction.hpp"
#include "Options.h"

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
        throw std::runtime_error("Invalid crystal direction");
    }
}

int main(int argc, char** argv) {
    std::cout << "EPSILON PROGRAM" << std::endl;
    TCLAP::CmdLine               cmd("Epsilon", ' ', "0.1");
    TCLAP::ValueArg<std::string> arg_yaml_config("c", "config", "YAML config file", true, "", "string");
    TCLAP::ValueArg<int>         arg_nb_threads("j", "threads", "Number of threads", false, 1, "int");
    TCLAP::ValueArg<int>         arg_crystal_dir("d", "dir", "Crystalographic direction (100, 110, 111)", false, 100, "int");

    cmd.add(arg_yaml_config);
    cmd.add(arg_nb_threads);
    cmd.add(arg_crystal_dir);
    cmd.parse(argc, argv);

    const std::string file_yaml_config = arg_yaml_config.getValue();
    const int         nb_threads       = arg_nb_threads.getValue();
    const int         crystal_dir      = arg_crystal_dir.getValue();

    YAML::Node config = YAML::LoadFile(file_yaml_config);
    // YAML::Node config = YAML::LoadFile("config.yaml");

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
    MyDielectricFunc.export_kpoints("kpoints.dat");

    bool                irreducible_wedge = (bz_sampling == 1) ? true : false;
    std::vector<double> list_energy;
    for (double energy = min_energy; energy <= max_energy; energy += energy_step) {
        list_energy.push_back(energy);
    }
    std::cout << "Number of energies: " << list_energy.size() << std::endl;
    std::cout << "Crystalo dir: " << arg_crystal_dir.getValue() << std::endl;
    double dqx = 1.0e-1;
    for (double qxyz = 1.0e-12; qxyz <= 1.0 + dqx; qxyz += dqx) {
        std::cout <<"\rScanning qxyz: " << qxyz * 100.0 << "%" << std::flush; 
        const Vector3D<double> q_vect       = get_q(qxyz, arg_crystal_dir.getValue());
        std::vector<double>    list_epsilon = MyDielectricFunc.compute_dielectric_function(q_vect, list_energy, eta_smearing, nb_threads);
        std::string export_dir = "Q" + std::to_string(crystal_dir) + "/";
        std::string export_filename = export_dir + "epsilon_dir" + std::to_string(arg_crystal_dir.getValue()) + "_Qx" + std::to_string(qxyz) + "_Nxyz" +
                                      std::to_string(Nkx) + ".csv";
        export_eps_result(export_filename, list_energy, list_epsilon, true);
    }
}