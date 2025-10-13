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
#include "integrals.hpp"

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
    TCLAP::ValueArg<int>         arg_nb_energies("e", "nenergy", "Number of energies to compute", false, 1000, "int");
    TCLAP::ValueArg<int>         arg_nb_conduction_bands("c", "cbands", "Number of conduction bands to consider", false, -1, "int");
    TCLAP::ValueArg<int>         arg_nb_valence_bands("v", "vbands", "Number of valence bands to consider", false, -1, "int");
    TCLAP::ValueArg<int>         arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    TCLAP::SwitchArg plot_with_python("P", "plot", "Call a python script after the computation to plot the band structure.", false);
    TCLAP::SwitchArg arg_test_interp("t", "test-interp", "Test the interpolation DOS.", false);
    cmd.add(plot_with_python);
    cmd.add(arg_mesh_file);
    cmd.add(arg_material);
    cmd.add(arg_nb_conduction_bands);
    cmd.add(arg_nb_valence_bands);
    cmd.add(arg_nb_energies);
    cmd.add(arg_nb_threads);
    cmd.add(arg_test_interp);

    cmd.parse(argc, argv);

    uepm::pseudopotential::Materials materials;
    const std::string                file_material_parameters = std::string(PROJECT_SRC_DIR) + "/parameter_files/materials-local.yaml";
    materials.load_material_parameters(file_material_parameters);

    Options my_options;
    my_options.materialName = arg_material.getValue();
    my_options.nrThreads    = arg_nb_threads.getValue();
    my_options.print_options();
    int  nb_conduction_bands = arg_nb_conduction_bands.getValue();
    int  nb_valence_bands    = arg_nb_valence_bands.getValue();
    bool use_interp      = arg_test_interp.getValue();
    auto start           = std::chrono::high_resolution_clock::now();

    uepm::pseudopotential::Material current_material = materials.materials.at(arg_material.getValue());

    const std::string mesh_band_input_file = arg_mesh_file.getValue();

    uepm::mesh_bz::MeshBZ my_bz_mesh{current_material};
    my_bz_mesh.read_mesh_geometry_from_msh_file(mesh_band_input_file);

    // Load bands: (nb_conduction_bands, nb_valence_bands)
    my_bz_mesh.read_mesh_bands_from_msh_file(mesh_band_input_file, nb_conduction_bands, nb_valence_bands);

    const int nb_bands = static_cast<int>(my_bz_mesh.get_number_bands_total());
    my_bz_mesh.print_band_info();

    std::cout << std::scientific;
    const double a     = current_material.get_lattice_constant_meter();
    const double Vcell = std::pow(a, 3) / 4.0;  // FCC primitive (m^3)
    const int    g_s   = 2;                     // spin degeneracy (Si)

    // (Keep symmetry-only factor inside your DOS implementation if using IBZ; do NOT do volume norm here.)
    std::cout << "Lattice constant a = " << a << " m\n";
    std::cout << "Mesh volume (mesh units): " << my_bz_mesh.compute_mesh_volume() << "\n";
    std::cout << "Vcell: " << Vcell << " m^3\n";

    std::vector<std::vector<double>> list_list_dos;
    std::vector<std::string>         list_header;
    list_list_dos.reserve(2 * nb_bands);
    list_header.reserve(2 * nb_bands);

    std::cout << "Compute DOS on " << nb_valence_bands << " valence bands and " << nb_conduction_bands << " conduction bands.\n";

    // Same number of energy samples per band keeps CSV rectangular:
    const int nE = arg_nb_energies.getValue();

    for (int band_index = 0; band_index < nb_bands; ++band_index) {
        auto lists_energies_dos = my_bz_mesh.compute_dos_band_at_band_auto(band_index, nE, my_options.nrThreads, use_interp);
        const auto  &E                  = lists_energies_dos[0];         // eV
        const auto  &G                  = lists_energies_dos[1];         // states / (m^3·eV)
        const double I                  = uepm::integrate::trapz(E, G);  // states / m^3
        const double target             = 2.0 / Vcell;                   // g_s / Vcell  (Si: g_s=2)
        std::cout << std::scientific << "raw ∫g(E)dE = " << I << "  target = " << target << "  ratio = " << (I / target)
                  << " × (g_s/Vcell)\n";

        // Append columns for this band (CSV expects pairs E_i, G_i)
        list_list_dos.push_back(lists_energies_dos[0]);
        list_list_dos.push_back(lists_energies_dos[1]);
        list_header.push_back("energy_band_" + std::to_string(band_index));
        list_header.push_back("dos_band_" + std::to_string(band_index));        
    }
    
    auto end              = std::chrono::high_resolution_clock::now();
    auto total_time_count = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Total computation time: " << total_time_count / double(1000) << std::endl;

    std::filesystem::path in_path(mesh_band_input_file);
    std::string           out_file_bands = "DOS_" + in_path.stem().replace_extension("").string();

    export_multiple_vector_to_csv(out_file_bands + ".csv", list_header, list_list_dos);

    const std::string python_plot_dos  = std::string(PROJECT_SRC_DIR) + "/python/plots/plot_density_of_states.py";
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