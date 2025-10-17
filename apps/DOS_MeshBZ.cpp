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
#include "fermi_level.hpp"
#include "integrals.hpp"
#include "csv_utils.hpp"


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
    TCLAP::SwitchArg arg_use_iw("", "iw", "Use the irreducible wedge only.", false);
    cmd.add(plot_with_python);
    cmd.add(arg_mesh_file);
    cmd.add(arg_material);
    cmd.add(arg_nb_conduction_bands);
    cmd.add(arg_nb_valence_bands);
    cmd.add(arg_nb_energies);
    cmd.add(arg_nb_threads);
    cmd.add(arg_test_interp);
    cmd.add(arg_use_iw);

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
    bool use_interp          = arg_test_interp.getValue();
    bool use_iw              = arg_use_iw.getValue();
    auto start               = std::chrono::high_resolution_clock::now();

    uepm::pseudopotential::Material current_material = materials.materials.at(arg_material.getValue());

    const std::string mesh_band_input_file = arg_mesh_file.getValue();

    uepm::mesh_bz::MeshBZ my_bz_mesh{current_material};
    my_bz_mesh.set_nb_threads(my_options.nrThreads);
    my_bz_mesh.read_mesh_geometry_from_msh_file(mesh_band_input_file);
    my_bz_mesh.load_kstar_ibz_to_bz();

    bool set_positive_valence_band = false;
    bool shift_conduction_band     = true;
    my_bz_mesh.read_mesh_bands_from_msh_file(mesh_band_input_file,
                                             nb_conduction_bands,
                                             nb_valence_bands,
                                             shift_conduction_band,
                                             set_positive_valence_band);

    my_bz_mesh.apply_scissor(1.12);  // eV

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
    std::cout << "Using " << my_options.nrThreads << " threads.\n";
    if (use_interp) {
        std::cout << "Using interpolation when computing DOS.\n";
    } else {
        std::cout << "NOT using interpolation when computing DOS.\n";
    }
    if (use_iw) {
        std::cout << "Using irreducible wedge only when computing DOS.\n";
    } else {
        std::cout << "NOT using irreducible wedge only when computing DOS.\n";
    }

    // Same number of energy samples per band keeps CSV rectangular:
    const int nE = arg_nb_energies.getValue();

    for (int band_index = 0; band_index < nb_bands; ++band_index) {
        auto lists_energies_dos = my_bz_mesh.compute_dos_band_at_band_auto(band_index, nE, my_options.nrThreads, use_interp, use_iw);
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

    uepm::utils::export_multiple_vector_to_csv(out_file_bands + ".csv", list_header, list_list_dos);

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

    // Solve for Fermi level and export CSV
    uepm::mesh_bz::fermi::Options fermi_options;
    fermi_options.nE         = 1000;  // number of energy points for DOS interpolation
    fermi_options.threads    = my_options.nrThreads;
    fermi_options.use_interp = use_interp;  // use interpolation when computing DOS at given energy
    fermi_options.T_K        = 300.0;       // temperature for Fermi-Dirac

    auto result = uepm::mesh_bz::fermi::solve_fermi(my_bz_mesh, fermi_options, use_iw);
    if (result.success) {
        std::cout << "Fermi level found: EF = " << result.EF_eV << " eV\n";
        std::cout << "  p = " << result.p_m3 * 1e-6 << " cm^-3\n";
        std::cout << "  n = " << result.n_m3 * 1e-6 << " cm^-3\n";
    } else {
        std::cout << "Fermi level not found.\n";
    }

    return 0;
}