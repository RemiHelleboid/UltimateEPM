/**
 * @file dielectric_mesh.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-05-17
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "dielectric_mesh.hpp"

#include "gmsh.h"
#include <regex>


namespace bz_mesh {

void DielectricMesh::read_dielectric_file(const std::string& filename) {
    std::cout << "Opening file " << filename << std::endl;
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 0);
    gmsh::open(filename);
    std::vector<int> viewTags;
    gmsh::view::getTags(viewTags);
    int nb_energies =  viewTags.size();
    if (nb_energies % 2 != 0) {
        // The number of viewTags should be even because we have the real and imaginary part of the dielectric function for each energy.
        std::cerr << "The number of energies is not even. The file is corrupted." << std::endl;
        return;
    }
    nb_energies /= 2;
    int count_energy = 0;
    std::vector<double> energies;
    std::vector<std::vector<double>> real_dielectric;
    std::vector<std::vector<double>> imag_dielectric;
    for (auto&& tag : viewTags) {
        const int   index_view  = gmsh::view::getIndex(tag);
        std::string name_object = "View[" + std::to_string(index_view) + "].Name";
        std::string name_view;
        try {
            gmsh::option::getString(name_object, name_view);
        } catch (const std::exception& e) {
            std::cerr << e.what() << '\n';
        }

        std::string energy_str = name_view.substr(name_view.find_last_of("_") + 1, name_view.size());
        double      energy     = std::stod(energy_str);
        if (count_energy % 2 == 0) {
            energies.push_back(energy);
        }

        // If the name of the view starts with "eps_r" we store the data in the real part of the dielectric function.
        // If the name of the view starts with "eps_i" we store the data in the imaginary part of the dielectric function.
        std::regex               re_real("eps_r");
        std::regex               re_imag("eps_i");

        std::smatch match_real;
        std::smatch match_imag;
        std::string              type;
        std::vector<std::size_t> tags;
        double                   time;
        int                      numComp;
        std::vector<double>      data_view;
        gmsh::view::getHomogeneousModelData(tag, 0, type, tags, data_view, time, numComp);
        // std::cout << "View " << name_view << " has " << data_view.size() << " values." << std::endl;
        if (std::regex_search(name_view, match_real, re_real)) {
            real_dielectric.push_back(data_view);
        } else if (std::regex_search(name_view, match_imag, re_imag)) {
            imag_dielectric.push_back(data_view);
        }
        count_energy++;
    }
    std::cout << "Number of energies: " << nb_energies << std::endl;
    std::cout << "Number of energies in the real part: " << energies.size() << std::endl;
    std::cout << "Number of dielectric functions in the real part: " << real_dielectric.size() << std::endl;
    std::cout << "Number of dielectric functions in the imaginary part: " << imag_dielectric.size() << std::endl;
    if (energies.size() != nb_energies) {
        std::cerr << "The number of energies is not the same in the real and imaginary part of the dielectric function." << std::endl;
        return;
    }
    gmsh::finalize();
    m_energies = energies;
    m_dielectric_function.resize(real_dielectric.size());
    for (std::size_t idx_node = 0; idx_node < real_dielectric.size(); ++idx_node) {
        m_dielectric_function[idx_node].resize(energies.size());
        for (std::size_t idx_energy = 0; idx_energy < energies.size(); ++idx_energy) {
            m_dielectric_function[idx_node][idx_energy] = std::complex<double>(real_dielectric[idx_node][idx_energy], imag_dielectric[idx_node][idx_energy]);
        }
    }
    
}

}  // namespace bz_mesh