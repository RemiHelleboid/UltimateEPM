/**
 * @file fermi_level.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief Fermi level calculation on a BZ mesh.
 * @version 0.1
 * @date 2025-10-13
 *
 *
 */

#pragma once

/**
 * @file bz_fermi_tool.cpp
 */

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <optional>
#include <string>
#include <vector>

#include "bz_mesh.hpp"

namespace uepm::mesh_bz::fermi {

struct Dopants {
    double Nd_cm3 = 0.0;    // donors [cm^-3]
    double Na_cm3 = 0.0;    // acceptors [cm^-3]
    double Ed_eV  = 0.045;  // donor depth below CBM (E_C - E_D)
    double Ea_eV  = 0.045;  // acceptor height above VBM (E_A - E_V)
    double gd     = 0.5;    // donor degeneracy
    double ga     = 4.0;    // acceptor degeneracy
};

struct Options {
    std::size_t nE         = 2000;  // energy samples per band
    int         threads    = 8;     // passed to MeshBZ DOS calls
    bool        use_interp = true;  // MeshBZ: use interpolated tetra DOS
    double      T_K        = 300.0;  // temperature (K)
    Dopants     dop;                 // dopant model (leave zeros for intrinsic)
};

struct Result {
    double EF_eV    = 0.0;
    double n_m3     = 0.0;
    double p_m3     = 0.0;
    double Nd_plus  = 0.0;
    double Na_minus = 0.0;
    int    iters    = 0;
    bool   success  = false;

    /**
     * @brief Energy levels for each band (in eV).
     * 
     */
    std::vector<std::vector<double>> energies_per_band;

    /**
     * @brief DOS for each band (in states / (m^3·eV)).
     * 
     */
    std::vector<std::vector<double>> dos_per_band;

};

/**
 * @brief Solve the Fermi level using charge neutrality.
 * 
 * @param mesh 
 * @param opt 
 * @return Result 
 */
Result solve_fermi(MeshBZ& mesh, const Options& opt);

}  // namespace uepm::mesh_bz::fermi
