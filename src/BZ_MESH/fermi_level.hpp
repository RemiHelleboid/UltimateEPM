/**
 * @file fermi_level.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
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
#pragma once
/**
 * @file bz_fermi_tool.hpp
 * @brief DOS per band (CSV) + Fermi level solver (charge neutrality).
 *
 * Usage:
 *   #include "bz_fermi_tool.hpp"
 *   using namespace uepm::mesh_bz::fermi;
 *   Options opt; opt.nE = 2000; opt.T_K = 300.0; opt.write_csv = true; opt.csv_path = "DOS_si_mesh.csv";
 *   Result r = solve_fermi_and_export_csv(my_bz_mesh, opt);
 */

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
    bool        write_csv  = true;
    std::string csv_path   = "";     // if empty, a default name is chosen by the caller
    double      T_K        = 300.0;  // temperature (K)
    Dopants     dop;                 // dopant model (leave zeros for intrinsic)
};

struct Result {
    // Fermi solution
    double EF_eV    = 0.0;
    double n_m3     = 0.0;
    double p_m3     = 0.0;
    double Nd_plus  = 0.0;
    double Na_minus = 0.0;
    int    iters    = 0;
    bool   success  = false;

    // DOS arrays used (for plotting or reuse)
    std::vector<std::vector<double>> energies_per_band;  // eV
    std::vector<std::vector<double>> dos_per_band;       // states/(m^3·eV)

    // CSV path if a file was written
    std::optional<std::string> csv_written;
};

/**
 * @brief Compute per-band DOS using MeshBZ, optionally export CSV, and solve EF by charge neutrality.
 *
 * Pre-conditions:
 *   - `mesh` already loaded with geometry + bands.
 *   - `MeshBZ::compute_dos_at_energy_and_band` returns DOS in states/(m^3·eV) with
 *     symmetry normalization and spin degeneracy handled consistently (as in your library).
 */
Result solve_fermi_and_export_csv(MeshBZ& mesh, const Options& opt, const std::string& csv_path_if_empty = "");

}  // namespace uepm::mesh_bz::fermi
