/**
 * @file fermi_level.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2025-10-13
 *
 *
 */

/**
 * @file bz_fermi_tool.cpp
 */

#include "fermi_level.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <stdexcept>

#include "integrals.hpp"  // uepm::integrate::trapz
#include "physical_constants.hpp"
#include "physical_functions.hpp"  // uepm::physics::fermi_dirac_distribution

namespace uepm::mesh_bz::fermi {

// ---------------- CSV helper ----------------
static void write_csv(const std::string& out_csv, const std::vector<std::string>& headers, const std::vector<std::vector<double>>& cols) {
    if (headers.size() != cols.size()) {
        throw std::runtime_error("CSV: headers/cols size mismatch");
    }
    std::ofstream out(out_csv);
    if (!out) {
        throw std::runtime_error("CSV: cannot open '" + out_csv + "'");
    }
    // header
    for (std::size_t j = 0; j < headers.size(); ++j) {
        out << headers[j] << (j + 1 < headers.size() ? "," : "\n");
    }
    // rows
    std::size_t rows = 0;
    for (auto& c : cols) {
        rows = std::max(rows, c.size());
    }
    out << std::setprecision(17);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols.size(); ++j) {
            if (i < cols[j].size()) {
                out << cols[j][i];
            }
            if (j + 1 < cols.size()) {
                out << ",";
            }
        }
        out << "\n";
    }
}

// -------------- carriers in a band --------------
static inline double electrons_in_band(const std::vector<double>& E_eV, const std::vector<double>& G_m3eV, double EF_eV, double T_K) {
    std::vector<double> w;
    w.reserve(E_eV.size());
    for (std::size_t i = 0; i < E_eV.size(); ++i) {
        w.push_back(G_m3eV[i] * uepm::physics::fermi_dirac_distribution(E_eV[i], EF_eV, T_K));
    }
    return uepm::integrate::trapz(E_eV, w);  // states / m^3
}

static inline double holes_in_band(const std::vector<double>& E_eV, const std::vector<double>& G_m3eV, double EF_eV, double T_K) {
    std::vector<double> w;
    w.reserve(E_eV.size());
    for (std::size_t i = 0; i < E_eV.size(); ++i) {
        w.push_back(G_m3eV[i] * (1.0 - uepm::physics::fermi_dirac_distribution(E_eV[i], EF_eV, T_K)));
    }
    return uepm::integrate::trapz(E_eV, w);  // states / m^3
}

// -------------- dopant ionization --------------
static inline double donors_ionized(double EF, double Ec, const Dopants& d, double T) {
    if (d.Nd_cm3 <= 0) {
        return 0.0;
    }
    const double kT  = uepm::Constants::k_b_eV * T;
    const double ED  = Ec - d.Ed_eV;
    const double occ = 1.0 / (1.0 + d.gd * std::exp((EF - ED) / kT));
    return d.Nd_cm3 * 1e6 * occ;  // cm^-3 -> m^-3
}

static inline double acceptors_ionized(double EF, double Ev, const Dopants& d, double T) {
    if (d.Na_cm3 <= 0) {
        return 0.0;
    }
    const double kT  = uepm::Constants::k_b_eV * T;
    const double EA  = Ev + d.Ea_eV;
    const double occ = 1.0 / (1.0 + d.ga * std::exp((EA - EF) / kT));
    return d.Na_cm3 * 1e6 * occ;  // cm^-3 -> m^-3
}

/**
 * @brief Compute the Fermi level and export DOS data to CSV.
 *
 * @param mesh
 * @param opt
 * @param csv_path_if_empty
 * @return Result
 */
Result solve_fermi_and_export_csv(MeshBZ& mesh, const Options& opt, const std::string& csv_path_if_empty) {
    Result res;

    // 1) Per-band DOS (ordered, thread-safe inside MeshBZ)
    const int nb_bands = static_cast<int>(mesh.get_number_bands_total());
    res.energies_per_band.reserve(nb_bands);
    res.dos_per_band.reserve(nb_bands);

    std::vector<std::vector<double>> csv_cols;
    std::vector<std::string>         csv_hdr;
    if (opt.write_csv) {
        csv_cols.reserve(2 * nb_bands);
        csv_hdr.reserve(2 * nb_bands);
    }

    const auto valence_count    = mesh.get_band_indices(MeshParticleType::valence).size();
    const auto conduction_count = mesh.get_band_indices(MeshParticleType::conduction).size();
    std::cout << "Compute DOS on " << valence_count << " valence bands and " << conduction_count << " conduction bands.\n";

    for (int b = 0; b < nb_bands; ++b) {
        auto lists = mesh.compute_dos_band_at_band_auto(b, static_cast<std::size_t>(opt.nE), opt.threads, opt.use_interp);

        res.energies_per_band.push_back(std::move(lists[0]));
        res.dos_per_band.push_back(std::move(lists[1]));

        if (opt.write_csv) {
            csv_cols.push_back(res.energies_per_band.back());
            csv_cols.push_back(res.dos_per_band.back());
            csv_hdr.push_back("energy_band_" + std::to_string(b));
            csv_hdr.push_back("dos_band_" + std::to_string(b));
        }
    }

    // 2) CSV export (optional)
    if (opt.write_csv) {
        const std::string out_csv =
            !opt.csv_path.empty() ? opt.csv_path : (!csv_path_if_empty.empty() ? csv_path_if_empty : "DOS_bz_mesh.csv");
        write_csv(out_csv, csv_hdr, csv_cols);
        std::cout << "Wrote CSV: " << out_csv << "\n";
        res.csv_written = out_csv;
    }

    // 3) Split bands via MeshBZ catalog
    std::vector<int> idx_val, idx_cond;
    {
        auto vidx = mesh.get_band_indices(MeshParticleType::valence);
        auto cidx = mesh.get_band_indices(MeshParticleType::conduction);
        idx_val.assign(vidx.begin(), vidx.end());
        idx_cond.assign(cidx.begin(), cidx.end());
    }

    // 4) Estimate edges for dopant references
    double Ev = -1e9, Ec = +1e9;
    for (int b : idx_val) {
        const auto& E = res.energies_per_band[b];
        if (!E.empty()) {
            Ev = std::max(Ev, *std::max_element(E.begin(), E.end()));
        }
    }
    for (int b : idx_cond) {
        const auto& E = res.energies_per_band[b];
        if (!E.empty()) {
            Ec = std::min(Ec, *std::min_element(E.begin(), E.end()));
        }
    }

    // 5) Neutrality function F(EF) = n - p + Nd+ - Na- (monotone ↑ in EF)
    auto F = [&](double EF) {
        double n = 0.0, p = 0.0;
        for (int b : idx_cond) {
            n += electrons_in_band(res.energies_per_band[b], res.dos_per_band[b], EF, opt.T_K);
        }
        for (int b : idx_val) {
            p += holes_in_band(res.energies_per_band[b], res.dos_per_band[b], EF, opt.T_K);
        }
        const double Nd = donors_ionized(EF, Ec, opt.dop, opt.T_K);
        const double Na = acceptors_ionized(EF, Ev, opt.dop, opt.T_K);
        return (n - p + Nd - Na);
    };

    // 6) Bracket EF and bisection
    double lo  = std::min(Ev, Ec) - 2.0;
    double hi  = std::max(Ev, Ec) + 2.0;
    double Flo = F(lo), Fhi = F(hi);
    int    expand = 0;
    while (Flo * Fhi > 0.0 && expand < 12) {
        lo -= 1.0;
        hi += 1.0;
        Flo = F(lo);
        Fhi = F(hi);
        ++expand;
    }
    if (Flo * Fhi > 0.0) {
        throw std::runtime_error("Fermi solve: could not bracket neutrality (same sign at ends).");
    }

    int    it    = 0;
    int    itMax = 120;
    double mid = 0.0, Fm = 0.0;
    while (it < itMax && (hi - lo) > 1e-8) {  // ~1e-8 eV tolerance
        mid = 0.5 * (lo + hi);
        Fm  = F(mid);
        if (Fm == 0.0) {
            lo = hi = mid;
            break;
        }
        if (Flo * Fm < 0.0) {
            hi  = mid;
            Fhi = Fm;
        } else {
            lo  = mid;
            Flo = Fm;
        }
        ++it;
    }
    res.iters = it;
    res.EF_eV = 0.5 * (lo + hi);
    if (it < itMax) {
        res.success = true;
        std::cout << "Fermi level found at EF = " << res.EF_eV << " eV after " << it << " iterations.\n";
    } else {
        res.success = false;
        std::cout << "Fermi level not converged after " << it << " iterations. Last EF = " << res.EF_eV << " eV.\n";
    }

    // 7) Final carriers
    res.n_m3 = res.p_m3 = 0.0;
    for (int b : idx_cond) {
        res.n_m3 += electrons_in_band(res.energies_per_band[b], res.dos_per_band[b], res.EF_eV, opt.T_K);
    }
    for (int b : idx_val) {
        res.p_m3 += holes_in_band(res.energies_per_band[b], res.dos_per_band[b], res.EF_eV, opt.T_K);
    }
    res.Nd_plus  = donors_ionized(res.EF_eV, Ec, opt.dop, opt.T_K);
    res.Na_minus = acceptors_ionized(res.EF_eV, Ev, opt.dop, opt.T_K);

    return res;
}

}  // namespace uepm::mesh_bz::fermi
