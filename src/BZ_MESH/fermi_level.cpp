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
#include <ostream>
#include <stdexcept>

#include "integrals.hpp"  // uepm::integrate::trapz
#include "physical_constants.hpp"
#include "physical_functions.hpp"  // uepm::physics::fermi_dirac_distribution

namespace uepm::mesh_bz::fermi {

/**
 * @brief Compute the number of electrons in a given band.
 *
 * @param E_eV
 * @param G_m3eV
 * @param EF_eV
 * @param T_K
 * @return double
 */
static inline double electrons_in_band(const std::vector<double>& E_eV, const std::vector<double>& G_m3eV, double EF_eV, double T_K) {
    std::vector<double> w;
    w.reserve(E_eV.size());
    for (std::size_t i = 0; i < E_eV.size(); ++i) {
        w.push_back(G_m3eV[i] * uepm::physics::fermi_dirac_distribution(E_eV[i], EF_eV, T_K));
    }
    return uepm::integrate::trapz(E_eV, w);  // states / m^3
}

/**
 * @brief Compute the number of holes in a given band.
 *
 * @param E_eV
 * @param G_m3eV
 * @param EF_eV
 * @param T_K
 * @return double
 */
static inline double holes_in_band(const std::vector<double>& E_eV, const std::vector<double>& G_m3eV, double EF_eV, double T_K) {
    std::vector<double> w;
    w.reserve(E_eV.size());
    for (std::size_t i = 0; i < E_eV.size(); ++i) {
        w.push_back(G_m3eV[i] * (1.0 - uepm::physics::fermi_dirac_distribution(E_eV[i], EF_eV, T_K)));
    }
    return uepm::integrate::trapz(E_eV, w);  // states / m^3
}

/**
 * @brief Compute the number of ionized donors at given EF.
 *
 * @param EF
 * @param Ec
 * @param d
 * @param T
 * @return double
 */
static inline double donors_ionized(double EF, double Ec, const Dopants& d, double T) {
    if (d.Nd_cm3 <= 0) {
        return 0.0;
    }
    const double kT  = uepm::Constants::k_b_eV * T;
    const double ED  = Ec - d.Ed_eV;
    const double occ = 1.0 / (1.0 + d.gd * std::exp((EF - ED) / kT));
    return d.Nd_cm3 * 1e6 * occ;  // cm^-3 -> m^-3
}

/**
 * @brief Compute the number of ionized acceptors at given EF.
 *
 * @param EF
 * @param Ev
 * @param d
 * @param T
 * @return double
 */
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
Result solve_fermi(MeshBZ& mesh, const Options& opt, bool use_iw) {
    Result results;

    // 1) Per-band DOS (ordered, thread-safe inside MeshBZ)
    const int nb_bands = static_cast<int>(mesh.get_number_bands_total());
    results.energies_per_band.reserve(nb_bands);
    results.dos_per_band.reserve(nb_bands);

    const auto valence_count    = mesh.get_band_indices(MeshParticleType::valence).size();
    const auto conduction_count = mesh.get_band_indices(MeshParticleType::conduction).size();
    std::cout << "Compute DOS on " << valence_count << " valence bands and " << conduction_count << " conduction bands.\n";
    std::cout << "Using " << opt.threads << " threads for DOS computation.\n";

    for (int b = 0; b < nb_bands; ++b) {
        auto lists = mesh.compute_dos_band_at_band_auto(b, static_cast<std::size_t>(opt.nE), opt.threads, opt.use_interp, use_iw);

        results.energies_per_band.push_back(std::move(lists[0]));
        results.dos_per_band.push_back(std::move(lists[1]));
    }

    const auto list_idx_val  = mesh.get_band_indices(MeshParticleType::valence);
    const auto list_idx_cond = mesh.get_band_indices(MeshParticleType::conduction);

    // 4) Estimate edges for dopant references
    double Ev = std::numeric_limits<double>::lowest();
    double Ec = std::numeric_limits<double>::max();
    for (int idx_band : list_idx_val) {
        const auto& E = results.energies_per_band[idx_band];
        if (!E.empty()) {
            Ev = std::max(Ev, *std::max_element(E.begin(), E.end()));
        }
    }
    for (int idx_band : list_idx_cond) {
        const auto& E = results.energies_per_band[idx_band];
        if (!E.empty()) {
            Ec = std::min(Ec, *std::min_element(E.begin(), E.end()));
        }
    }

    // 5) Neutrality function F(EF) = n - p + Nd+ - Na- (monotone ↑ in EF)
    auto F = [&](double EF) {
        double n = 0.0, p = 0.0;
        for (int idx_band : list_idx_cond) {
            n += electrons_in_band(results.energies_per_band[idx_band], results.dos_per_band[idx_band], EF, opt.T_K);
        }
        for (int idx_band : list_idx_val) {
            p += holes_in_band(results.energies_per_band[idx_band], results.dos_per_band[idx_band], EF, opt.T_K);
        }
        const double Nd = donors_ionized(EF, Ec, opt.dop, opt.T_K);
        const double Na = acceptors_ionized(EF, Ev, opt.dop, opt.T_K);
        return (n - p + Nd - Na);
    };

    // 6) Bracket EF and bisection
    double                low        = std::min(Ev, Ec) - 2.0;
    double                high       = std::max(Ev, Ec) + 2.0;
    double                Flo        = F(low);
    double                Fhi        = F(high);
    int                   expand     = 0;
    constexpr std::size_t max_expand = 12;
    while (Flo * Fhi > 0.0 && expand < max_expand) {
        low -= 1.0;
        high += 1.0;
        Flo = F(low);
        Fhi = F(high);
        ++expand;
    }
    if (Flo * Fhi > 0.0) {
        throw std::runtime_error("Fermi solve: could not bracket neutrality (same sign at ends).");
    }

    std::size_t           iter  = 0;
    constexpr std::size_t itMax = 100;
    double                mid   = 0.0;
    double                Fm    = 0.0;
    constexpr double      tol   = 1e-12;
    while (iter < itMax && (high - low) > tol) {  // ~1e-12 eV tolerance
        mid = 0.5 * (low + high);
        Fm  = F(mid);
        if (Fm == 0.0) {
            low = high = mid;
            break;
        }
        if (Flo * Fm < 0.0) {
            high = mid;
            Fhi  = Fm;
        } else {
            low = mid;
            Flo = Fm;
        }
        ++iter;
    }
    results.iters = iter;
    results.EF_eV = 0.5 * (low + high);
    if (iter < itMax) {
        results.success = true;
        std::cout << "Fermi level found at EF = " << results.EF_eV << " eV after " << iter << " iterations.\n";
    } else {
        results.success = false;
        std::cout << "Fermi level not converged after " << iter << " iterations. Last EF = " << results.EF_eV << " eV.\n";
    }

    // 7) Final carriers
    results.n_m3 = results.p_m3 = 0.0;
    for (int idx_band : list_idx_cond) {
        results.n_m3 += electrons_in_band(results.energies_per_band[idx_band], results.dos_per_band[idx_band], results.EF_eV, opt.T_K);
    }
    for (int idx_band : list_idx_val) {
        results.p_m3 += holes_in_band(results.energies_per_band[idx_band], results.dos_per_band[idx_band], results.EF_eV, opt.T_K);
    }
    results.Nd_plus  = donors_ionized(results.EF_eV, Ec, opt.dop, opt.T_K);
    results.Na_minus = acceptors_ionized(results.EF_eV, Ev, opt.dop, opt.T_K);

    return results;
}

}  // namespace uepm::mesh_bz::fermi
