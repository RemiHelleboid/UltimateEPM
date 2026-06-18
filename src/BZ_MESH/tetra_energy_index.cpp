/**
 * @file tetra_energy_index.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-06-18
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#include "tetra_energy_index.hpp"

#include <fmt/core.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "numerical_helper.hpp"

namespace uepm::mesh_bz {

std::span<const std::size_t> BandTetraEnergyIndex::candidate_indices(double minimum_energy,
                                                                     double maximum_energy) const noexcept {
    const auto first       = std::lower_bound(ordered_min_energies.begin(),
                                        ordered_min_energies.end(),
                                        minimum_energy - maximum_energy_spread);
    const auto last        = std::upper_bound(ordered_min_energies.begin(), ordered_min_energies.end(), maximum_energy);
    const auto begin_index = static_cast<std::size_t>(first - ordered_min_energies.begin());
    const auto end_index   = static_cast<std::size_t>(last - ordered_min_energies.begin());
    return std::span<const std::size_t>(ordered_tetra_indices).subspan(begin_index, end_index - begin_index);
}

void TetraEnergyIndex::rebuild(const std::vector<Tetra>& tetrahedra,
                               std::size_t               number_of_bands,
                               double                    maximum_energy,
                               int                       number_of_threads) {
    if (number_of_threads <= 0) {
        throw std::invalid_argument("TetraEnergyIndex thread count must be positive");
    }

    fmt::print("Recomputing tetra ordered energies ...\n");
    m_bands.clear();
    m_bands.resize(number_of_bands);

    for (std::size_t band_index = 0; band_index < number_of_bands; ++band_index) {
        std::vector<double> minimum_energies(tetrahedra.size());
        std::vector<double> maximum_energies(tetrahedra.size());
#pragma omp parallel for schedule(dynamic) num_threads(number_of_threads)
        for (std::size_t tetra_index = 0; tetra_index < tetrahedra.size(); ++tetra_index) {
            minimum_energies[tetra_index] = tetrahedra[tetra_index].get_min_energy_at_band(band_index);
            maximum_energies[tetra_index] = tetrahedra[tetra_index].get_max_energy_at_band(band_index);
        }

        const auto sorted_indices  = uepm::numerical::argsort(minimum_energies);
        auto&      band            = m_bands[band_index];
        band.ordered_tetra_indices = sorted_indices;
        band.ordered_min_energies.resize(sorted_indices.size());
        for (std::size_t index = 0; index < sorted_indices.size(); ++index) {
            band.ordered_min_energies[index] = minimum_energies[sorted_indices[index]];
        }

        const auto last =
            std::upper_bound(band.ordered_min_energies.begin(), band.ordered_min_energies.end(), maximum_energy);
        const std::size_t retained = static_cast<std::size_t>(last - band.ordered_min_energies.begin());
        fmt::print(
            "\nBand {}: {} tetras have min energy <= {:.3f} eV (out of {} = {:.3f} %)\n",
            band_index,
            retained,
            maximum_energy,
            tetrahedra.size(),
            tetrahedra.empty() ? 0.0 : 100.0 * static_cast<double>(retained) / static_cast<double>(tetrahedra.size()));

        band.ordered_tetra_indices.resize(retained);
        band.ordered_min_energies.resize(retained);
        if (retained == 0) {
            band.maximum_energy_spread = 0.0;
            continue;
        }

        std::vector<double> energy_spreads(retained);
        for (std::size_t index = 0; index < retained; ++index) {
            const std::size_t tetra_index = band.ordered_tetra_indices[index];
            energy_spreads[index]         = std::abs(maximum_energies[tetra_index] - minimum_energies[tetra_index]);
        }

        const auto [minimum_spread, maximum_spread] = std::minmax_element(energy_spreads.begin(), energy_spreads.end());
        band.maximum_energy_spread                  = *maximum_spread;
        fmt::print("Band {}: min deltaE among tetras = {:.6f} eV, max deltaE = {:.6f} eV, at tetra #{}\n",
                   band_index,
                   *minimum_spread,
                   *maximum_spread,
                   static_cast<std::size_t>(maximum_spread - energy_spreads.begin()));
    }
    fmt::print("Done recomputing tetra ordered energies.\n");
}

}  // namespace uepm::mesh_bz
