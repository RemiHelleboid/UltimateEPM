/**
 * @file tetra_energy_index.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-06-18
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "mesh_tetra.hpp"

namespace uepm::mesh_bz {

struct BandTetraEnergyIndex {
    std::vector<std::size_t> ordered_tetra_indices;
    std::vector<double>      ordered_min_energies;
    double                   maximum_energy_spread = 0.0;

    std::span<const std::size_t> candidate_indices(double minimum_energy, double maximum_energy) const noexcept;
};

class TetraEnergyIndex {
 private:
    std::vector<BandTetraEnergyIndex> m_bands;

 public:
    void rebuild(const std::vector<Tetra>& tetrahedra,
                 std::size_t               number_of_bands,
                 double                    maximum_energy,
                 int                       number_of_threads);

    void        clear() noexcept { m_bands.clear(); }
    std::size_t size() const noexcept { return m_bands.size(); }

    const BandTetraEnergyIndex& at(std::size_t band_index) const { return m_bands.at(band_index); }
    BandTetraEnergyIndex&       at(std::size_t band_index) { return m_bands.at(band_index); }
};

}  // namespace uepm::mesh_bz
