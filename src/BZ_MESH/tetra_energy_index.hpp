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

#include <algorithm>
#include <cstddef>
#include <span>
#include <vector>

#include "mesh_tetra.hpp"

namespace uepm::mesh_bz {

struct BandTetraEnergyIndex {
    std::vector<std::size_t> ordered_tetra_indices;
    std::vector<double>      ordered_min_energies;
    std::vector<double>      maximum_energy_tree;
    double                   maximum_energy_spread = 0.0;

    std::span<const std::size_t> candidate_indices(double minimum_energy, double maximum_energy) const noexcept;

    template <typename Function>
    void for_each_candidate(double minimum_energy, double maximum_energy, Function&& function) const {
        const auto upper =
            std::upper_bound(ordered_min_energies.begin(), ordered_min_energies.end(), maximum_energy);
        const std::size_t end_index = static_cast<std::size_t>(upper - ordered_min_energies.begin());
        if (end_index == 0 || maximum_energy_tree.empty()) {
            return;
        }

        const auto visit = [&](auto&& self,
                               std::size_t node,
                               std::size_t begin,
                               std::size_t end) -> void {
            if (begin >= end_index || maximum_energy_tree[node] < minimum_energy) {
                return;
            }
            if (end - begin == 1) {
                function(ordered_tetra_indices[begin]);
                return;
            }
            const std::size_t middle = begin + (end - begin) / 2;
            self(self, 2 * node + 1, begin, middle);
            self(self, 2 * node + 2, middle, end);
        };
        visit(visit, 0, 0, ordered_tetra_indices.size());
    }
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
