/**
 * @file band_catalog.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-06-18
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#include "band_catalog.hpp"

#include <numeric>

namespace uepm::mesh_bz {

void BandCatalog::register_band(MeshParticleType type, double minimum_energy, double maximum_energy) {
    classify_existing_band(type);
    m_minimum_energies.push_back(minimum_energy);
    m_maximum_energies.push_back(maximum_energy);
}

void BandCatalog::classify_existing_band(MeshParticleType type) {
    if (type == MeshParticleType::valence && m_conduction.count != 0) {
        throw std::runtime_error("Valence and conduction bands must form contiguous ranges");
    }

    BandRange& selected_range = range(type);
    if (selected_range.count == 0) {
        selected_range.global_start_index = static_cast<int>(m_total);
    }
    m_info.push_back(BandInfo{type, selected_range.count});
    ++selected_range.count;
    ++m_total;
}

std::vector<std::size_t> BandCatalog::indices(MeshParticleType type) const {
    const BandRange&         selected_range = range(type);
    std::vector<std::size_t> result(selected_range.count);
    std::iota(result.begin(), result.end(), static_cast<std::size_t>(selected_range.global_start_index));
    return result;
}

std::size_t BandCatalog::local_index(int global_band_index) const {
    if (global_band_index < 0 || global_band_index >= static_cast<int>(m_total)) {
        throw std::out_of_range("Global band index out of range");
    }
    return m_info.at(static_cast<std::size_t>(global_band_index)).local_index;
}

std::size_t BandCatalog::global_index(std::size_t local_band_index, MeshParticleType type) const {
    const BandRange& selected_range = range(type);
    if (local_band_index >= selected_range.count) {
        throw std::out_of_range("Local band index out of range");
    }
    return static_cast<std::size_t>(selected_range.global_start_index) + local_band_index;
}

std::pair<double, double> BandCatalog::extrema(std::size_t band_index) const {
    return {m_minimum_energies.at(band_index), m_maximum_energies.at(band_index)};
}

}  // namespace uepm::mesh_bz
