/**
 * @file band_catalog.hpp
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
#include <stdexcept>
#include <utility>
#include <vector>

namespace uepm::mesh_bz {

enum class MeshParticleType { valence, conduction };

struct BandInfo {
    MeshParticleType type{MeshParticleType::conduction};
    std::size_t      local_index{0};
};

struct BandRange {
    int         global_start_index{-1};
    std::size_t count{0};
};

class BandCatalog {
 private:
    std::size_t           m_total = 0;
    std::vector<BandInfo> m_info;
    BandRange             m_valence;
    BandRange             m_conduction;
    std::vector<double>   m_minimum_energies;
    std::vector<double>   m_maximum_energies;

 public:
    void clear() noexcept {
        m_total = 0;
        m_info.clear();
        m_valence    = {};
        m_conduction = {};
        m_minimum_energies.clear();
        m_maximum_energies.clear();
    }
    void clear_classification() noexcept {
        m_total = 0;
        m_info.clear();
        m_valence    = {};
        m_conduction = {};
    }

    void register_band(MeshParticleType type, double minimum_energy, double maximum_energy);
    void classify_existing_band(MeshParticleType type);

    std::size_t total() const noexcept { return m_total; }
    void        set_total(std::size_t total) noexcept { m_total = total; }

    const BandRange& range(MeshParticleType type) const noexcept {
        return type == MeshParticleType::valence ? m_valence : m_conduction;
    }
    BandRange& range(MeshParticleType type) noexcept {
        return type == MeshParticleType::valence ? m_valence : m_conduction;
    }

    const std::vector<BandInfo>& info() const noexcept { return m_info; }
    std::vector<BandInfo>&       info() noexcept { return m_info; }

    const std::vector<double>& minima() const noexcept { return m_minimum_energies; }
    std::vector<double>&       minima() noexcept { return m_minimum_energies; }
    const std::vector<double>& maxima() const noexcept { return m_maximum_energies; }
    std::vector<double>&       maxima() noexcept { return m_maximum_energies; }

    std::vector<std::size_t>  indices(MeshParticleType type) const;
    std::size_t               local_index(int global_band_index) const;
    std::size_t               global_index(std::size_t local_band_index, MeshParticleType type) const;
    std::pair<double, double> extrema(std::size_t band_index) const;
};

}  // namespace uepm::mesh_bz
