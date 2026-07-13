/**
 * @file mmmc_particle_pools.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-07-10
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "device_admc_simulation.hpp"
#include "mmmc_policy.hpp"
#include "pbmc_particle.hpp"
#include "pbmc_transport_kernel.hpp"

namespace uepm::MMMC {

struct transfer_counters {
    std::size_t m_pbmc_to_admc = 0;
    std::size_t m_admc_to_pbmc = 0;

    void reset() noexcept;
};

class particle_pools {
 public:
    std::vector<std::unique_ptr<PBMC::pbmc_particle>> m_pbmc_particles;
    std::vector<ADMC::device_admc_particle>           m_admc_particles;
    transfer_counters                                 m_last_transfer_counters{};
    transfer_counters                                 m_total_transfer_counters{};

    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] std::size_t pbmc_size() const noexcept;
    [[nodiscard]] std::size_t admc_size() const noexcept;

    void apply_policy(const bbox_transport_policy& policy,
                      PBMC::pbmc_transport_kernel& electron_transport,
                      PBMC::pbmc_transport_kernel& hole_transport);
};

}  // namespace uepm::MMMC
