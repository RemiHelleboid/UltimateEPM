/**
 * @file mmmc_particle_transfer.hpp
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

#include "admc_carrier.hpp"
#include "device_admc_simulation.hpp"
#include "pbmc_particle.hpp"
#include "pbmc_transport_kernel.hpp"

namespace uepm::MMMC {

[[nodiscard]] ADMC::carrier_type  to_admc_carrier_type(PBMC::particle_type type);
[[nodiscard]] PBMC::particle_type to_pbmc_particle_type(ADMC::carrier_type type);

[[nodiscard]] ADMC::device_admc_particle convert_pbmc_to_admc(const PBMC::pbmc_particle& particle,
                                                              std::size_t                new_index);

[[nodiscard]] PBMC::pbmc_particle convert_admc_to_pbmc(const ADMC::device_admc_particle& particle,
                                                       std::size_t                       new_index,
                                                       PBMC::pbmc_transport_kernel&      transport);

}  // namespace uepm::MMMC
