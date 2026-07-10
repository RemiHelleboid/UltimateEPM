/**
 * @file mmmc_particle_transfer.hpp
 * @brief Particle-state handoff helpers between PBMC and ADMC pools.
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
