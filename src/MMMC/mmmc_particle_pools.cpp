/**
 * @file mmmc_particle_pools.cpp
 * @brief Two-pool particle ownership for Mixed Method Monte Carlo.
 */

#include "mmmc_particle_pools.hpp"

#include "mmmc_particle_transfer.hpp"
#include "unit_conversion.hpp"

namespace uepm::MMMC {
namespace {

mesh::vector3 admc_position_um(const ADMC::device_admc_particle& particle) {
    const auto& position_m = particle.particle.state().position_m;
    return {position_m.x() * units::meter_to_micron,
            position_m.y() * units::meter_to_micron,
            position_m.z() * units::meter_to_micron};
}

PBMC::pbmc_transport_kernel& pbmc_transport_for(ADMC::carrier_type           type,
                                                PBMC::pbmc_transport_kernel& electron_transport,
                                                PBMC::pbmc_transport_kernel& hole_transport) {
    return type == ADMC::carrier_type::electron ? electron_transport : hole_transport;
}

}  // namespace

void transfer_counters::reset() noexcept {
    m_pbmc_to_admc = 0;
    m_admc_to_pbmc = 0;
}

std::size_t particle_pools::size() const noexcept { return pbmc_size() + admc_size(); }

std::size_t particle_pools::pbmc_size() const noexcept { return m_pbmc_particles.size(); }

std::size_t particle_pools::admc_size() const noexcept { return m_admc_particles.size(); }

void particle_pools::apply_policy(const bbox_transport_policy& policy,
                                  PBMC::pbmc_transport_kernel& electron_transport,
                                  PBMC::pbmc_transport_kernel& hole_transport) {
    m_last_transfer_counters.reset();

    for (std::size_t i = 0; i < m_pbmc_particles.size();) {
        const auto& particle = *m_pbmc_particles[i];
        if (policy.method_for_position(particle.state().position) == transport_method::admc) {
            m_admc_particles.push_back(convert_pbmc_to_admc(particle, particle.index()));
            m_pbmc_particles[i] = std::move(m_pbmc_particles.back());
            m_pbmc_particles.pop_back();
            ++m_last_transfer_counters.m_pbmc_to_admc;
            continue;
        }
        ++i;
    }

    for (std::size_t i = 0; i < m_admc_particles.size();) {
        const auto& particle = m_admc_particles[i];
        if (policy.method_for_position(admc_position_um(particle)) == transport_method::pbmc) {
            auto& transport = pbmc_transport_for(particle.particle.type(), electron_transport, hole_transport);
            auto  converted = convert_admc_to_pbmc(particle, particle.particle.index(), transport);
            m_pbmc_particles.push_back(std::make_unique<PBMC::pbmc_particle>(std::move(converted)));
            m_admc_particles[i] = std::move(m_admc_particles.back());
            m_admc_particles.pop_back();
            ++m_last_transfer_counters.m_admc_to_pbmc;
            continue;
        }
        ++i;
    }

    m_total_transfer_counters.m_pbmc_to_admc += m_last_transfer_counters.m_pbmc_to_admc;
    m_total_transfer_counters.m_admc_to_pbmc += m_last_transfer_counters.m_admc_to_pbmc;
}

}  // namespace uepm::MMMC
