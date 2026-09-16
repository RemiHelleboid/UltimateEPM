/**
 * @file mmmc_policy.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-07-10
 * 
 * @copyright Copyright (c) 2026
 * 
 */

 
#include "mmmc_policy.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace uepm::MMMC {
namespace {

std::uint64_t splitmix64(std::uint64_t value) noexcept {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

}  // namespace

std::string_view transport_method_name(transport_method method) noexcept {
    switch (method) {
        case transport_method::pbmc:
            return "PBMC";
        case transport_method::admc:
            return "ADMC";
    }
    return "unknown";
}

void bbox_transport_policy::validate() const {
    if (!std::isfinite(m_pbmc_region_um.get_x_min()) || !std::isfinite(m_pbmc_region_um.get_x_max()) ||
        !std::isfinite(m_pbmc_region_um.get_y_min()) || !std::isfinite(m_pbmc_region_um.get_y_max()) ||
        !std::isfinite(m_pbmc_region_um.get_z_min()) || !std::isfinite(m_pbmc_region_um.get_z_max())) {
        throw std::invalid_argument("MMMC PBMC bbox coordinates must be finite.");
    }
    if (!m_pbmc_region_um.check_order()) {
        throw std::invalid_argument("MMMC PBMC bbox min/max coordinates are not ordered.");
    }
    if (!std::isfinite(m_buffer_width_um) || m_buffer_width_um < 0.0) {
        throw std::invalid_argument("MMMC PBMC bbox buffer width must be finite and non-negative.");
    }
}

double bbox_transport_policy::admc_probability_at(const mesh::vector3& position_um) const noexcept {
    if (!std::isfinite(position_um.x()) || !std::isfinite(position_um.y()) || !std::isfinite(position_um.z())) {
        return 1.0;
    }
    if (m_pbmc_region_um.is_inside(position_um)) {
        return 0.0;
    }
    if (m_buffer_width_um <= 0.0) {
        return 1.0;
    }

    const std::array<double, 6> distances_outside{
        m_pbmc_region_um.get_x_min() - position_um.x(),
        position_um.x() - m_pbmc_region_um.get_x_max(),
        m_pbmc_region_um.get_y_min() - position_um.y(),
        position_um.y() - m_pbmc_region_um.get_y_max(),
        m_pbmc_region_um.get_z_min() - position_um.z(),
        position_um.z() - m_pbmc_region_um.get_z_max(),
    };
    const double distance_into_buffer_um = std::max(0.0, *std::max_element(distances_outside.begin(),
                                                                           distances_outside.end()));
    const double buffer_coordinate = std::clamp(distance_into_buffer_um / m_buffer_width_um, 0.0, 1.0);
    return buffer_coordinate * buffer_coordinate * (3.0 - 2.0 * buffer_coordinate);
}

double bbox_transport_policy::particle_threshold(std::size_t particle_index, std::uint64_t seed) const noexcept {
    const std::uint64_t mixed = splitmix64(static_cast<std::uint64_t>(particle_index) ^ splitmix64(seed));
    constexpr double    inverse_2_pow_53 = 1.0 / 9007199254740992.0;
    return (static_cast<double>(mixed >> 11U) + 0.5) * inverse_2_pow_53;
}

transport_method bbox_transport_policy::method_for_position(const mesh::vector3& position_um) const {
    return m_pbmc_region_um.is_inside(position_um) ? transport_method::pbmc : transport_method::admc;
}

transport_method bbox_transport_policy::method_for_particle(const mesh::vector3& position_um,
                                                             std::size_t          particle_index,
                                                             std::uint64_t        seed) const noexcept {
    const double admc_probability = admc_probability_at(position_um);
    if (admc_probability <= 0.0) {
        return transport_method::pbmc;
    }
    if (admc_probability >= 1.0) {
        return transport_method::admc;
    }
    return admc_probability >= particle_threshold(particle_index, seed) ? transport_method::admc
                                                                        : transport_method::pbmc;
}

bool bbox_transport_policy::uses_pbmc_at(const mesh::vector3& position_um) const {
    return method_for_position(position_um) == transport_method::pbmc;
}

}  // namespace uepm::MMMC
