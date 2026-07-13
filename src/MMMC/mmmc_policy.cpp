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

#include <cmath>
#include <stdexcept>

namespace uepm::MMMC {

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
}

transport_method bbox_transport_policy::method_for_position(const mesh::vector3& position_um) const {
    return m_pbmc_region_um.is_inside(position_um) ? transport_method::pbmc : transport_method::admc;
}

bool bbox_transport_policy::uses_pbmc_at(const mesh::vector3& position_um) const {
    return method_for_position(position_um) == transport_method::pbmc;
}

}  // namespace uepm::MMMC
