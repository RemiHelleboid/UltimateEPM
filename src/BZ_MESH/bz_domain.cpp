/**
 * @file bz_domain.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-06-18
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#include "bz_domain.hpp"


#include <cmath>
namespace uepm::mesh_bz {

std::string_view bz_domain_mode_name(BZDomainMode mode) noexcept {
    return mode == BZDomainMode::positive_octant ? "positive-octant" : "full";
}

CanonicalK canonicalize_k(const vector3& physical_k, BZDomainMode mode) noexcept {
    if (mode == BZDomainMode::full) {
        return {physical_k, {1, 1, 1}};
    }

    const std::array<int, 3> signs{
        physical_k.x() < 0.0 ? -1 : 1,
        physical_k.y() < 0.0 ? -1 : 1,
        physical_k.z() < 0.0 ? -1 : 1,
    };
    return {vector3{std::abs(physical_k.x()), std::abs(physical_k.y()), std::abs(physical_k.z())}, signs};
}

vector3 apply_sign_image(const vector3& representative, const std::array<int, 3>& signs) noexcept {
    return {static_cast<double>(signs[0]) * representative.x(),
            static_cast<double>(signs[1]) * representative.y(),
            static_cast<double>(signs[2]) * representative.z()};
}

}  // namespace uepm::mesh_bz
