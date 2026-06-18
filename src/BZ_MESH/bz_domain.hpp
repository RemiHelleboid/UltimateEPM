/**
 * @file bz_domain.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-06-18
 * 
 * @copyright Copyright (c) 2026
 * 
 */

 
#pragma once

#include <array>
#include <cstdint>
#include <string_view>

#include "vector_bz.hpp"

namespace uepm::mesh_bz {

enum class BZDomainMode : std::uint8_t { full, positive_octant };

std::string_view bz_domain_mode_name(BZDomainMode mode) noexcept;

struct CanonicalK {
    vector3            representative;
    std::array<int, 3> signs{1, 1, 1};
};

CanonicalK canonicalize_k(const vector3& physical_k, BZDomainMode mode) noexcept;
vector3    apply_sign_image(const vector3& representative, const std::array<int, 3>& signs) noexcept;

inline constexpr std::array<std::array<int, 3>, 8> positive_octant_images = {
    std::array<int, 3>{1, 1, 1},
    std::array<int, 3>{1, 1, -1},
    std::array<int, 3>{1, -1, 1},
    std::array<int, 3>{1, -1, -1},
    std::array<int, 3>{-1, 1, 1},
    std::array<int, 3>{-1, 1, -1},
    std::array<int, 3>{-1, -1, 1},
    std::array<int, 3>{-1, -1, -1},
};

}  // namespace uepm::mesh_bz
