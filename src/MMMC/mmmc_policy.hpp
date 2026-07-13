/**
 * @file mmmc_policy.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-07-10
 * 
 * @copyright Copyright (c) 2026
 * 
 */
#pragma once

#include <string_view>

#include "bbox.hpp"
#include "vector_mesh.hpp"

namespace uepm::MMMC {

enum class transport_method { pbmc, admc };

std::string_view transport_method_name(transport_method method) noexcept;

struct bbox_transport_policy {
    mesh::bbox m_pbmc_region_um{};

    void                           validate() const;
    [[nodiscard]] transport_method method_for_position(const mesh::vector3& position_um) const;
    [[nodiscard]] bool             uses_pbmc_at(const mesh::vector3& position_um) const;
};

}  // namespace uepm::MMMC
