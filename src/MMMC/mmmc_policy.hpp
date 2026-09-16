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

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <stdexcept>
#include <string_view>

#include "bbox.hpp"
#include "vector_mesh.hpp"

namespace uepm::MMMC {

enum class transport_method { pbmc, admc };

enum class transport_policy { bbox, electric_field };

std::string_view transport_method_name(transport_method method) noexcept;

struct bbox_transport_policy {
    mesh::bbox m_pbmc_region_um{};
    double     m_buffer_width_um{0.0};

    void validate() const;
    [[nodiscard]] double admc_probability_at(const mesh::vector3& position_um) const noexcept;
    [[nodiscard]] double particle_threshold(std::size_t particle_index, std::uint64_t seed) const noexcept;
    [[nodiscard]] transport_method method_for_position(const mesh::vector3& position_um) const;
    [[nodiscard]] transport_method method_for_particle(const mesh::vector3& position_um,
                                                       std::size_t          particle_index,
                                                       std::uint64_t        seed) const noexcept;
    [[nodiscard]] bool uses_pbmc_at(const mesh::vector3& position_um) const;
};

struct electric_field_transport_policy {
    double m_pbmc_electric_field_threshold_V_per_cm{0.0};

    void validate() const {
        if (!std::isfinite(m_pbmc_electric_field_threshold_V_per_cm)) {
            throw std::invalid_argument("MMMC PBMC electric field threshold must be finite.");
        }
        if (m_pbmc_electric_field_threshold_V_per_cm < 0.0) {
            throw std::invalid_argument("MMMC PBMC electric field threshold must be non-negative.");
        }
    }

    [[nodiscard]] transport_method method_for_electric_field(double electric_field_V_per_cm) const {
        return (electric_field_V_per_cm <= m_pbmc_electric_field_threshold_V_per_cm) ? transport_method::pbmc : transport_method::admc;
    }

    [[nodiscard]] bool uses_pbmc_at(double electric_field_V_per_cm) const { return method_for_electric_field(electric_field_V_per_cm) == transport_method::pbmc; }
};

struct mmmc_policy {
    transport_policy m_transport_policy{transport_policy::bbox};
    bbox_transport_policy m_bbox_policy{};
    electric_field_transport_policy m_electric_field_policy{};

    void validate() const {
        switch (m_transport_policy) {
            case transport_policy::bbox:
                m_bbox_policy.validate();
                break;
            case transport_policy::electric_field:
                m_electric_field_policy.validate();
                break;
        }
    }
};

}  // namespace uepm::MMMC
