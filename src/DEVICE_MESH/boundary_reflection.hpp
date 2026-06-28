/**
 * @file boundary_reflection.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-27
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>

#include "element.hpp"

namespace uepm::mesh {

enum class boundary_reflection_model { reverse, specular, diffuse };

inline std::string_view boundary_reflection_model_name(boundary_reflection_model model) {
    switch (model) {
        case boundary_reflection_model::reverse:
            return "reverse";
        case boundary_reflection_model::specular:
            return "specular";
        case boundary_reflection_model::diffuse:
            return "diffuse";
    }
    return "unknown";
}

inline boundary_reflection_model parse_boundary_reflection_model(const std::string& value) {
    if (value == "reverse" || value == "legacy" || value == "old") {
        return boundary_reflection_model::reverse;
    }
    if (value == "specular") {
        return boundary_reflection_model::specular;
    }
    if (value == "diffuse" || value == "diffusive") {
        return boundary_reflection_model::diffuse;
    }
    throw std::invalid_argument("boundary reflection must be one of: reverse, specular, diffuse.");
}

struct boundary_exit_hit {
    vector3 position{};
    vector3 inward_normal{};
    double  segment_fraction = 0.0;
};

inline vector3 normalized_or_zero(vector3 value) {
    const double norm = value.norm();
    if (norm == 0.0 || !std::isfinite(norm)) {
        return {};
    }
    return value / norm;
}

inline vector3 reflect_vector_specular(const vector3& value, const vector3& unit_normal) {
    return value - 2.0 * value.dot(unit_normal) * unit_normal;
}

inline vector3 place_reflected_position_inside(const element&            old_element,
                                               const vector3&            previous_position,
                                               const vector3&            trial_position,
                                               const boundary_exit_hit&  hit,
                                               const vector3&            outgoing_remaining_displacement,
                                               int                       dimension) {
    const double inward_epsilon = std::max(1.0e-9, 1.0e-9 * (trial_position - previous_position).norm());

    vector3 candidate = hit.position + outgoing_remaining_displacement + inward_epsilon * hit.inward_normal;
    if (dimension == 2) {
        candidate.to_2d_inplace();
    }
    if (old_element.is_location_inside_element(candidate)) {
        return candidate;
    }

    vector3 fallback = hit.position + inward_epsilon * hit.inward_normal;
    if (dimension == 2) {
        fallback.to_2d_inplace();
    }
    return fallback;
}

inline vector3 align_displacement_with_direction(const vector3& displacement, const vector3& direction) {
    const double displacement_norm = displacement.norm();
    const double direction_norm    = direction.norm();
    if (displacement_norm == 0.0 || direction_norm == 0.0 || !std::isfinite(direction_norm)) {
        return {};
    }
    return displacement_norm * direction / direction_norm;
}


inline std::optional<boundary_exit_hit> find_boundary_exit_hit(const element& old_element,
                                                               const vector3& previous_position,
                                                               const vector3& trial_position,
                                                               int            dimension) {
    const vector3 segment              = trial_position - previous_position;
    const double  segment_norm_squared = segment.norm_squared();
    if (segment_norm_squared == 0.0 || !std::isfinite(segment_norm_squared)) {
        return std::nullopt;
    }

    const auto intersections = old_element.compute_element_line_intersection(previous_position, trial_position);
    if (intersections.empty()) {
        return std::nullopt;
    }

    constexpr double min_tolerance = -1.0e-12;
    constexpr double max_tolerance = 1.0 + 1.0e-12;

    std::optional<boundary_exit_hit> best_hit;
    for (const auto& [face, intersection] : intersections) {
        const double t = (intersection - previous_position).dot(segment) / segment_norm_squared;
        if (t < min_tolerance || t > max_tolerance) {
            continue;
        }

        vector3 normal = face->compute_surface_normal();
        if (dimension == 2) {
            normal.to_2d_inplace();
        }
        normal = normalized_or_zero(normal);
        if (normal.norm_squared() == 0.0) {
            continue;
        }

        const vector3 old_element_interior = old_element.get_barycenter() - intersection;
        if (normal.dot(old_element_interior) < 0.0) {
            normal *= -1.0;
        }

        const boundary_exit_hit candidate{
            .position         = intersection,
            .inward_normal    = normal,
            .segment_fraction = std::clamp(t, 0.0, 1.0),
        };
        if (!best_hit.has_value() || candidate.segment_fraction < best_hit->segment_fraction) {
            best_hit = candidate;
        }
    }

    return best_hit;
}

template <typename UniformRandomBitGenerator>
vector3 draw_diffuse_reflection_vector(const vector3&             incoming,
                                       const vector3&             inward_unit_normal,
                                       int                        dimension,
                                       UniformRandomBitGenerator& rng) {
    const double magnitude = incoming.norm();
    if (magnitude == 0.0 || !std::isfinite(magnitude)) {
        return {};
    }

    std::uniform_real_distribution<double> unit_distribution(0.0, 1.0);

    if (dimension == 2) {
        vector3 normal = inward_unit_normal;
        normal.to_2d_inplace();
        normal = normalized_or_zero(normal);
        if (normal.norm_squared() == 0.0) {
            return -1.0 * incoming;
        }
        const vector3 tangent{-normal.y(), normal.x(), 0.0};
        const double  sin_theta = 2.0 * unit_distribution(rng) - 1.0;
        const double  cos_theta = std::sqrt(std::max(0.0, 1.0 - sin_theta * sin_theta));
        return magnitude * (cos_theta * normal + sin_theta * tangent);
    }

    vector3 helper   = std::abs(inward_unit_normal.x()) < 0.9 ? vector3{1.0, 0.0, 0.0} : vector3{0.0, 1.0, 0.0};
    vector3 tangent1 = normalized_or_zero(cross_product(inward_unit_normal, helper));
    if (tangent1.norm_squared() == 0.0) {
        tangent1 = {0.0, 0.0, 1.0};
    }
    const vector3 tangent2 = normalized_or_zero(cross_product(inward_unit_normal, tangent1));

    constexpr double two_pi = 6.283185307179586476925286766559;
    const double     u1     = unit_distribution(rng);
    const double     u2     = unit_distribution(rng);
    const double     r      = std::sqrt(u1);
    const double     phi    = two_pi * u2;
    const double     z      = std::sqrt(std::max(0.0, 1.0 - u1));

    return magnitude * (z * inward_unit_normal + r * std::cos(phi) * tangent1 + r * std::sin(phi) * tangent2);
}

}  // namespace uepm::mesh
