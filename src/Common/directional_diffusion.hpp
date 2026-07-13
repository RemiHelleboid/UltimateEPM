#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace uepm::statistics {

struct directional_diffusion {
    std::array<double, 3> mean_displacement_m{};
    std::array<double, 3> displacement_variance_m2{};
    std::array<double, 3> coefficient_m2_per_s{};
    std::size_t           sample_count = 0;
    double                observation_time_s = 0.0;
};

template <typename PositionAccessor>
directional_diffusion estimate_directional_diffusion(std::size_t      sample_count,
                                                     double           observation_time_s,
                                                     PositionAccessor displacement_m) {
    directional_diffusion result;
    result.sample_count       = sample_count;
    result.observation_time_s = observation_time_s;

    const double undefined = std::numeric_limits<double>::quiet_NaN();
    if (sample_count < 2 || !std::isfinite(observation_time_s) || observation_time_s <= 0.0) {
        result.displacement_variance_m2.fill(undefined);
        result.coefficient_m2_per_s.fill(undefined);
        return result;
    }

    std::array<double, 3> displacement_sum_m{};
    std::array<double, 3> displacement_squared_sum_m2{};
    for (std::size_t index = 0; index < sample_count; ++index) {
        const auto displacement = displacement_m(index);
        for (std::size_t axis = 0; axis < 3; ++axis) {
            displacement_sum_m[axis] += displacement[axis];
            displacement_squared_sum_m2[axis] += displacement[axis] * displacement[axis];
        }
    }

    for (std::size_t axis = 0; axis < 3; ++axis) {
        const double mean = displacement_sum_m[axis] / static_cast<double>(sample_count);
        const double centered_sum =
            displacement_squared_sum_m2[axis] - static_cast<double>(sample_count) * mean * mean;
        const double variance = std::max(0.0, centered_sum / static_cast<double>(sample_count - 1));
        result.mean_displacement_m[axis]     = mean;
        result.displacement_variance_m2[axis] = variance;
        result.coefficient_m2_per_s[axis]    = variance / (2.0 * observation_time_s);
    }
    return result;
}

}  // namespace uepm::statistics
