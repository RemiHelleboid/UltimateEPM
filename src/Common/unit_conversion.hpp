/**
 * @file unit_conversion.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-27
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <cmath>

namespace uepm::units {

inline constexpr double electric_field_V_per_cm_to_V_per_m = 1.0e2;
inline constexpr double meter_to_micron                    = 1.0e6;
inline constexpr double micron3_to_cm3                     = 1.0e-12;
inline constexpr double cm_to_micron                       = 1.0e4;

}  // namespace uepm::units