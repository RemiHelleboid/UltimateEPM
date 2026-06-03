/**
 * @file scattering_events.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-04-16
 *
 *
 */

#pragma once

#include <cstddef>
#include <string>

namespace uepm::amc {

enum class scattering_event : std::size_t {
    acoustic = 0,
    intervalley_absorption,
    intervalley_emission,
    impurity,
    impact_ionization,
    self_scattering,
    count
};

// std::string to_string(scattering_event event) {
//     switch (event) {
//         case scattering_event::acoustic:
//             return "acoustic";
//         case scattering_event::intervalley_absorption:
//             return "intervalley_absorption";
//         case scattering_event::intervalley_emission:
//             return "intervalley_emission";
//         case scattering_event::impact_ionization:
//             return "impact_ionization";
//         case scattering_event::self_scattering:
//             return "self_scattering";
//         default:
//             return "unknown";
//     }
// }

}  // namespace uepm::amc