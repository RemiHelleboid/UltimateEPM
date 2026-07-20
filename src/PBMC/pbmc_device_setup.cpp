/**
 * @file pbmc_device_setup.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-09
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "pbmc_device_setup.hpp"

#include <fmt/core.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <optional>
#include <stdexcept>

namespace uepm::PBMC {

impurity_scattering_model parse_impurity_model(const std::string& text) {
    if (text == "mobility") {
        return impurity_scattering_model::mobility_empirical;
    }
    if (text == "screened-coulomb") {
        return impurity_scattering_model::screened_coulomb;
    }
    throw std::invalid_argument("--impurity-model must be either mobility or screened-coulomb.");
}

impurity_screening_model parse_impurity_screening_model(const std::string& text) {
    if (text == "debye") {
        return impurity_screening_model::debye_analytic;
    }
    if (text == "full") {
        return impurity_screening_model::finite_temperature_full;
    }
    throw std::invalid_argument("--impurity-screening must be either debye or full.");
}

std::string_view impurity_screening_model_name(impurity_screening_model model) {
    switch (model) {
        case impurity_screening_model::debye_analytic:
            return "debye-analytic";
        case impurity_screening_model::finite_temperature_full:
            return "finite-temperature-full";
    }
    throw std::runtime_error("unknown impurity screening model");
}

contact_injection_distribution parse_contact_injection_distribution(const std::string& text) {
    if (text == "maxwellian") {
        return contact_injection_distribution::maxwellian;
    }
    if (text == "velocity_weighted_maxwellian") {
        return contact_injection_distribution::velocity_weighted_maxwellian;
    }
    throw std::invalid_argument(
        "particles.contact_injection_distribution must be either maxwellian or "
        "velocity_weighted_maxwellian.");
}

std::string_view contact_injection_distribution_name(contact_injection_distribution distribution) {
    switch (distribution) {
        case contact_injection_distribution::maxwellian:
            return "maxwellian";
        case contact_injection_distribution::velocity_weighted_maxwellian:
            return "velocity_weighted_maxwellian";
    }
    throw std::runtime_error("unknown contact injection distribution");
}

particle_type parse_particle_type(const std::string& text) {
    if (text == "electron" || text == "e") {
        return particle_type::electron;
    }
    if (text == "hole" || text == "h") {
        return particle_type::hole;
    }
    throw std::invalid_argument("Invalid --inject-type. Expected one of: electron, e, hole, h.");
}

void validate_material_symbol(const std::string& material_symbol) {
    if (material_symbol != "Si") {
        throw std::invalid_argument("Only Si is currently supported by the analytical PBMC transport model.");
    }
}

std::string make_default_output_directory(const std::string& mesh_file) {
    return fmt::format("self_consistent_pbmc_{}", std::filesystem::path(mesh_file).stem().string());
}

void add_collecting_contacts(uepm::device::device&           simulation_device,
                             uepm::mesh::mesh&               mesh,
                             const std::vector<std::string>& contact_names) {
    constexpr double contact_margin   = 10.0;  // µm
    constexpr double ohmic_resistance = 0.0;

    const mesh::bbox device_bbox             = mesh.get_bounding_box();
    const double     tolerance               = 1.0e-8 * std::max(device_bbox.get_diagonal_size(), 1.0);
    const double     contact_collection_depth = tolerance;
    const auto       bulk_elements           = mesh.get_list_bulk_element();

    for (const auto& contact_name : contact_names) {
        const auto* region = mesh.get_p_region(contact_name);
        if (region == nullptr) {
            throw std::runtime_error("Collecting contact region '" + contact_name + "' does not exist in the mesh.");
        }
        if (region->get_region_type() != mesh::RegionType::contact) {
            throw std::runtime_error("Collecting contact '" + contact_name + "' is not a mesh contact region.");
        }

        const mesh::bbox            region_bbox = region->compute_bounding_box();
        const std::array<double, 3> region_min{region_bbox.get_x_min(),
                                               region_bbox.get_y_min(),
                                               region_bbox.get_z_min()};
        const std::array<double, 3> region_max{region_bbox.get_x_max(),
                                               region_bbox.get_y_max(),
                                               region_bbox.get_z_max()};
        std::array<double, 3>       box_min{};
        std::array<double, 3>       box_max{};
        std::optional<std::size_t>  normal_axis;

        for (std::size_t axis = 0; axis < 3; ++axis) {
            const double region_size = region_max[axis] - region_min[axis];
            const double device_size = axis == 0   ? device_bbox.get_x_size()
                                       : axis == 1 ? device_bbox.get_y_size()
                                                   : device_bbox.get_z_size();
            if (device_size > tolerance && region_size <= tolerance) {
                if (normal_axis.has_value()) {
                    throw std::runtime_error("Cannot determine a unique normal for collecting contact '" +
                                             contact_name + "'.");
                }
                normal_axis = axis;
            } else if (device_size <= tolerance) {
                box_min[axis] = region_min[axis] - contact_margin;
                box_max[axis] = region_max[axis] + contact_margin;
            } else {
                box_min[axis] = region_min[axis] - tolerance;
                box_max[axis] = region_max[axis] + tolerance;
            }
        }

        if (!normal_axis.has_value()) {
            throw std::runtime_error("Cannot determine the normal of collecting contact '" + contact_name + "'.");
        }

        const auto adjacent_indices = mesh.get_idx_bulk_elements_adjacent_to_contact_region(contact_name);
        if (adjacent_indices.empty()) {
            throw std::runtime_error("Collecting contact '" + contact_name + "' has no adjacent bulk element.");
        }

        double mean_adjacent_coordinate = 0.0;
        for (const auto element_index : adjacent_indices) {
            if (element_index >= bulk_elements.size()) {
                throw std::runtime_error("Contact-adjacent bulk element index is out of range.");
            }
            const auto barycenter = bulk_elements[element_index]->get_barycenter();
            mean_adjacent_coordinate += (*normal_axis == 0   ? barycenter.x()
                                         : *normal_axis == 1 ? barycenter.y()
                                                             : barycenter.z());
        }
        mean_adjacent_coordinate /= static_cast<double>(adjacent_indices.size());

        const double contact_plane = 0.5 * (region_min[*normal_axis] + region_max[*normal_axis]);
        if (mean_adjacent_coordinate < contact_plane) {
            box_min[*normal_axis] = contact_plane - contact_collection_depth;
            box_max[*normal_axis] = contact_plane + contact_margin;
        } else {
            box_min[*normal_axis] = contact_plane - contact_margin;
            box_max[*normal_axis] = contact_plane + contact_collection_depth;
        }

        const mesh::vector3 corner1{box_min[0], box_min[1], box_min[2]};
        const mesh::vector3 corner2{box_max[0], box_max[1], box_max[2]};
        simulation_device.add_contact(contact_name, corner1, corner2, ohmic_resistance);
        fmt::print("Added collecting contact '{}': ({:.6e}, {:.6e}, {:.6e}) to ({:.6e}, {:.6e}, {:.6e})\n",
                   contact_name,
                   corner1.x(),
                   corner1.y(),
                   corner1.z(),
                   corner2.x(),
                   corner2.y(),
                   corner2.z());
    }
}

}  // namespace uepm::PBMC
