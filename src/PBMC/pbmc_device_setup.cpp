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

void add_collecting_contacts(uepm::device::device&          simulation_device,
                             uepm::mesh::mesh&              mesh,
                             const std::vector<std::string>& contact_names) {
    constexpr double contact_collection_depth = 0.001;  // µm
    constexpr double contact_margin           = 10.0;   // µm
    constexpr double ohmic_resistance         = 0.0;

    const mesh::bbox device_bbox = mesh.get_bounding_box();
    const std::array<double, 3> device_min{
        device_bbox.get_x_min(), device_bbox.get_y_min(), device_bbox.get_z_min()};
    const std::array<double, 3> device_max{
        device_bbox.get_x_max(), device_bbox.get_y_max(), device_bbox.get_z_max()};
    const double tolerance = 1.0e-8 * std::max(device_bbox.get_diagonal_size(), 1.0);

    for (const auto& contact_name : contact_names) {
        const auto* region = mesh.get_p_region(contact_name);
        if (region == nullptr) {
            throw std::runtime_error("Collecting contact region '" + contact_name + "' does not exist in the mesh.");
        }
        if (region->get_region_type() != mesh::RegionType::contact) {
            throw std::runtime_error("Collecting contact '" + contact_name + "' is not a mesh contact region.");
        }

        const mesh::bbox region_bbox = region->compute_bounding_box();
        const std::array<double, 3> region_min{
            region_bbox.get_x_min(), region_bbox.get_y_min(), region_bbox.get_z_min()};
        const std::array<double, 3> region_max{
            region_bbox.get_x_max(), region_bbox.get_y_max(), region_bbox.get_z_max()};
        std::array<double, 3> box_min{};
        std::array<double, 3> box_max{};
        bool                  found_boundary_normal = false;

        for (std::size_t axis = 0; axis < 3; ++axis) {
            const double device_size = device_max[axis] - device_min[axis];
            const bool at_min = std::abs(region_min[axis] - device_min[axis]) <= tolerance &&
                                std::abs(region_max[axis] - device_min[axis]) <= tolerance;
            const bool at_max = std::abs(region_min[axis] - device_max[axis]) <= tolerance &&
                                std::abs(region_max[axis] - device_max[axis]) <= tolerance;

            if (device_size > tolerance && at_min) {
                box_min[axis] = device_min[axis] - contact_margin;
                box_max[axis] = device_min[axis] + contact_collection_depth;
                found_boundary_normal = true;
            } else if (device_size > tolerance && at_max) {
                box_min[axis] = device_max[axis] - contact_collection_depth;
                box_max[axis] = device_max[axis] + contact_margin;
                found_boundary_normal = true;
            } else if (device_size <= tolerance) {
                box_min[axis] = region_min[axis] - contact_margin;
                box_max[axis] = region_max[axis] + contact_margin;
            } else if (std::abs(region_min[axis] - device_min[axis]) <= tolerance &&
                       std::abs(region_max[axis] - device_max[axis]) <= tolerance) {
                box_min[axis] = device_min[axis] - contact_margin;
                box_max[axis] = device_max[axis] + contact_margin;
            } else {
                box_min[axis] = region_min[axis] - tolerance;
                box_max[axis] = region_max[axis] + tolerance;
            }
        }

        if (!found_boundary_normal) {
            throw std::runtime_error("Collecting contact '" + contact_name +
                                     "' is not located on an exterior face of the mesh.");
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
