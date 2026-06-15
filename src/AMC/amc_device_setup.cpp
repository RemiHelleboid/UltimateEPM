/**
 * @file amc_device_setup.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-09
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "amc_device_setup.hpp"

#include <fmt/core.h>

#include <filesystem>
#include <stdexcept>

namespace uepm::amc {

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
        throw std::invalid_argument("Only Si is currently supported by the analytical AMC transport model.");
    }
}

std::string make_default_output_directory(const std::string& mesh_file) {
    return fmt::format("self_consistent_amc_{}", std::filesystem::path(mesh_file).stem().string());
}

void add_default_contacts(uepm::device::device& simulation_device, uepm::mesh::mesh& mesh) {
    constexpr double contact_collection_depth = 0.001;  // µm
    constexpr double contact_margin           = 10.0;   // µm
    constexpr double ohmic_resistance         = 0.0;

    const double min_x = mesh.get_bounding_box().get_x_min();
    const double max_x = mesh.get_bounding_box().get_x_max();
    const double min_y = mesh.get_bounding_box().get_y_min();
    const double max_y = mesh.get_bounding_box().get_y_max();
    const double min_z = mesh.get_bounding_box().get_z_min();
    const double max_z = mesh.get_bounding_box().get_z_max();

    const mesh::vector3 cathode_corner1(min_x - contact_margin, min_y - contact_margin, min_z - contact_margin);
    const mesh::vector3 cathode_corner2(min_x + contact_collection_depth,
                                        max_y + contact_margin,
                                        max_z + contact_margin);

    simulation_device.add_contact("cathode", cathode_corner1, cathode_corner2, ohmic_resistance);
    const mesh::vector3 anode_corner1(max_x - contact_collection_depth, min_y - contact_margin, min_z - contact_margin);

    const mesh::vector3 anode_corner2(max_x + contact_margin, max_y + contact_margin, max_z + contact_margin);
    simulation_device.add_contact("anode", anode_corner1, anode_corner2, ohmic_resistance);

    fmt::print("Added device contacts:\n");
    fmt::print("  anode   x in [{:.6e}, {:.6e}]\n", anode_corner1.x(), anode_corner2.x());
    fmt::print("  cathode x in [{:.6e}, {:.6e}]\n", cathode_corner1.x(), cathode_corner2.x());
}

}  // namespace uepm::amc
