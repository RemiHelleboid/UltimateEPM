/**
 * @file amc_device_setup.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-09
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <string>
#include <string_view>

#include "amc_transport_kernel.hpp"
#include "particle_amc.hpp"
#include "device.hpp"
#include "mesh.hpp"

namespace uepm::amc {

impurity_scattering_model parse_impurity_model(const std::string& text);
impurity_screening_model  parse_impurity_screening_model(const std::string& text);
std::string_view          impurity_screening_model_name(impurity_screening_model model);
particle_type             parse_particle_type(const std::string& text);

void validate_material_symbol(const std::string& material_symbol);

std::string make_default_output_directory(const std::string& mesh_file);

void add_default_contacts(uepm::device::device& simulation_device, uepm::mesh::mesh& mesh);

}  // namespace uepm::amc
