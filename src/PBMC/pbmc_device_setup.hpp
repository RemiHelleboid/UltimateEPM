/**
 * @file pbmc_device_setup.hpp
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
#include <vector>

#include "device.hpp"
#include "mesh.hpp"
#include "pbmc_particle.hpp"
#include "pbmc_transport_kernel.hpp"

namespace uepm::PBMC {

impurity_scattering_model parse_impurity_model(const std::string& text);
impurity_screening_model  parse_impurity_screening_model(const std::string& text);
std::string_view          impurity_screening_model_name(impurity_screening_model model);
contact_injection_distribution parse_contact_injection_distribution(const std::string& text);
std::string_view contact_injection_distribution_name(contact_injection_distribution distribution);
particle_type             parse_particle_type(const std::string& text);

void validate_material_symbol(const std::string& material_symbol);

std::string make_default_output_directory(const std::string& mesh_file);

void add_collecting_contacts(uepm::device::device&           simulation_device,
                             uepm::mesh::mesh&               mesh,
                             const std::vector<std::string>& contact_names);

}  // namespace uepm::PBMC
