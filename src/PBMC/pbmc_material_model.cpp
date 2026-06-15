/**
 * @file pbmc_material_model.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-26
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "pbmc_material_model.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

#include "physical_constants.hpp"
#include "yaml-cpp/yaml.h"

namespace uepm::PBMC {
namespace {

template <typename T>
T required_value(const YAML::Node& node, const char* key, const std::string& filename) {
    const auto value = node[key];
    if (!value) {
        throw std::runtime_error("PBMC material file '" + filename + "' does not define '" + key + "'.");
    }
    return value.as<T>();
}

YAML::Node required_map(const YAML::Node& node, const char* key, const std::string& filename) {
    const auto value = node[key];
    if (!value || !value.IsMap()) {
        throw std::runtime_error("PBMC material file '" + filename + "' must define map '" + key + "'.");
    }
    return value;
}

YAML::Node required_sequence(const YAML::Node& node, const char* key, const std::string& filename) {
    const auto value = node[key];
    if (!value || !value.IsSequence()) {
        throw std::runtime_error("PBMC material file '" + filename + "' must define sequence '" + key + "'.");
    }
    return value;
}

valley_model::band_type parse_dispersion(const std::string& value, const std::string& filename) {
    if (value == "parabolic") {
        return valley_model::band_type::parabolic;
    }
    if (value == "kane") {
        return valley_model::band_type::kane;
    }
    throw std::runtime_error("PBMC material file '" + filename + "' has invalid dispersion '" + value + "'.");
}

valley_model::mat3 parse_rotation(const YAML::Node& node, const std::string& filename) {
    if (!node || !node.IsSequence() || node.size() != 3) {
        throw std::runtime_error("PBMC material file '" + filename + "' contains an invalid rotation matrix.");
    }

    valley_model::mat3 rotation{};
    for (std::size_t row = 0; row < 3; ++row) {
        if (!node[row].IsSequence() || node[row].size() != 3) {
            throw std::runtime_error("PBMC material file '" + filename + "' contains an invalid rotation matrix.");
        }
        for (std::size_t column = 0; column < 3; ++column) {
            rotation[row][column] = node[row][column].as<double>();
        }
    }
    return rotation;
}

std::vector<valley_model> parse_bands(const YAML::Node& node, const std::string& filename) {
    std::vector<valley_model> bands;
    bands.reserve(node.size());
    for (const auto& entry : node) {
        valley_model::parameters parameters;
        parameters.name = required_value<std::string>(entry, "name", filename);
        parameters.transverse_effective_mass =
            required_value<double>(entry, "transverse_effective_mass_m0", filename) * uepm::constants::m_e;
        parameters.longitudinal_effective_mass =
            required_value<double>(entry, "longitudinal_effective_mass_m0", filename) * uepm::constants::m_e;
        parameters.non_parabolicity        = required_value<double>(entry, "non_parabolicity_eV_1", filename);
        parameters.energy_offset           = required_value<double>(entry, "energy_offset_eV", filename);
        parameters.phonon_reference_energy = required_value<double>(entry, "phonon_reference_energy_eV", filename);
        parameters.degeneracy              = required_value<std::size_t>(entry, "degeneracy", filename);
        parameters.dispersion = parse_dispersion(required_value<std::string>(entry, "dispersion", filename), filename);
        parameters.rotation   = parse_rotation(entry["rotation"], filename);
        bands.emplace_back(parameters);
    }
    return bands;
}

acoustic_scattering_parameters parse_acoustic_parameters(const YAML::Node&  node,
                                                         double             mass_density_kg_m3,
                                                         const std::string& filename) {
    return {
        .mass_density_kg_per_m3   = mass_density_kg_m3,
        .sound_velocity_m_per_s   = required_value<double>(node, "sound_velocity_m_per_s", filename),
        .deformation_potential_eV = required_value<double>(node, "deformation_potential_eV", filename),
        .overlap_factor           = required_value<double>(node, "overlap_factor", filename),
    };
}

intervalley_family parse_intervalley_family(const std::string& value, const std::string& filename) {
    if (value == "f") {
        return intervalley_family::f;
    }
    if (value == "g") {
        return intervalley_family::g;
    }
    throw std::runtime_error("PBMC material file '" + filename + "' has invalid intervalley family '" + value + "'.");
}

intervalley_order parse_intervalley_order(const std::string& value, const std::string& filename) {
    if (value == "zeroth") {
        return intervalley_order::zeroth;
    }
    if (value == "first") {
        return intervalley_order::first;
    }
    throw std::runtime_error("PBMC material file '" + filename + "' has invalid intervalley order '" + value + "'.");
}

std::vector<intervalley_phonon_branch> parse_intervalley_phonons(const YAML::Node& node, const std::string& filename) {
    std::vector<intervalley_phonon_branch> branches;
    branches.reserve(node.size());
    for (const auto& entry : node) {
        branches.push_back({
            .m_name   = required_value<std::string>(entry, "name", filename),
            .m_family = parse_intervalley_family(required_value<std::string>(entry, "family", filename), filename),
            .m_order  = parse_intervalley_order(required_value<std::string>(entry, "order", filename), filename),
            .m_phonon_energy_eV = required_value<double>(entry, "phonon_energy_eV", filename),
            .m_deformation_potential_0 =
                required_value<double>(entry, "zeroth_order_deformation_potential_eV_per_m", filename),
            .m_deformation_potential_1 =
                required_value<double>(entry, "first_order_deformation_potential_eV", filename),
            .m_final_valley_count = required_value<std::size_t>(entry, "final_valley_count", filename),
        });
    }
    return branches;
}

std::vector<hole_optical_transition> parse_hole_optical_phonons(const YAML::Node& node, const std::string& filename) {
    std::vector<hole_optical_transition> transitions;
    transitions.reserve(node.size());
    for (const auto& entry : node) {
        transitions.push_back({
            .name                           = required_value<std::string>(entry, "name", filename),
            .initial_band                   = required_value<std::size_t>(entry, "initial_band", filename),
            .final_band                     = required_value<std::size_t>(entry, "final_band", filename),
            .phonon_energy_eV               = required_value<double>(entry, "phonon_energy_eV", filename),
            .deformation_potential_eV_per_m = required_value<double>(entry, "deformation_potential_eV_per_m", filename),
            .overlap_factor                 = required_value<double>(entry, "overlap_factor", filename),
        });
    }
    return transitions;
}

impurity_mobility_parameters parse_impurity_mobility(const YAML::Node& node, const std::string& filename) {
    return {
        .m_mu0_cm2_per_V_s    = required_value<double>(node, "mu0_cm2_per_V_s", filename),
        .m_mu_min_cm2_per_V_s = required_value<double>(node, "mu_min_cm2_per_V_s", filename),
        .m_n_ref_cm_3         = required_value<double>(node, "reference_density_cm_3", filename),
        .m_alpha              = required_value<double>(node, "alpha", filename),
    };
}

impact_ionization_parameters parse_impact_ionization(const YAML::Node& node, const std::string& filename) {
    return {
        .m_threshold_eV  = required_value<double>(node, "threshold_eV", filename),
        .m_prefactor_s_1 = required_value<double>(node, "prefactor_s_1", filename),
        .m_exponent      = required_value<double>(node, "exponent", filename),
    };
}

}  // namespace

dielectric_properties make_dielectric_properties(const uepm::physics::material_info& material) {
    if (material.static_relative_permittivity <= 0.0) {
        throw std::invalid_argument("PBMC material requires a positive common static relative permittivity.");
    }
    return {.epsilon_r = material.static_relative_permittivity};
}

void pbmc_material_model::validate() const {
    if (m_dielectric.epsilon_r <= 0.0) {
        throw std::invalid_argument("PBMC material relative permittivity must be positive");
    }
    const auto validate_acoustic = [](const acoustic_scattering_parameters& parameters) {
        if (parameters.mass_density_kg_per_m3 <= 0.0) {
            throw std::invalid_argument("PBMC material mass density must be positive");
        }
        if (parameters.sound_velocity_m_per_s <= 0.0) {
            throw std::invalid_argument("PBMC material sound velocity must be positive");
        }
        if (parameters.deformation_potential_eV < 0.0) {
            throw std::invalid_argument("PBMC acoustic deformation potential must be non-negative");
        }
        if (parameters.overlap_factor < 0.0) {
            throw std::invalid_argument("PBMC acoustic overlap factor must be non-negative");
        }
    };
    const auto validate_mobility = [](const impurity_mobility_parameters& parameters) {
        if (parameters.m_mu0_cm2_per_V_s <= 0.0 || parameters.m_mu_min_cm2_per_V_s <= 0.0 ||
            parameters.m_n_ref_cm_3 <= 0.0 || parameters.m_alpha <= 0.0) {
            throw std::invalid_argument("PBMC impurity mobility parameters must be positive");
        }
    };
    const auto validate_impact = [](const impact_ionization_parameters& parameters) {
        if (parameters.m_threshold_eV <= 0.0 || parameters.m_prefactor_s_1 < 0.0 || parameters.m_exponent <= 0.0) {
            throw std::invalid_argument("PBMC impact ionization parameters are invalid");
        }
    };

    validate_acoustic(m_electron_acoustic);
    validate_acoustic(m_hole_acoustic);
    validate_mobility(m_impurity_mobility.m_electron);
    validate_mobility(m_impurity_mobility.m_hole);
    validate_impact(m_impact_ionization.m_electron);
    validate_impact(m_impact_ionization.m_hole);
    if (m_electron_valleys.empty() || m_hole_bands.empty()) {
        throw std::invalid_argument("PBMC material must define electron valleys and hole bands");
    }
    for (const auto& transition : m_electron_intervalley_transitions) {
        if (transition.m_name.empty() || transition.m_phonon_energy_eV <= 0.0 || transition.m_final_valley_count == 0) {
            throw std::invalid_argument("PBMC electron intervalley transition is invalid");
        }
        if (transition.is_zeroth_order() && transition.m_deformation_potential_0 <= 0.0) {
            throw std::invalid_argument("PBMC zeroth-order intervalley deformation potential must be positive");
        }
        if (transition.is_first_order() && transition.m_deformation_potential_1 <= 0.0) {
            throw std::invalid_argument("PBMC first-order intervalley deformation potential must be positive");
        }
    }
    for (const auto& transition : m_hole_optical_transitions) {
        if (transition.name.empty() || transition.phonon_energy_eV <= 0.0 ||
            transition.deformation_potential_eV_per_m <= 0.0 || transition.overlap_factor < 0.0) {
            throw std::invalid_argument("PBMC hole optical transition is invalid");
        }
        if (transition.initial_band >= m_hole_bands.size() || transition.final_band >= m_hole_bands.size()) {
            throw std::invalid_argument("PBMC hole optical transition references an invalid band");
        }
    }
}

pbmc_material_model load_pbmc_material_model(const uepm::physics::material_repository& repository,
                                             const uepm::physics::material_info&       common_material,
                                             const std::string&                        parameter_set) {
    if (common_material.id != uepm::physics::material_id::silicon || common_material.symbol != "Si") {
        throw std::invalid_argument("The analytical PBMC model currently requires common material 'Si'.");
    }
    if (common_material.mass_density_kg_m3 <= 0.0) {
        throw std::invalid_argument("The PBMC material requires a positive common mass density.");
    }

    const auto filename = repository.parameter_file(common_material.symbol, "pbmc", parameter_set);
    const auto config   = YAML::LoadFile(filename.string());
    if (!config.IsMap() || !config["schema_version"] || config["schema_version"].as<int>() != 1 ||
        !config["material"] || config["material"].as<std::string>() != common_material.symbol || !config["model"] ||
        config["model"].as<std::string>() != "pbmc" || !config["parameter_set"] ||
        config["parameter_set"].as<std::string>() != parameter_set) {
        throw std::runtime_error("Invalid PBMC material file '" + filename.string() + "'.");
    }

    const auto bands      = required_map(config, "bands", filename.string());
    const auto scattering = required_map(config, "scattering", filename.string());
    const auto acoustic   = required_map(scattering, "acoustic", filename.string());
    const auto mobility   = required_map(scattering, "impurity_mobility", filename.string());
    const auto impact     = required_map(scattering, "impact_ionization", filename.string());

    pbmc_material_model material;
    material.m_id                = common_material.id;
    material.m_dielectric        = make_dielectric_properties(common_material);
    material.m_electron_acoustic = parse_acoustic_parameters(required_map(acoustic, "electron", filename.string()),
                                                             common_material.mass_density_kg_m3,
                                                             filename.string());
    material.m_hole_acoustic     = parse_acoustic_parameters(required_map(acoustic, "hole", filename.string()),
                                                         common_material.mass_density_kg_m3,
                                                         filename.string());
    material.m_electron_valleys =
        parse_bands(required_sequence(bands, "electron_valleys", filename.string()), filename.string());
    material.m_hole_bands = parse_bands(required_sequence(bands, "hole_bands", filename.string()), filename.string());
    material.m_electron_intervalley_transitions =
        parse_intervalley_phonons(required_sequence(scattering, "electron_intervalley_phonons", filename.string()),
                                  filename.string());
    material.m_hole_optical_transitions =
        parse_hole_optical_phonons(required_sequence(scattering, "hole_optical_phonons", filename.string()),
                                   filename.string());
    material.m_impurity_mobility = {
        .m_electron = parse_impurity_mobility(required_map(mobility, "electron", filename.string()), filename.string()),
        .m_hole     = parse_impurity_mobility(required_map(mobility, "hole", filename.string()), filename.string()),
    };
    material.m_impact_ionization = {
        .m_electron = parse_impact_ionization(required_map(impact, "electron", filename.string()), filename.string()),
        .m_hole     = parse_impact_ionization(required_map(impact, "hole", filename.string()), filename.string()),
    };
    material.validate();
    return material;
}

pbmc_material_model load_pbmc_material_model(const uepm::physics::material_repository& repository,
                                             const std::string&                        material_symbol,
                                             const std::string&                        parameter_set) {
    return load_pbmc_material_model(repository, repository.load_material(material_symbol), parameter_set);
}

pbmc_material_model load_pbmc_material_model(const std::string& material_symbol, const std::string& parameter_set) {
    return load_pbmc_material_model(uepm::physics::material_repository{}, material_symbol, parameter_set);
}

carrier_impurity_mobility_parameters make_silicon_impurity_mobility_parameters() {
    return load_pbmc_material_model().m_impurity_mobility;
}

carrier_impact_ionization_parameters make_silicon_impact_ionization_parameters() {
    return load_pbmc_material_model().m_impact_ionization;
}

std::vector<valley_model> make_silicon_delta_valleys() { return load_pbmc_material_model().m_electron_valleys; }

std::vector<valley_model> make_silicon_hole_bands() { return load_pbmc_material_model().m_hole_bands; }

std::vector<hole_optical_transition> make_silicon_hole_optical_transitions() {
    return load_pbmc_material_model().m_hole_optical_transitions;
}

pbmc_material_model make_silicon_pbmc_material_model(const uepm::physics::material_info& common_material) {
    return load_pbmc_material_model(uepm::physics::material_repository{}, common_material);
}

pbmc_material_model make_silicon_pbmc_material_model() { return load_pbmc_material_model(); }

}  // namespace uepm::PBMC
