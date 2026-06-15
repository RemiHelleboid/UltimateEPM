#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "physical_constants.hpp"

namespace uepm::physics {

enum class material_id {
    custom,
    silicon,
    germanium,
    tin,
    gallium_phosphide,
    gallium_arsenide,
    aluminium_antimonide,
    aluminium_arsenide,
    indium_phosphide,
    indium_arsenide,
    indium_antimonide,
    zinc_selenide,
    cadmium_telluride,
    silicon_dioxide,
    gas
};

std::string_view to_string(material_id id);
material_id      material_id_from_name_or_symbol(std::string_view name_or_symbol);

struct material_info {
    material_id id = material_id::custom;

    std::string name;
    std::string symbol;

    double lattice_constant_m           = 0.0;
    double mass_density_kg_m3           = 0.0;
    double static_relative_permittivity = 0.0;

    double absolute_permittivity_F_m() const { return uepm::constants::eps_0 * static_relative_permittivity; }
};

struct material_state {
    material_id id = material_id::custom;

    double temperature_K        = 300.0;
    double donor_density_cm3    = 0.0;
    double acceptor_density_cm3 = 0.0;
};

class material_database {
 public:
    material_database() = default;
    explicit material_database(const std::string& filename) { load_from_file(filename); }

    void load_from_file(const std::string& filename);
    void add(material_info material);

    bool contains(material_id id) const;
    bool contains(const std::string& name_or_symbol) const;

    const material_info& require(material_id id) const;
    const material_info& require(const std::string& name_or_symbol) const;

    const std::vector<material_info>& materials() const noexcept { return m_materials; }
    void                              print_materials() const;

    // Transitional names retained for callers using the previous API.
    void load_materials_from_file(const std::string& filename) { load_from_file(filename); }
    bool is_material_available(const std::string& name_or_symbol) const { return contains(name_or_symbol); }
    const material_info& get_material(const std::string& name_or_symbol) const { return require(name_or_symbol); }

 private:
    std::vector<material_info> m_materials;
};

material_info silicon_material_info();

}  // namespace uepm::physics

namespace uepm::physic::material {

using material       = uepm::physics::material_info;
using list_materials = uepm::physics::material_database;

}  // namespace uepm::physic::material
