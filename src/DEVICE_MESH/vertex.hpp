#pragma once

#include <array>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "geometry_entity.hpp"
#include "vector_mesh.hpp"

namespace uepm {

namespace mesh {

class vertex : public geometry_entity, public vector3 {
 private:
    double  m_doping_concentration{0.0};
    double  m_space_charge{0.0};
    double  m_charge_density{0.0};
    double  m_electron_density{0.0};
    double  m_hole_density{0.0};
    vector3 m_electric_field{0.0, 0.0, 0.0};
    vector3 m_e_drift_velocity{0.0, 0.0, 0.0};
    vector3 m_h_drift_velocity{0.0, 0.0, 0.0};
    double  m_local_temperature{300.0};

 public:
    // Constructors
    vertex() = default;
    vertex(std::size_t index, double x, double y, double z) : geometry_entity{index}, vector3(x, y, z) {}

    // Electric field accessors
    inline vector3 get_electric_field() const noexcept { return m_electric_field; }
    inline void    set_electric_field(const vector3& new_electric_field) noexcept {
        m_electric_field = new_electric_field;
    }

    // Temperature accessors
    inline double get_temperature() const noexcept { return m_local_temperature; }
    inline void   set_temperature(double new_temperature) noexcept { m_local_temperature = new_temperature; }

    // Space charge accessors
    inline double get_space_charge() const noexcept { return m_space_charge; }
    inline void   set_space_charge(double new_space_charge) noexcept { m_space_charge = new_space_charge; }

    // Charge density accessors
    inline double get_charge_density() const noexcept { return m_charge_density; }
    inline void   set_charge_density(double new_density) noexcept { m_charge_density = new_density; }
    inline void   add_charge_to_vertex(double charge_quantity) noexcept { m_charge_density += charge_quantity; }
    inline double get_electron_density() const noexcept { return m_electron_density; }
    inline double get_hole_density() const noexcept { return m_hole_density; }
    inline void   set_mobile_carrier_densities(double electron_density, double hole_density) noexcept {
        m_electron_density = electron_density;
        m_hole_density     = hole_density;
        m_charge_density   = hole_density - electron_density;
    }

    // Doping concentration accessors
    inline double get_doping_concentration() const noexcept { return m_doping_concentration; }
    inline void   set_doping_concentration(double new_doping_concentration) noexcept {
        m_doping_concentration = new_doping_concentration;
    }

    // Friend functions to compute the distance between two vertices.
    friend inline double distance_between(const vertex& vtx_1, const vertex& vtx_2) noexcept {
        return (static_cast<const vector3&>(vtx_1) - static_cast<const vector3&>(vtx_2)).norm();
    }

    friend inline double distance_between(const vertex* vtx_1, const vertex* vtx_2) noexcept {
        return (static_cast<const vector3&>(*vtx_1) - static_cast<const vector3&>(*vtx_2)).norm();
    }

    // Function to print information about the vertex.
    void print_info() const;
};

}  // namespace mesh

}  // namespace uepm
