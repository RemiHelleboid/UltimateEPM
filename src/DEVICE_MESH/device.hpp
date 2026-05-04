/**
 * @file device.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-05-04
 * 
 * 
 */

/**
 * @file device.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-11-15
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once

#include <algorithm>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "contact.hpp"
#include "mesh.hpp"

namespace uepm::device {

struct device_metadata {
    std::string m_name;
    std::string m_technology;
    std::string m_node_number;
    double      m_voltage_bias = std::numeric_limits<double>::quiet_NaN();

    device_metadata() = default;
    explicit device_metadata(const std::string &device_file_name);
};

class device {
 protected:
    mesh::mesh                 *m_mesh;
    std::vector<device_contact> m_list_contacts;
    device_metadata             m_metadata;
    double                      m_lattice_temperature = 300.0;

 public:
    explicit device(mesh::mesh *device_mesh) : m_mesh(device_mesh) {}
    device(mesh::mesh *device_mesh, const std::string &filename) : m_mesh(device_mesh), m_metadata(filename) {}

    void add_contact(const std::string &contact_name, const mesh::bbox contact_box, double ohmic_resistance = 0.0);
    void add_contact(const std::string &contact_name, mesh::vector3 corner1, mesh::vector3 corner2, double ohmic_resistance = 0.0);

    /**
     * @brief Return true if the segment [point_A, point_B] is crossing one of the device contact.
     *
     * @param point_A
     * @param point_B
     * @return true
     * @return false
     */
    bool check_crossing_contact(const mesh::vector3 &point_A, const mesh::vector3 &point_B) const;

    bool check_enters_contact(const mesh::vector3 &point);
    bool check_enters_contact(mesh::sp_element crossing_element) const;

    mesh::mesh            *get_p_mesh() const { return m_mesh; }
    int                    get_dimension() const { return m_mesh->get_dimension(); }
    const device_metadata &get_metadata() const { return m_metadata; }

    double interpolate_scalar_at_location(const std::string &fieldname, const mesh::vector3 &location) const {
        return m_mesh->interpolate_scalar_at_location(fieldname, location);
    }
    mesh::vector3 interpolate_vector_at_location(const std::string &fieldname, const mesh::vector3 &location) const {
        return m_mesh->interpolate_vector_at_location(fieldname, location);
    }
    mesh::element *find_element_at_location(const mesh::vector3 &location) const { return m_mesh->find_element_at_location(location); }
    const mesh::region_bulk *get_p_region_at_location(const mesh::vector3 &location) { return m_mesh->get_p_region_at_location(location); }
    std::string get_material_name_at_location(const mesh::vector3 &location) { return m_mesh->get_material_name_at_location(location); }
    std::string get_material_name_at_element(mesh::element *query_element) { return m_mesh->get_material_name_at_element(query_element); }

    double get_lattice_temperature() const { return m_lattice_temperature; }
    void   set_lattice_temperature(double new_temperature) { m_lattice_temperature = new_temperature; }

    std::vector<device_contact> get_list_contacts() const { return m_list_contacts; }

    std::vector<double> get_electrode_currents() const {
        std::vector<double> electrode_currents;
        for (auto &&contact : m_list_contacts) {
            electrode_currents.push_back(contact.get_contact_current());
        }
        return electrode_currents;
    }
    void reset_contacts_current() {
        for (auto &&contact : m_list_contacts) {
            contact.reset_contact_current();
        }
    }
};

}  // namespace uepm::device