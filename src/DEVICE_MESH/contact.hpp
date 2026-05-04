/**
 * @file contact.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2026-05-04
 * 
 * 
 */


#pragma once
#include "bbox.hpp"
#include "mesh.hpp"

namespace uepm::device {

/**
 * @brief This class represent a ohmic (metal-semiconductor) contact.
 * Its purpose is to be used in the main device class.
 * A contact is defined by a its name and by a box for its geometry.
 *
 */
class device_contact {
 protected:
    /**
     * @brief The name of the contact, e.g. "p-plus".
     *
     */
    std::string m_contact_name;

    /**
     * @brief The box that represent geometrically the contact.
     *
     */
    mesh::bbox m_contact_box;

    /**
     * @brief List of elements that define the contact.
     *
     */
    std::vector<mesh::element *> m_contact_elements_list;

    /**
     * @brief List of bulk elements that are adjacent to the contact.
     *
     */
    std::vector<mesh::element *> m_contact_bulk_adjacent_elements_list;

    /**
     * @brief The reference ohmic resistance of the contact.
     *
     */
    double m_ohmic_resistance;

    /**
     * @brief The current flowing through the contact.
     *
     */
    double m_contact_current = 0.0;

 public:
    /**
     * @brief Construct a new device contact object
     *
     * @param contact_name
     * @param contact_box
     */
    device_contact(const std::string &contact_name, mesh::vector3 corner1, mesh::vector3 corner2, double ohmic_resistance)
        : m_contact_name{contact_name},
          m_contact_box{corner1, corner2},
          m_ohmic_resistance{ohmic_resistance} {}

    /**
     * @brief Construct a new device contact object
     *
     * @param contact_name
     * @param contact_box
     */
    device_contact(const std::string &contact_name, const mesh::bbox &contact_box, double ohmic_resistance)
        : m_contact_name{contact_name},
          m_contact_box{contact_box},
          m_ohmic_resistance{ohmic_resistance} {}

    /**
     * @brief Set the element list object by getting the elements that intersect the contact box.
     *
     * @param p_mesh
     */
    void set_element_list(mesh::mesh *p_mesh) {
        for (auto &&p_elem : p_mesh->get_list_p_bulk_element()) {
            if (p_elem == nullptr) {
                continue;
            }
            if (m_contact_box.is_overlapping_tetra(*p_elem)) {
                m_contact_elements_list.push_back(p_elem);
            }
        }
        std::cout << "The contact is made of " << m_contact_elements_list.size() << " elements." << std::endl;
    }

    /**
     * @brief Return true if the segment [A, B] intersect the contact.
     *
     * @param point_A
     * @param point_B
     * @return true
     * @return false
     */
    bool line_intersect_contact(const mesh::vector3 &point_A, const mesh::vector3 &point_B) const {
        for (auto &&p_elem : m_contact_elements_list) {
            auto intersection = p_elem->compute_element_line_intersection(point_A, point_B);
            if (!intersection.empty()) {
                // std::cout << "Intersection contact at : " << intersection.value().m_intersection_location << std::endl;
                // std::cout << "with line: " << point_A << "  ,   " << point_B << std::endl;
                return true;
            }
        }
        return false;
    }

    bool is_element_in_contact(mesh::sp_element check_element) const {
        auto it = std::find_if(m_contact_elements_list.begin(), m_contact_elements_list.end(), [&](auto *p_element) {
            return p_element->has_same_vertices(*check_element);
        });
        return it != m_contact_elements_list.end();
    }

    bool point_in_contact_box(const mesh::vector3 &point) const { return m_contact_box.is_inside(point); }

    double get_ohmic_resistance() const { return m_ohmic_resistance; }
    void   set_ohmic_resistance(double ohmic_resistance) { m_ohmic_resistance = ohmic_resistance; }

    std::string get_contact_name() const { return m_contact_name; }
    void        set_contact_name(const std::string &contact_name) { m_contact_name = contact_name; }

    mesh::bbox get_contact_box() const { return m_contact_box; }
    void       set_contact_box(const mesh::bbox &contact_box) { m_contact_box = contact_box; }

    double get_contact_current() const { return m_contact_current; }
    void   set_contact_current(double contact_current) { m_contact_current = contact_current; }
    void   reset_contact_current() { m_contact_current = 0.0; }
    void   add_contact_current(double contact_current) { m_contact_current += contact_current; }
    void   get_contact_current(double contact_current, double mult_factor) { m_contact_current += contact_current * mult_factor; }
};

}  // namespace uepm::device