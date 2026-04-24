/**
 * @file element.hpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief Header file for implementation of element, element1d, element2d and element3d classes.
 * @version 0.1
 * @date 2021-07-10
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <array>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "bbox.hpp"
#include "geometry_entity.hpp"
#include "vertex.hpp"

namespace uepm {

namespace mesh {

/**
 * @brief Small struct to carry the result of "compute_element_line_intersection" function.
 *
 */
struct element_line_intersection {
    /**
     * @brief Position of the intersection
     *
     */
    vector3 m_intersection_location;

    /**
     * @brief Pointer on the element face on which the intersection was found.
     *
     */
    std::shared_ptr<element> m_p_element_face_intersection;

    element_line_intersection(const vector3 &intersection_location, std::shared_ptr<element> element_face_intersection)
        : m_intersection_location(intersection_location),
          m_p_element_face_intersection(element_face_intersection){};
};

enum class element_type { edge = 1, triangle = 2, rectange = 3, tetrahedron = 4 };

class element : public geometry_entity {
 protected:
    std::vector<vertex *> m_vertices;
    double                m_n_charge_density = 0.0;
    double                m_p_charge_density = 0.0;
    bbox                  m_bounding_box;
    int                   m_region_index = -1;

    /**
     * @brief Flag to indicate if the velocity diffusion is up to date, meaning that 
     * it has been computed at the current time step (so with the right electric field).
     * Used to avoid computing it multiple times.
     *
     */
    bool m_velocity_diffusion_up_to_date = false;

 public:
    // Constructors
    element() {}
    explicit element(const std::vector<vertex *> &vertices)
        : geometry_entity(0),
          m_vertices(vertices),
          m_bounding_box(compute_bounding_box()) {}
    element(std::size_t index_element, const std::vector<vertex *> &vertices)
        : geometry_entity(index_element),
          m_vertices(vertices),
          m_bounding_box(compute_bounding_box()) {}
    virtual ~element() {}
    // Setters and Getters
    void                                    set_region_index(unsigned int new_index) { m_region_index = new_index; }
    int                                     get_region_index() { return m_region_index; }
    const std::vector<vertex *>            &get_vertices() const { return m_vertices; }
    vertex                                 *get_vertex(const std::size_t index) const;
    virtual std::vector<std::size_t>        get_vertices_index() const;
    bool                                    contains_vertex_with_index(std::size_t index_vtx) const;
    std::optional<std::vector<std::size_t>> is_adjacent_to_element(const element &second_element) const;
    // Others
    virtual std::vector<std::array<std::size_t, 2>> get_edges_as_index_pair() const { return {}; }
    unsigned int                                    get_element_type() const;
    vector3                                         get_barycenter() const;
    virtual double                                  get_measure() const = 0;
    const bbox                                     &get_bounding_box() const { return m_bounding_box; }
    bbox                                            compute_bounding_box() const;
    virtual bool                                    is_location_inside_element(const vector3 &location) const = 0;
    virtual vector3                                 compute_surface_normal() const                            = 0;

    double get_n_charge() const { return m_n_charge_density; }
    void   set_n_charge(double new_density) { m_n_charge_density = new_density; }
    void   add_n_charge(double additional_charge) { m_n_charge_density += additional_charge; }
    double get_p_charge() const { return m_p_charge_density; }
    void   set_p_charge(double new_density) { m_p_charge_density = new_density; }
    void   add_p_charge(double additional_charge) { m_p_charge_density += additional_charge; }
    

    virtual std::vector<double> compute_barycentric_coordinate(const vector3 &location) const = 0;
    virtual double              interpolate_scalar_at_location(const std::string &name, const vector3 &location) const;
    virtual vector3             interpolate_vector_at_location(const std::string &name, const vector3 &location) const;
    virtual double              interpolate_doping_at_location(const vector3 &location) const;
    virtual vector3             interpolate_electric_field_at_location(const vector3 &location) const;
    virtual vector3             interpolate_e_grad_diffusion_at_location(const vector3 &location) const;
    virtual vector3             interpolate_h_grad_diffusion_at_location(const vector3 &location) const;
    virtual vector3             compute_gradient(const std::string &scalar_field_name) const                            = 0;
    virtual std::map<std::shared_ptr<element>, vector3> compute_element_line_intersection(const vector3 &point_A,
                                                                                          const vector3 &point_B) const = 0;
    virtual void                                        distribute_charge_quantity_over_vertices(double charge_quantity);

    virtual double  integrate_scalar(const std::string &scalar_field_name) const = 0;
    virtual vector3 integrate_vector(const std::string &vector_field_name) const = 0;

    bool has_same_vertices(const element &second_element) const {
        return std::is_permutation(m_vertices.begin(), m_vertices.end(), second_element.get_vertices().begin());
    }

    virtual vector3 draw_uniform_random_point_inside_element() const = 0;
    virtual vector3 draw_uniform_random_point_inside_element(std::minstd_rand &random_generator) const = 0;
};

}  // namespace mesh

}  // namespace uepm