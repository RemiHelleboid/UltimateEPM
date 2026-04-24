/**
 * @file element1d.hpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief
 * @version 0.1
 * @date 2021-09-01
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <array>
#include <vector>

#include "element.hpp"

namespace uepm {

namespace mesh {

class element1d : public element {
 private:
 public:
    element1d() = delete;
    element1d(vertex *VtxA, vertex *VtxB) : element({VtxA, VtxB}) {}
    element1d(std::size_t index, vertex *VtxA, vertex *VtxB) : element(index, {VtxA, VtxB}) {}
    // Setters and Getters
    std::vector<std::size_t> get_vertices_index() const override { return {m_vertices[0]->get_index(), m_vertices[1]->get_index()}; }
    std::vector<std::array<std::size_t, 2>> get_edges_as_index_pair() const override;

    double get_length() const { return (*m_vertices[1] - *m_vertices[0]).norm(); }
    double compute_measure() const { return get_length(); }
    double get_measure() const override { return get_length(); }
    bool   is_location_inside_element(const vector3 &location) const override;

    std::vector<double> compute_barycentric_coordinate(const vector3 &location) const override;
    vector3             compute_gradient(const std::string &scalar_field_name) const override;

    double  integrate_scalar(const std::string &scalar_field_name) const override;
    vector3 integrate_vector(const std::string &vector_field_name) const override;


    vector3 compute_surface_normal() const override;
    virtual std::map<std::shared_ptr<element>, vector3> compute_element_line_intersection(const vector3 &point_A, const vector3 &point_B) const override;
    vector3 draw_uniform_random_point_inside_element() const override;
    vector3 draw_uniform_random_point_inside_element(std::minstd_rand &random_generator) const override;


};

}  // namespace mesh

}  // namespace uepm