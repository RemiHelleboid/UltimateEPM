/**
 * @file element3d.hpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief
 * @version 0.1
 * @date 2021-09-01
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <optional>
#include <vector>

#include "bbox.hpp"
#include "element.hpp"
#include "element1d.hpp"
#include "element2d.hpp"
#include "vector.hpp"
#include "vertex.hpp"

namespace uepm {

namespace mesh {

class bbox;

class element3d : public element {
 private:
    const vector3 m_v0;
    const vector3 m_v1;
    const vector3 m_v2;
    const vector3 m_v3;
    const vector3 m_v4;
    const double  m_signed_volume;

 public:
    element3d(vertex *VtxA, vertex *VtxB, vertex *VtxC, vertex *VtxD)
        : element(0, {VtxA, VtxB, VtxC, VtxD}),
          m_v0{*m_vertices[1] - *m_vertices[0]},
          m_v1{*m_vertices[2] - *m_vertices[0]},
          m_v2{*m_vertices[3] - *m_vertices[0]},
          m_v3{*m_vertices[2] - *m_vertices[1]},
          m_v4{*m_vertices[3] - *m_vertices[1]},
          m_signed_volume{(1.0 / 6.0) * scalar_triple_product(m_v0, m_v1, m_v2)} {}

    element3d(std::size_t index, vertex *VtxA, vertex *VtxB, vertex *VtxC, vertex *VtxD)
        : element(index, {VtxA, VtxB, VtxC, VtxD}),
          m_v0{*m_vertices[1] - *m_vertices[0]},
          m_v1{*m_vertices[2] - *m_vertices[0]},
          m_v2{*m_vertices[3] - *m_vertices[0]},
          m_v3{*m_vertices[2] - *m_vertices[1]},
          m_v4{*m_vertices[3] - *m_vertices[1]},
          m_signed_volume{(1.0 / 6.0) * scalar_triple_product(m_v0, m_v1, m_v2)} {}
    // Setters and Getters
    std::vector<std::size_t> get_vertices_index() const override {
        return {m_vertices[0]->get_index(), m_vertices[1]->get_index(), m_vertices[2]->get_index(), m_vertices[3]->get_index()};
    }
    std::vector<std::array<std::size_t, 2>> get_edges_as_index_pair() const override;
    std::vector<element1d>                  get_list_edges() const;
    std::vector<element2d>                  get_faces_as_element2d() const;

    vector3 compute_surface_normal() const override;

    double                   get_signed_volume() const { return m_signed_volume; }
    double                   get_volume() const { return fabs(m_signed_volume); }
    double                   compute_measure() const { return fabs(m_signed_volume); }
    double                   get_measure() const override { return fabs(get_signed_volume()); }
    bool                     is_location_inside_element(const vector3 &location) const override;
    std::optional<element2d> is_adjacent_to_element(const element3d &second_element) const;

    std::vector<double> compute_barycentric_coordinate(const vector3 &location) const override;
    vector3             compute_gradient(const std::string &scalar_field_name) const override;

    double  integrate_scalar(const std::string &scalar_field_name) const override;
    vector3 integrate_vector(const std::string &vector_field_name) const override;

    std::map<std::shared_ptr<element>, vector3> compute_element_line_intersection(const vector3 &point_A,
                                                                                  const vector3 &point_B) const override;

    vector3 draw_uniform_random_point_inside_element() const override;
    vector3 draw_uniform_random_point_inside_element(std::minstd_rand &random_generator) const override;
};

}  // namespace mesh

}  // namespace uepm