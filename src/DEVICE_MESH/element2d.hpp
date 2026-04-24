/**
 * @file element2d.hpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief
 * @version 0.1
 * @date 2021-09-01
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once

#include <cmath>
#include <optional>
#include <vector>

#include "bbox.hpp"
#include "element.hpp"
#include "element1d.hpp"
#include "vector.hpp"
#include "vertex.hpp"

namespace uepm {

namespace mesh {

class element2d : public element {
 private:
    double        m_measure = 0.0;
    const vector3 m_precomputed_v0;
    const vector3 m_precomputed_v1;
    double        m_bary_coord_precomputed_v00     = 0.0;
    double        m_bary_coord_precomputed_v01     = 0.0;
    double        m_bary_coord_precomputed_v11     = 0.0;
    double        m_inverse_bary_coord_precomputed_divisor = 0.0;

 public:
    element2d() = default;
    element2d(vertex *VtxA, vertex *VtxB, vertex *VtxC)
        : element(0, {VtxA, VtxB, VtxC}),
          m_precomputed_v0{*m_vertices[1] - *m_vertices[0]},
          m_precomputed_v1{*m_vertices[2] - *m_vertices[0]},
          m_measure{compute_area()} {
        compute_precomputed_field_for_barycentric_coordinates();
    }

    element2d(std::size_t index, vertex *VtxA, vertex *VtxB, vertex *VtxC)
        : element(index, {VtxA, VtxB, VtxC}),
          m_precomputed_v0{*m_vertices[1] - *m_vertices[0]},
          m_precomputed_v1{*m_vertices[2] - *m_vertices[0]},
          m_measure{compute_area()} {
        compute_precomputed_field_for_barycentric_coordinates();
    }

    void compute_precomputed_field_for_barycentric_coordinates();

    // Setters and Getters
    std::vector<std::size_t> get_vertices_index() const override {
        return {m_vertices[0]->get_index(), m_vertices[1]->get_index(), m_vertices[2]->get_index()};
    }
    std::vector<std::pair<vertex *, vertex *>> get_edges_as_vertex_pair() const;
    std::vector<std::array<std::size_t, 2>>    get_edges_as_index_pair() const override;
    std::vector<element1d>                     get_list_edges() const;

    vector3 compute_surface_normal() const override;
    double  get_perimeters() const;
    double  compute_area() const;
    double  get_area() const { return m_measure; }
    double  compute_measure() const { return m_measure; }
    double  get_measure() const override { return m_measure; }
    bool    is_location_inside_element(const vector3 &location) const override;

    std::vector<double> compute_barycentric_coordinate(const vector3 &location) const override;
    vector3             compute_gradient(const std::string &scalar_field_name) const override;

    double  integrate_scalar(const std::string &scalar_field_name) const override;
    vector3 integrate_vector(const std::string &vector_field_name) const override;

    std::map<std::shared_ptr<element>, vector3> compute_element_line_intersection(const vector3 &point_A,
                                                                                  const vector3 &point_B) const override;
    std::optional<vector3>                      compute_line_triangle_intersection_3d(const vector3 &point_A, const vector3 &point_B) const;

    vector3 draw_uniform_random_point_inside_element() const override;
    vector3 draw_uniform_random_point_inside_element(std::minstd_rand &random_generator) const override;
};

}  // namespace mesh

}  // namespace uepm