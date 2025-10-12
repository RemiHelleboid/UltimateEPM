/**
 * @file test_tetra_gradient.cpp
 * @brief Doctest for uepm::mesh_bz::Tetra::compute_gradient_scalar_field
 *
 * Assumptions:
 *  - Tetra constructor computes/stores m_signed_volume (or otherwise makes it available
 *    so compute_gradient_scalar_field can divide by 6*volume).
 *  - Vertex has a constructor taking (index, position) or adjust the factory below.
 *  - vector3 has a (x,y,z) constructor and x(), y(), z() getters (as used in your code).
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include <array>
#include <random>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

#include "mesh_tetra.hpp"
#include "mesh_vertex.hpp"

using uepm::mesh_bz::Tetra;
using uepm::mesh_bz::Vertex;
using uepm::mesh_bz::vector3;

// ---- Small helpers ----------------------------------------------------------

static inline vector3 v3(double x, double y, double z) {
    // If your vector3 uses another constructor API, adjust here.
    return vector3(x, y, z);
}

static inline double dot_xyz(double ax, double ay, double az, const vector3& p) {
    return ax * p.x() + ay * p.y() + az * p.z();
}

// Make a tetra and keep its Vertex storage alive in the calling scope.
struct TetraWithVerts {
    std::array<Vertex, 4> verts;
    std::array<Vertex*, 4> ptrs;
    Tetra tet;

    // NOTE: Adjust Vertex constructor if your API differs.
    TetraWithVerts(const std::array<vector3,4>& P)
        : verts{ Vertex(0, P[0]), Vertex(1, P[1]), Vertex(2, P[2]), Vertex(3, P[3]) }
        , ptrs{ &verts[0], &verts[1], &verts[2], &verts[3] }
        , tet(0, ptrs) {}
};

static inline std::array<double,4> values_from_linear_field(const Tetra& t, double ax, double ay, double az, double d) {
    const auto& V = t.get_list_vertices();
    std::array<double,4> vals{};
    for (int i=0;i<4;++i) {
        const auto& p = V[i]->get_position();
        vals[i] = dot_xyz(ax, ay, az, p) + d;
    }
    return vals;
}

static inline void CHECK_VEC3_CLOSE(const vector3& g, double gx, double gy, double gz, double eps=1e-12) {
    CHECK(g.x() == doctest::Approx(gx).epsilon(eps).scale(1.0));
    CHECK(g.y() == doctest::Approx(gy).epsilon(eps).scale(1.0));
    CHECK(g.z() == doctest::Approx(gz).epsilon(eps).scale(1.0));
}

