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
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>
#include <vector>

#include "doctest/doctest.h"
#include "mesh_tetra.hpp"
#include "mesh_vertex.hpp"

using uepm::mesh_bz::Tetra;
using uepm::mesh_bz::vector3;
using uepm::mesh_bz::Vertex;

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
    std::array<Vertex, 4>  verts;
    std::array<Vertex*, 4> ptrs;
    Tetra                  tet;

    // NOTE: Adjust Vertex constructor if your API differs.
    TetraWithVerts(const std::array<vector3, 4>& P)
        : verts{Vertex(0, P[0]), Vertex(1, P[1]), Vertex(2, P[2]), Vertex(3, P[3])},
          ptrs{&verts[0], &verts[1], &verts[2], &verts[3]},
          tet(0, ptrs) {}
};

static inline std::array<double, 4> values_from_linear_field(const Tetra& t,
                                                             double       ax,
                                                             double       ay,
                                                             double       az,
                                                             double       d) {
    const auto&           V = t.get_list_vertices();
    std::array<double, 4> vals{};
    for (int i = 0; i < 4; ++i) {
        const auto& p = V[i]->get_position();
        vals[i]       = dot_xyz(ax, ay, az, p) + d;
    }
    return vals;
}

static inline void CHECK_VEC3_CLOSE(const vector3& g, double gx, double gy, double gz, double eps = 1e-12) {
    CHECK(g.x() == doctest::Approx(gx).epsilon(eps).scale(1.0));
    CHECK(g.y() == doctest::Approx(gy).epsilon(eps).scale(1.0));
    CHECK(g.z() == doctest::Approx(gz).epsilon(eps).scale(1.0));
}

TEST_CASE("tetra location test includes faces and vertices") {
    TetraWithVerts fixture({v3(0.0, 0.0, 0.0), v3(1.0, 0.0, 0.0), v3(0.0, 1.0, 0.0), v3(0.0, 0.0, 1.0)});

    CHECK(fixture.tet.is_location_inside(v3(0.25, 0.25, 0.25)));
    CHECK(fixture.tet.is_location_inside(v3(0.0, 0.25, 0.25)));
    CHECK(fixture.tet.is_location_inside(v3(0.0, 0.0, 0.0)));
    CHECK_FALSE(fixture.tet.is_location_inside(v3(-1e-6, 0.25, 0.25)));
}

TEST_CASE("degenerate tetrahedra reject barycentric interpolation") {
    TetraWithVerts fixture({v3(0.0, 0.0, 0.0), v3(1.0, 0.0, 0.0), v3(0.0, 1.0, 0.0), v3(1.0, 1.0, 0.0)});

    CHECK_FALSE(fixture.tet.is_location_inside(v3(0.25, 0.25, 0.0)));
    CHECK_THROWS_AS(fixture.tet.compute_barycentric_coordinates(v3(0.25, 0.25, 0.0)), std::domain_error);
}

TEST_CASE("tetra DOS matches an analytic linear band") {
    std::array<Vertex, 4> vertices = {
        Vertex(0, v3(0.0, 0.0, 0.0)),
        Vertex(1, v3(1.0, 0.0, 0.0)),
        Vertex(2, v3(0.0, 1.0, 0.0)),
        Vertex(3, v3(0.0, 0.0, 1.0)),
    };
    for (Vertex& vertex : vertices) {
        vertex.add_band_energy_value(vertex.get_position().x());
    }
    std::array<Vertex*, 4> pointers = {&vertices[0], &vertices[1], &vertices[2], &vertices[3]};
    Tetra                  tetra(0, pointers);
    tetra.compute_min_max_energies_at_bands();

    constexpr double energy = 0.25;
    const double expected_area = 0.5 * (1.0 - energy) * (1.0 - energy);
    const double expected_dos =
        expected_area / (8.0 * std::numbers::pi * std::numbers::pi * std::numbers::pi);

    CHECK(tetra.compute_tetra_dos_energy_band(energy, 0) ==
          doctest::Approx(expected_dos).epsilon(1e-12).scale(1.0));
}

TEST_CASE("equal vertex energies produce a finite quadrilateral iso-surface") {
    std::array<Vertex, 4> vertices = {
        Vertex(0, v3(0.0, 0.0, 0.0)),
        Vertex(1, v3(1.0, 0.0, 0.0)),
        Vertex(2, v3(0.0, 1.0, 0.0)),
        Vertex(3, v3(0.0, 0.0, 1.0)),
    };
    const std::array<double, 4> energies = {0.0, 0.0, 1.0, 1.0};
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        vertices[i].add_band_energy_value(energies[i]);
    }
    std::array<Vertex*, 4> pointers = {&vertices[0], &vertices[1], &vertices[2], &vertices[3]};
    Tetra                  tetra(0, pointers);

    const auto surface = tetra.compute_band_iso_energy_surface(0.5, 0);
    REQUIRE(surface.size() == 4);
    for (const vector3& point : surface) {
        CHECK(std::isfinite(point.x()));
        CHECK(std::isfinite(point.y()));
        CHECK(std::isfinite(point.z()));
    }
}
