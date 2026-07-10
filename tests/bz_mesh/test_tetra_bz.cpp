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
#include <memory>
#include <numbers>
#include <random>
#include <stdexcept>
#include <vector>

#include "doctest/doctest.h"
#include "elph_deformation_potential.hpp"
#include "mesh_tetra.hpp"
#include "mesh_vertex.hpp"
#include "physical_constants.hpp"
#include "tetra_energy_index.hpp"
#include "vector_bz.hpp"

using uepm::mesh_bz::Tetra;
using uepm::mesh_bz::vector3;
using uepm::mesh_bz::Vertex;

TEST_CASE("YAML deformation-potential model uses initial carrier energy and phonon mode") {
    const vector3 q(3.0, 4.0, 0.0);

    const uepm::mesh_bz::DeformationPotential acoustic(uepm::mesh_bz::PhononMode::acoustic, 4.0, 5.0, 2.0);
    CHECK(acoustic.get_deformation_potential(q, 1.0) == doctest::Approx(15.0));
    CHECK(acoustic.get_deformation_potential(q, 3.0) == doctest::Approx(5.0 * std::sqrt(14.0)));

    const uepm::mesh_bz::DeformationPotential optical(uepm::mesh_bz::PhononMode::optical, 9.0, 7.0, 0.0);
    CHECK(optical.get_deformation_potential(q, 4.0) == doctest::Approx(3.0));
}

TEST_CASE("deformation-potential model rejects invalid configured values") {
    const vector3                             q(1.0, 0.0, 0.0);
    const uepm::mesh_bz::DeformationPotential invalid(uepm::mesh_bz::PhononMode::optical, 1.0, -2.0, 1.0);

    CHECK_THROWS_AS(invalid.get_deformation_potential(q, 1.0), std::domain_error);
    CHECK_THROWS_AS(invalid.get_deformation_potential(q, std::numeric_limits<double>::infinity()),
                    std::invalid_argument);
}

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

    constexpr double energy        = 0.25;
    const double     expected_area = 0.5 * (1.0 - energy) * (1.0 - energy);
    const double     expected_dos  = expected_area / (8.0 * std::numbers::pi * std::numbers::pi * std::numbers::pi);

    CHECK(tetra.compute_tetra_dos_energy_band(energy, 0) == doctest::Approx(expected_dos).epsilon(1e-12).scale(1.0));
}

TEST_CASE("reusing an iso-energy polygon preserves tetra DOS") {
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

    for (int sample = 1; sample < 20; ++sample) {
        const double energy  = static_cast<double>(sample) / 20.0;
        const auto   polygon = tetra.compute_band_iso_energy_polygon(energy, 0);
        const double direct  = tetra.compute_tetra_dos_energy_band(energy, 0);
        const double reused  = tetra.compute_tetra_dos_energy_band(energy, 0, polygon);
        CHECK(reused == doctest::Approx(direct).epsilon(1e-12));
    }
}

TEST_CASE("tetra energy interval index exactly matches brute-force overlap") {
    constexpr std::array<std::array<double, 4>, 4> energy_ranges = {
        std::array<double, 4>{0.0, 0.2, 0.4, 0.6},
        std::array<double, 4>{0.5, 0.8, 1.0, 1.2},
        std::array<double, 4>{1.1, 1.4, 1.8, 2.0},
        std::array<double, 4>{2.5, 2.7, 2.9, 3.0},
    };
    std::array<std::array<Vertex, 4>, energy_ranges.size()> vertices;
    std::vector<Tetra>                                      tetrahedra;
    tetrahedra.reserve(energy_ranges.size());

    for (std::size_t tetra_index = 0; tetra_index < energy_ranges.size(); ++tetra_index) {
        vertices[tetra_index] = {
            Vertex(4 * tetra_index + 0, v3(0.0, 0.0, 0.0)),
            Vertex(4 * tetra_index + 1, v3(1.0, 0.0, 0.0)),
            Vertex(4 * tetra_index + 2, v3(0.0, 1.0, 0.0)),
            Vertex(4 * tetra_index + 3, v3(0.0, 0.0, 1.0)),
        };
        std::array<Vertex*, 4> pointers{};
        for (std::size_t vertex_index = 0; vertex_index < 4; ++vertex_index) {
            vertices[tetra_index][vertex_index].add_band_energy_value(energy_ranges[tetra_index][vertex_index]);
            pointers[vertex_index] = &vertices[tetra_index][vertex_index];
        }
        tetrahedra.emplace_back(tetra_index, pointers);
        tetrahedra.back().compute_min_max_energies_at_bands();
    }

    uepm::mesh_bz::TetraEnergyIndex index;
    index.rebuild(tetrahedra, 1, 10.0, 1);

    for (const auto [minimum_energy, maximum_energy] :
         {std::pair{0.1, 0.3}, std::pair{0.55, 1.15}, std::pair{1.5, 2.6}, std::pair{3.1, 4.0}}) {
        std::vector<std::size_t> indexed;
        index.at(0).for_each_candidate(minimum_energy, maximum_energy, [&](std::size_t tetra_index) {
            indexed.push_back(tetra_index);
        });

        std::vector<std::size_t> brute_force;
        for (std::size_t tetra_index = 0; tetra_index < tetrahedra.size(); ++tetra_index) {
            if (tetrahedra[tetra_index].get_min_energy_at_band(0) <= maximum_energy &&
                tetrahedra[tetra_index].get_max_energy_at_band(0) >= minimum_energy) {
                brute_force.push_back(tetra_index);
            }
        }
        std::sort(indexed.begin(), indexed.end());
        std::sort(brute_force.begin(), brute_force.end());
        CHECK(indexed == brute_force);
    }
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

TEST_CASE("allocation-free tetra DOS matches the reference geometry implementation") {
    std::mt19937_64                        rng(0x5eed1234ULL);
    std::uniform_real_distribution<double> coordinate(-2.0, 2.0);
    std::uniform_real_distribution<double> energy_distribution(-3.0, 3.0);

    constexpr std::size_t tetra_count       = 250;
    constexpr std::size_t samples_per_tetra = 41;
    for (std::size_t tetra_index = 0; tetra_index < tetra_count; ++tetra_index) {
        std::array<vector3, 4>          positions;
        TetraWithVerts*                 fixture = nullptr;
        std::unique_ptr<TetraWithVerts> storage;
        do {
            positions = {
                v3(coordinate(rng), coordinate(rng), coordinate(rng)),
                v3(coordinate(rng), coordinate(rng), coordinate(rng)),
                v3(coordinate(rng), coordinate(rng), coordinate(rng)),
                v3(coordinate(rng), coordinate(rng), coordinate(rng)),
            };
            storage = std::make_unique<TetraWithVerts>(positions);
            fixture = storage.get();
        } while (std::abs(fixture->tet.get_signed_volume()) < 1e-3);

        std::array<double, 4> energies;
        for (std::size_t vertex_index = 0; vertex_index < energies.size(); ++vertex_index) {
            energies[vertex_index] = energy_distribution(rng);
            fixture->verts[vertex_index].add_band_energy_value(energies[vertex_index]);
        }
        fixture->tet.compute_min_max_energies_at_bands();

        const auto   minmax = std::minmax_element(energies.begin(), energies.end());
        const double span   = *minmax.second - *minmax.first;
        for (std::size_t sample = 0; sample < samples_per_tetra; ++sample) {
            const double fraction = static_cast<double>(sample) / static_cast<double>(samples_per_tetra - 1);
            const double energy   = *minmax.first + fraction * span;
            const double expected = fixture->tet.compute_tetra_dos_energy_band_reference(energy, 0);
            const double actual   = fixture->tet.compute_tetra_dos_energy_band(energy, 0);

            CHECK(actual == doctest::Approx(expected).epsilon(2e-12).scale(1e-30));
        }
    }
}

TEST_CASE("allocation-free tetra DOS preserves repeated-energy edge cases") {
    const std::array<std::array<double, 4>, 5> energy_cases = {
        std::array<double, 4>{0.0, 0.0, 1.0, 1.0},
        std::array<double, 4>{0.0, 0.0, 0.0, 1.0},
        std::array<double, 4>{0.0, 1.0, 1.0, 1.0},
        std::array<double, 4>{-1.0, 0.0, 0.0, 2.0},
        std::array<double, 4>{-1.0, -1.0, 2.0, 2.0},
    };

    for (const auto& energies : energy_cases) {
        std::array<Vertex, 4> vertices = {
            Vertex(0, v3(0.0, 0.0, 0.0)),
            Vertex(1, v3(1.0, 0.0, 0.0)),
            Vertex(2, v3(0.0, 1.0, 0.0)),
            Vertex(3, v3(0.0, 0.0, 1.0)),
        };
        for (std::size_t index = 0; index < vertices.size(); ++index) {
            vertices[index].add_band_energy_value(energies[index]);
        }
        std::array<Vertex*, 4> pointers = {&vertices[0], &vertices[1], &vertices[2], &vertices[3]};
        Tetra                  tetra(0, pointers);
        tetra.compute_min_max_energies_at_bands();

        const auto minmax = std::minmax_element(energies.begin(), energies.end());
        for (std::size_t sample = 0; sample <= 20; ++sample) {
            const double fraction = static_cast<double>(sample) / 20.0;
            const double energy   = *minmax.first + fraction * (*minmax.second - *minmax.first);
            CHECK(
                tetra.compute_tetra_dos_energy_band(energy, 0) ==
                doctest::Approx(tetra.compute_tetra_dos_energy_band_reference(energy, 0)).epsilon(2e-12).scale(1e-30));
        }
    }
}
