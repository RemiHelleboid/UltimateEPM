#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>

#include <random>

#include "boundary_reflection.hpp"
#include "element2d.hpp"
#include "element3d.hpp"
#include "pbmc_transport_kernel.hpp"
#include "vertex.hpp"

TEST_CASE("boundary reflection detects a 2D triangle exit edge and specularly reflects") {
    uepm::mesh::vertex          v0{0, 0.0, 0.0, 0.0};
    uepm::mesh::vertex          v1{1, 1.0, 0.0, 0.0};
    uepm::mesh::vertex          v2{2, 0.0, 1.0, 0.0};
    const uepm::mesh::element2d triangle{&v0, &v1, &v2};

    const uepm::mesh::vector3 previous{0.25, 0.25, 0.0};
    const uepm::mesh::vector3 trial{0.45, -0.25, 0.0};
    const auto                hit = uepm::mesh::find_boundary_exit_hit(triangle, previous, trial, 2);

    REQUIRE(hit.has_value());
    CHECK(hit->position.x() == doctest::Approx(0.35));
    CHECK(hit->position.y() == doctest::Approx(0.0));
    CHECK(hit->inward_normal.x() == doctest::Approx(0.0));
    CHECK(hit->inward_normal.y() == doctest::Approx(1.0));

    const uepm::mesh::vector3 reflected = uepm::mesh::reflect_vector_specular({2.0, -3.0, 0.0}, hit->inward_normal);
    CHECK(reflected.x() == doctest::Approx(2.0));
    CHECK(reflected.y() == doctest::Approx(3.0));

    const auto reflected_position = uepm::mesh::place_reflected_position_inside(
        triangle,
        previous,
        trial,
        *hit,
        uepm::mesh::reflect_vector_specular(trial - hit->position, hit->inward_normal),
        2);
    CHECK(reflected_position.x() > previous.x());
    CHECK(reflected_position.x() == doctest::Approx(0.45));
    CHECK(reflected_position.y() > 0.0);
    CHECK(reflected_position.y() == doctest::Approx(previous.y()));
    CHECK(triangle.is_location_inside_element(reflected_position));
}

TEST_CASE("boundary reflection detects a 3D tetrahedron exit face") {
    uepm::mesh::vertex          v0{0, 0.0, 0.0, 0.0};
    uepm::mesh::vertex          v1{1, 1.0, 0.0, 0.0};
    uepm::mesh::vertex          v2{2, 0.0, 1.0, 0.0};
    uepm::mesh::vertex          v3{3, 0.0, 0.0, 1.0};
    const uepm::mesh::element3d tetra{&v0, &v1, &v2, &v3};

    const uepm::mesh::vector3 previous{0.1, 0.1, 0.1};
    const uepm::mesh::vector3 trial{-0.1, 0.1, 0.1};
    const auto                hit = uepm::mesh::find_boundary_exit_hit(tetra, previous, trial, 3);

    REQUIRE(hit.has_value());
    CHECK(hit->position.x() == doctest::Approx(0.0));
    CHECK(hit->position.y() == doctest::Approx(0.1));
    CHECK(hit->position.z() == doctest::Approx(0.1));
    CHECK(hit->inward_normal.x() > 0.0);
    CHECK(hit->inward_normal.y() == doctest::Approx(0.0));
    CHECK(hit->inward_normal.z() == doctest::Approx(0.0));
}

TEST_CASE("diffuse reflection preserves magnitude and points inward") {
    std::minstd_rand          rng{7};
    const uepm::mesh::vector3 inward_normal{0.0, 1.0, 0.0};
    const uepm::mesh::vector3 incoming{3.0, -4.0, 0.0};

    const auto reflected = uepm::mesh::draw_diffuse_reflection_vector(incoming, inward_normal, 2, rng);

    CHECK(reflected.norm() == doctest::Approx(incoming.norm()));
    CHECK(reflected.dot(inward_normal) >= 0.0);
    CHECK(reflected.z() == doctest::Approx(0.0));
}

TEST_CASE("PBMC boundary direction update preserves kinetic energy") {
    uepm::PBMC::pbmc_transport_config config;
    config.m_carrier_type             = uepm::PBMC::particle_type::electron;
    config.m_lattice_temperature      = 300.0;
    config.m_enable_impact_ionization = false;

    uepm::PBMC::pbmc_transport_kernel transport{config, 17};
    transport.initialize();

    uepm::PBMC::pbmc_particle particle{0, uepm::PBMC::particle_type::electron};
    transport.initialize_particle_state(particle);

    const double initial_energy = particle.state().kinetic_energy;
    transport.set_particle_velocity_direction_preserving_energy(particle, {1.0, 2.0, 0.0});

    CHECK(particle.state().kinetic_energy == doctest::Approx(initial_energy).epsilon(1.0e-12));
    CHECK(particle.state().velocity.dot({1.0, 2.0, 0.0}) > 0.0);
}
