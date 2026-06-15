#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "bz_mesh.hpp"
#include "doctest/doctest.h"
#include "epm_material.hpp"

using uepm::mesh_bz::MeshBZ;
using uepm::mesh_bz::vector3;

namespace {

uepm::pseudopotential::epm_material make_test_material() {
    uepm::physics::material_info material_info{
        .id                 = uepm::physics::material_id::custom,
        .name               = "Test material",
        .symbol             = "test",
        .lattice_constant_m = 5.0e-10,
    };
    return {material_info, 0.0, 0.0, 0.0, 0.0};
}

}  // namespace

TEST_CASE("reduced and SI k-space conversions round trip") {
    const auto    material = make_test_material();
    const MeshBZ  mesh(material);
    const vector3 reduced(0.25, -0.5, 1.0);

    const vector3 si         = mesh.reduced_to_si_k(reduced);
    const vector3 round_trip = mesh.si_to_reduced_k(si);

    CHECK(round_trip.x() == doctest::Approx(reduced.x()));
    CHECK(round_trip.y() == doctest::Approx(reduced.y()));
    CHECK(round_trip.z() == doctest::Approx(reduced.z()));
}

TEST_CASE("reciprocal lattice vectors fold back into the first BZ") {
    const auto            material = make_test_material();
    MeshBZ                mesh(material);
    const double          reciprocal_scale = 1.0 / mesh.si_to_reduced_scale();
    const Eigen::Vector3d b1               = reciprocal_scale * Eigen::Vector3d{-1.0, 1.0, 1.0};
    const Eigen::Vector3d b2               = reciprocal_scale * Eigen::Vector3d{1.0, -1.0, 1.0};
    const Eigen::Vector3d b3               = reciprocal_scale * Eigen::Vector3d{1.0, 1.0, -1.0};
    mesh.init_reciprocal_basis(b1, b2, b3, 1.0, mesh.si_to_reduced_scale());

    const vector3 reciprocal_vector(b1.x() + b2.x(), b1.y() + b2.y(), b1.z() + b2.z());
    const vector3 folded = mesh.fold_ws_bcc(reciprocal_vector);

    CHECK(mesh.inside_ws_bcc(folded));
    CHECK(folded.norm() == doctest::Approx(0.0).scale(reciprocal_scale).epsilon(1e-12));
}
