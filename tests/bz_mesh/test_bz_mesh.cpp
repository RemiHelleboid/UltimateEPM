#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <Eigen/LU>

#include "bz_mesh.hpp"
#include "band_catalog.hpp"
#include "bz_domain.hpp"
#include "doctest/doctest.h"
#include "epm_material.hpp"

using uepm::mesh_bz::MeshBZ;
using uepm::mesh_bz::MeshParticleType;
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

TEST_CASE("band catalog owns contiguous global and local band indexing") {
    uepm::mesh_bz::BandCatalog catalog;
    catalog.register_band(MeshParticleType::valence, -1.0, 0.0);
    catalog.register_band(MeshParticleType::valence, -2.0, -0.1);
    catalog.register_band(MeshParticleType::conduction, 1.1, 3.0);
    catalog.register_band(MeshParticleType::conduction, 1.3, 4.0);

    CHECK(catalog.total() == 4);
    CHECK(catalog.range(MeshParticleType::valence).global_start_index == 0);
    CHECK(catalog.range(MeshParticleType::valence).count == 2);
    CHECK(catalog.range(MeshParticleType::conduction).global_start_index == 2);
    CHECK(catalog.range(MeshParticleType::conduction).count == 2);
    CHECK(catalog.local_index(3) == 1);
    CHECK(catalog.global_index(1, MeshParticleType::conduction) == 3);
    CHECK(catalog.extrema(2).first == doctest::Approx(1.1));
}

TEST_CASE("band catalog rejects interleaved valence and conduction ranges") {
    uepm::mesh_bz::BandCatalog catalog;
    catalog.register_band(MeshParticleType::conduction, 1.0, 2.0);
    CHECK_THROWS_AS(catalog.register_band(MeshParticleType::valence, -1.0, 0.0), std::runtime_error);
}

TEST_CASE("positive-octant canonicalization preserves physical sign orientation") {
    const vector3 physical{-1.0, 2.0, -3.0};
    const auto canonical =
        uepm::mesh_bz::canonicalize_k(physical, uepm::mesh_bz::BZDomainMode::positive_octant);

    CHECK(canonical.representative.x() == doctest::Approx(1.0));
    CHECK(canonical.representative.y() == doctest::Approx(2.0));
    CHECK(canonical.representative.z() == doctest::Approx(3.0));
    CHECK(canonical.signs == std::array<int, 3>{-1, 1, -1});

    const vector3 reconstructed = uepm::mesh_bz::apply_sign_image(canonical.representative, canonical.signs);
    CHECK(reconstructed.x() == doctest::Approx(physical.x()));
    CHECK(reconstructed.y() == doctest::Approx(physical.y()));
    CHECK(reconstructed.z() == doctest::Approx(physical.z()));
}

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

TEST_CASE("Wigner-Seitz folding preserves momentum modulo the reciprocal lattice") {
    const auto            material = make_test_material();
    MeshBZ                mesh(material);
    const double          reciprocal_scale = 1.0 / mesh.si_to_reduced_scale();
    const Eigen::Vector3d b1               = reciprocal_scale * Eigen::Vector3d{-1.0, 1.0, 1.0};
    const Eigen::Vector3d b2               = reciprocal_scale * Eigen::Vector3d{1.0, -1.0, 1.0};
    const Eigen::Vector3d b3               = reciprocal_scale * Eigen::Vector3d{1.0, 1.0, -1.0};
    Eigen::Matrix3d       reciprocal_basis;
    reciprocal_basis.col(0) = b1;
    reciprocal_basis.col(1) = b2;
    reciprocal_basis.col(2) = b3;
    mesh.init_reciprocal_basis(b1, b2, b3, 1.0, mesh.si_to_reduced_scale());

    const std::array<vector3, 4> inputs = {
        vector3{0.9 * reciprocal_scale, 0.8 * reciprocal_scale, 0.7 * reciprocal_scale},
        vector3{-2.4 * reciprocal_scale, 1.3 * reciprocal_scale, 3.1 * reciprocal_scale},
        vector3{4.25 * reciprocal_scale, -3.75 * reciprocal_scale, 1.1 * reciprocal_scale},
        vector3{-7.2 * reciprocal_scale, 5.6 * reciprocal_scale, -4.4 * reciprocal_scale},
    };

    for (const vector3& input : inputs) {
        const vector3 folded = mesh.fold_ws_bcc(input);
        REQUIRE(mesh.inside_ws_bcc(folded));

        const Eigen::Vector3d delta(input.x() - folded.x(), input.y() - folded.y(), input.z() - folded.z());
        const Eigen::Vector3d lattice_coordinates = reciprocal_basis.inverse() * delta;
        for (Eigen::Index i = 0; i < lattice_coordinates.size(); ++i) {
            CHECK(lattice_coordinates[i] ==
                  doctest::Approx(std::round(lattice_coordinates[i])).epsilon(1e-11).scale(1.0));
        }
    }
}
