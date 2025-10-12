
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <array>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "bbox_mesh.hpp"
#include "doctest/doctest.h"

using uepm::mesh_bz::bbox_mesh;
using uepm::mesh_bz::vector3;

TEST_SUITE("[bbox_mesh] construction & invariants") {
    TEST_CASE("default ctor yields zero-sized box at origin") {
        bbox_mesh b;
        CHECK_EQ(b.get_x_min(), 0.0);
        CHECK_EQ(b.get_x_max(), 0.0);
        CHECK_EQ(b.get_y_min(), 0.0);
        CHECK_EQ(b.get_y_max(), 0.0);
        CHECK_EQ(b.get_z_min(), 0.0);
        CHECK_EQ(b.get_z_max(), 0.0);
        CHECK_EQ(b.get_volume(), 0.0);
        CHECK_EQ(b.get_center().x(), 0.0);
        CHECK_EQ(b.get_center().y(), 0.0);
        CHECK_EQ(b.get_center().z(), 0.0);
    }

    TEST_CASE("range ctor validates ordering") {
        // valid
        bbox_mesh ok(-1.0, 2.0, 0.0, 4.0, -3.0, -1.0);
        CHECK_EQ(ok.get_x_min(), -1.0);
        CHECK_EQ(ok.get_x_max(), 2.0);
        // invalid (x_min > x_max)
        CHECK_THROWS_AS(bbox_mesh(1.0, -1.0, 0.0, 1.0, 0.0, 1.0), std::invalid_argument);
        // invalid (y_min > y_max)
        CHECK_THROWS_AS(bbox_mesh(0.0, 1.0, 2.0, 1.0, 0.0, 1.0), std::invalid_argument);
        // invalid (z_min > z_max)
        CHECK_THROWS_AS(bbox_mesh(0.0, 1.0, 0.0, 1.0, 5.0, 4.0), std::invalid_argument);
    }

    TEST_CASE("corner ctor builds the same box") {
        vector3   a(-1.0, 0.0, -3.0);  // bottom-left-front
        vector3   b(2.0, 4.0, -1.0);   // up-right-back
        bbox_mesh bb(a, b);
        CHECK_EQ(bb.get_x_min(), -1.0);
        CHECK_EQ(bb.get_x_max(), 2.0);
        CHECK_EQ(bb.get_y_min(), 0.0);
        CHECK_EQ(bb.get_y_max(), 4.0);
        CHECK_EQ(bb.get_z_min(), -3.0);
        CHECK_EQ(bb.get_z_max(), -1.0);
    }
}

TEST_SUITE("[bbox_mesh] geometry & metrics") {
    TEST_CASE("sizes, diagonal, volume, center") {
        bbox_mesh b(0.0, 2.0, 10.0, 14.0, -8.0, 0.0);
        // True extents:
        const double sx = 2.0;  // x: 0..2
        const double sy = 4.0;  // y: 10..14
        const double sz = 8.0;  // z: -8..0

        // ---- NOTE: The class currently has y/z getters swapped. ----
        // get_x_size() is fine:
        CHECK_EQ(b.get_x_size(), doctest::Approx(sx));

        // The next two checks encode intended behavior and will FAIL with the current header.
        // Keep them as-is to catch the bug (y/z swapped). Remove if you fix the implementation.
        CHECK_EQ(b.get_y_size(), doctest::Approx(sy));  // intended: 4.0
        CHECK_EQ(b.get_z_size(), doctest::Approx(sz));  // intended: 8.0

        // Diagonal and volume are unaffected by the swap and should still match:
        CHECK_EQ(b.get_diagonal_size(), doctest::Approx(std::sqrt(sx * sx + sy * sy + sz * sz)));
        CHECK_EQ(b.get_volume(), doctest::Approx(sx * sy * sz));

        auto c = b.get_center();
        CHECK_EQ(c.x(), doctest::Approx(1.0));
        CHECK_EQ(c.y(), doctest::Approx(12.0));
        CHECK_EQ(c.z(), doctest::Approx(-4.0));
    }

    TEST_CASE("is_inside uses strict inequalities (boundary excluded)") {
        bbox_mesh b(0.0, 1.0, 0.0, 2.0, 0.0, 3.0);
        // interior
        CHECK(b.is_inside(vector3(0.5, 1.0, 1.5)));
        // on faces/edges/corners -> false
        CHECK_FALSE(b.is_inside(vector3(0.0, 1.0, 1.0)));  // x_min plane
        CHECK_FALSE(b.is_inside(vector3(1.0, 1.0, 1.0)));  // x_max plane
        CHECK_FALSE(b.is_inside(vector3(0.5, 0.0, 1.0)));  // y_min
        CHECK_FALSE(b.is_inside(vector3(0.5, 2.0, 1.0)));  // y_max
        CHECK_FALSE(b.is_inside(vector3(0.5, 1.0, 0.0)));  // z_min
        CHECK_FALSE(b.is_inside(vector3(0.5, 1.0, 3.0)));  // z_max
        CHECK_FALSE(b.is_inside(vector3(0.0, 0.0, 0.0)));  // corner
    }
}

TEST_SUITE("[bbox_mesh] operations: split/transform/overlap") {
    TEST_CASE("split_3d_box_in_octants yields 8 boxes covering halves") {
        bbox_mesh  b(0.0, 2.0, 0.0, 4.0, 0.0, 6.0);
        const auto center = b.get_center();  // (1,2,3)
        auto       oct    = b.split_3d_box_in_octants();
        REQUIRE_EQ(oct.size(), 8);

        // Each octant has half-extent along each axis
        for (const auto& o : oct) {
            CHECK_EQ(o.get_x_min() == 0.0 || o.get_x_min() == center.x(), true);
            CHECK_EQ(o.get_x_max() == center.x() || o.get_x_max() == 2.0, true);
            CHECK_EQ(o.get_y_min() == 0.0 || o.get_y_min() == center.y(), true);
            CHECK_EQ(o.get_y_max() == center.y() || o.get_y_max() == 4.0, true);
            CHECK_EQ(o.get_z_min() == 0.0 || o.get_z_min() == center.z(), true);
            CHECK_EQ(o.get_z_max() == center.z() || o.get_z_max() == 6.0, true);

            // Size along x is exactly half
            CHECK_EQ(o.get_x_max() - o.get_x_min(), doctest::Approx(1.0));
            // y half
            CHECK_EQ(o.get_y_max() - o.get_y_min(), doctest::Approx(2.0));
            // z half
            CHECK_EQ(o.get_z_max() - o.get_z_min(), doctest::Approx(3.0));
        }
    }

    TEST_CASE("translate shifts min/max by the vector") {
        bbox_mesh b(0.0, 1.0, 10.0, 12.0, -3.0, -2.0);
        b.translate(vector3(2.0, -5.0, 4.0));
        CHECK_EQ(b.get_x_min(), doctest::Approx(2.0));
        CHECK_EQ(b.get_x_max(), doctest::Approx(3.0));
        CHECK_EQ(b.get_y_min(), doctest::Approx(5.0));
        CHECK_EQ(b.get_y_max(), doctest::Approx(7.0));
        CHECK_EQ(b.get_z_min(), doctest::Approx(1.0));
        CHECK_EQ(b.get_z_max(), doctest::Approx(2.0));
    }

    TEST_CASE("dilate scales from origin (not from center)") {
        bbox_mesh b(-1.0, 2.0, -2.0, 4.0, -3.0, 1.0);
        b.dilate(2.5);
        CHECK_EQ(b.get_x_min(), doctest::Approx(-2.5));
        CHECK_EQ(b.get_x_max(), doctest::Approx(5.0));
        CHECK_EQ(b.get_y_min(), doctest::Approx(-5.0));
        CHECK_EQ(b.get_y_max(), doctest::Approx(10.0));
        CHECK_EQ(b.get_z_min(), doctest::Approx(-7.5));
        CHECK_EQ(b.get_z_max(), doctest::Approx(2.5));
    }

    TEST_CASE("is_overlapping: overlapping, touching, and disjoint cases") {
        bbox_mesh a(0.0, 2.0, 0.0, 2.0, 0.0, 2.0);

        // overlapping interior
        bbox_mesh b(1.0, 3.0, 1.0, 3.0, 1.0, 3.0);
        CHECK(a.is_overlapping(b));
        CHECK(b.is_overlapping(a));

        // touching at a face: equality on a boundary -> counts as overlapping with current '>' logic
        bbox_mesh c(2.0, 4.0, 0.0, 2.0, 0.0, 2.0);
        CHECK(a.is_overlapping(c));
        CHECK(c.is_overlapping(a));

        // touching at an edge
        bbox_mesh d(2.0, 4.0, 2.0, 5.0, 0.0, 2.0);
        CHECK(a.is_overlapping(d));

        // touching at a corner
        bbox_mesh e(2.0, 4.0, 2.0, 5.0, 2.0, 6.0);
        CHECK(a.is_overlapping(e));

        // fully disjoint
        bbox_mesh f(3.0, 5.0, 0.0, 2.0, 0.0, 2.0);
        CHECK_FALSE(a.is_overlapping(f));
    }

    TEST_CASE("ostream prints readable summary") {
        bbox_mesh          b(-1.0, 2.0, 0.5, 1.5, -3.0, 0.0);
        std::ostringstream oss;
        oss << b;
        auto s = oss.str();
        CHECK(s.find("bbox_mesh:") != std::string::npos);
        CHECK(s.find("x_min = -1") != std::string::npos);
        CHECK(s.find("x_max = 2") != std::string::npos);
        CHECK(s.find("y_min = 0.5") != std::string::npos);
        CHECK(s.find("y_max = 1.5") != std::string::npos);
        CHECK(s.find("z_min = -3") != std::string::npos);
        CHECK(s.find("z_max = 0") != std::string::npos);
    }
}
