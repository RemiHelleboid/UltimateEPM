#include <optional>
#include <sstream>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"
#include "vector.hpp"

using uepm::mesh_bz::permutaion_type;
using uepm::mesh_bz::vector3;

TEST_SUITE("[vector3] basics") {
    TEST_CASE("default and value construction + accessors") {
        vector3 a;
        CHECK_EQ(a.x(), 0.0);
        CHECK_EQ(a.y(), 0.0);
        CHECK_EQ(a.z(), 0.0);

        vector3 b(1.5, -2.0);
        CHECK_EQ(b.x(), 1.5);
        CHECK_EQ(b.y(), -2.0);
        CHECK_EQ(b.z(), 0.0);

        vector3 c(1.0, 2.0, 3.0);
        CHECK_EQ(c.x(), 1.0);
        CHECK_EQ(c.y(), 2.0);
        CHECK_EQ(c.z(), 3.0);
    }

    TEST_CASE("mutators & to_2d") {
        vector3 v;
        v.set_x(4.0);
        v.set_y(-5.0);
        v.set_z(6.0);
        CHECK_EQ(v.x(), 4.0);
        CHECK_EQ(v.y(), -5.0);
        CHECK_EQ(v.z(), 6.0);

        v.set_coordinates(1.0, 2.0, 3.0);
        CHECK_EQ(v.x(), 1.0);
        CHECK_EQ(v.y(), 2.0);
        CHECK_EQ(v.z(), 3.0);

        auto v2 = v.to_2d();
        CHECK_EQ(v2.x(), 1.0);
        CHECK_EQ(v2.y(), 2.0);
        CHECK_EQ(v2.z(), 0.0);
    }

    TEST_CASE("stream output") {
        vector3            v(1.0, 2.0, 3.0);
        std::ostringstream oss;
        oss << v;
        CHECK_EQ(oss.str(), "1,2,3");
    }
}

TEST_SUITE("[vector3] arithmetic") {
    TEST_CASE("+=, -=, *=, /=") {
        vector3 v(1.0, 2.0, 3.0);
        vector3 w(-2.0, 4.0, 0.5);

        v += w;  // (-1, 6, 3.5)
        CHECK_EQ(v.x(), doctest::Approx(-1.0));
        CHECK_EQ(v.y(), doctest::Approx(6.0));
        CHECK_EQ(v.z(), doctest::Approx(3.5));

        v -= w;  // back to (1,2,3)
        CHECK_EQ(v.x(), doctest::Approx(1.0));
        CHECK_EQ(v.y(), doctest::Approx(2.0));
        CHECK_EQ(v.z(), doctest::Approx(3.0));

        v *= 2.0;  // (2,4,6)
        CHECK_EQ(v.x(), doctest::Approx(2.0));
        CHECK_EQ(v.y(), doctest::Approx(4.0));
        CHECK_EQ(v.z(), doctest::Approx(6.0));

        v /= 2.0;  // (1,2,3)
        CHECK_EQ(v.x(), doctest::Approx(1.0));
        CHECK_EQ(v.y(), doctest::Approx(2.0));
        CHECK_EQ(v.z(), doctest::Approx(3.0));
    }

    TEST_CASE("+, -, scalar *, /") {
        vector3 a(1.0, 2.0, 3.0);
        vector3 b(4.0, -1.0, 0.5);

        vector3 s = a + b;
        CHECK_EQ(s.x(), doctest::Approx(5.0));
        CHECK_EQ(s.y(), doctest::Approx(1.0));
        CHECK_EQ(s.z(), doctest::Approx(3.5));

        vector3 d = a - b;
        CHECK_EQ(d.x(), doctest::Approx(-3.0));
        CHECK_EQ(d.y(), doctest::Approx(3.0));
        CHECK_EQ(d.z(), doctest::Approx(2.5));

        vector3 m1 = a * 2.0;
        vector3 m2 = 2.0 * a;
        CHECK_EQ(m1.x(), doctest::Approx(2.0));
        CHECK_EQ(m1.y(), doctest::Approx(4.0));
        CHECK_EQ(m1.z(), doctest::Approx(6.0));
        CHECK_EQ(m2.x(), doctest::Approx(2.0));
        CHECK_EQ(m2.y(), doctest::Approx(4.0));
        CHECK_EQ(m2.z(), doctest::Approx(6.0));

        vector3 q = a / 2.0;
        CHECK_EQ(q.x(), doctest::Approx(0.5));
        CHECK_EQ(q.y(), doctest::Approx(1.0));
        CHECK_EQ(q.z(), doctest::Approx(1.5));
    }

    TEST_CASE("norm and norm_squared") {
        vector3 v(3.0, 4.0, 12.0);
        CHECK_EQ(v.norm_squared(), doctest::Approx(9.0 + 16.0 + 144.0));
        CHECK_EQ(v.norm(), doctest::Approx(std::sqrt(169.0)));
    }

    TEST_CASE("re_normalize: unit length for non-zero, no-op for zero") {
        vector3 v(3.0, 4.0, 0.0);
        v.re_normalize();
        CHECK_EQ(v.norm(), doctest::Approx(1.0).epsilon(1e-12));
        CHECK_EQ(v.x() / v.y(), doctest::Approx(3.0 / 4.0));

        vector3 z(0.0, 0.0, 0.0);
        z.re_normalize();  // should not crash or change values
        CHECK_EQ(z.x(), 0.0);
        CHECK_EQ(z.y(), 0.0);
        CHECK_EQ(z.z(), 0.0);
    }
}

TEST_SUITE("[vector3] geometry helpers") {
    TEST_CASE("dot & cross & scalar triple") {
        vector3 i(1.0, 0.0, 0.0);
        vector3 j(0.0, 1.0, 0.0);
        vector3 k(0.0, 0.0, 1.0);

        CHECK_EQ(i.dot(j), doctest::Approx(0.0));
        CHECK_EQ(i.dot(i), doctest::Approx(1.0));

        vector3 c = cross_product(i, j);
        CHECK_EQ(c.x(), doctest::Approx(0.0));
        CHECK_EQ(c.y(), doctest::Approx(0.0));
        CHECK_EQ(c.z(), doctest::Approx(1.0));

        CHECK_EQ(scalar_triple_product(i, j, k), doctest::Approx(1.0));
    }

    TEST_CASE("middle, point_pair_to_vector, distance") {
        vector3 a(1.0, 2.0, 3.0);
        vector3 b(3.0, 6.0, 7.0);

        vector3 mid = middle(a, b);
        CHECK_EQ(mid.x(), doctest::Approx(2.0));
        CHECK_EQ(mid.y(), doctest::Approx(4.0));
        CHECK_EQ(mid.z(), doctest::Approx(5.0));

        vector3 ab = point_pair_to_vector(a, b);
        CHECK_EQ(ab.x(), doctest::Approx(2.0));
        CHECK_EQ(ab.y(), doctest::Approx(4.0));
        CHECK_EQ(ab.z(), doctest::Approx(4.0));

        CHECK_EQ(distance(a, b), doctest::Approx(std::sqrt(2.0 * 2.0 + 4.0 * 4.0 + 4.0 * 4.0)));
    }

    TEST_CASE("compte_cos_angle normal case and zero-norm behavior") {
        vector3 a(1.0, 0.0, 0.0);
        vector3 b(1.0, 1.0, 0.0);
        CHECK_EQ(compute_cos_angle(a, b), doctest::Approx(1.0 / std::sqrt(2.0)));

        vector3 z(0.0, 0.0, 0.0);
        CHECK_EQ(compute_cos_angle(a, z), doctest::Approx(1.0));
        CHECK_EQ(compute_cos_angle(z, z), doctest::Approx(1.0));
    }

    TEST_CASE("is_point_between_two_others respects epsilon") {
        vector3 A(0.0, 0.0, 0.0);
        vector3 B(2.0, 0.0, 0.0);
        vector3 P_on(1.0, 0.0, 0.0);
        vector3 P_off(3.0, 0.0, 0.0);
        vector3 P_near(1.0, 1e-10, 0.0);

        CHECK(is_point_between_two_others(A, B, P_on));
        CHECK_FALSE(is_point_between_two_others(A, B, P_off));
        CHECK(is_point_between_two_others(A, B, P_near));  // within epsilon
    }

    TEST_CASE("double_scalar_product_2d is 2D determinant") {
        vector3 u(2.0, 1.0, 7.0);
        vector3 v(3.0, 5.0, 9.0);
        CHECK_EQ(double_scalar_product_2d(u, v), doctest::Approx(2.0 * 5.0 - 1.0 * 3.0));  // 10 - 3 = 7
    }
}

TEST_SUITE("[vector3] transformations") {
    TEST_CASE("apply_reflection") {
        vector3 v(1.0, -2.0, 3.0);
        v.apply_reflection(-1, 1, -1);
        CHECK_EQ(v.x(), doctest::Approx(-1.0));
        CHECK_EQ(v.y(), doctest::Approx(-2.0));
        CHECK_EQ(v.z(), doctest::Approx(-3.0));
    }

    TEST_CASE("apply_permutation XY/XZ/YZ/XYZ") {
        vector3 base(1.0, 2.0, 3.0);

        {
            vector3 v = base;
            v.apply_permutation(permutaion_type::XY);
            CHECK_EQ(v.x(), doctest::Approx(2.0));
            CHECK_EQ(v.y(), doctest::Approx(1.0));
            CHECK_EQ(v.z(), doctest::Approx(3.0));
        }
        {
            vector3 v = base;
            v.apply_permutation(permutaion_type::XZ);
            CHECK_EQ(v.x(), doctest::Approx(3.0));
            CHECK_EQ(v.y(), doctest::Approx(2.0));
            CHECK_EQ(v.z(), doctest::Approx(1.0));
        }
        {
            vector3 v = base;
            v.apply_permutation(permutaion_type::YZ);
            CHECK_EQ(v.x(), doctest::Approx(1.0));
            CHECK_EQ(v.y(), doctest::Approx(3.0));
            CHECK_EQ(v.z(), doctest::Approx(2.0));
        }
        {
            vector3 v = base;
            v.apply_permutation(permutaion_type::XYZ);
            CHECK_EQ(v.x(), doctest::Approx(1.0));
            CHECK_EQ(v.y(), doctest::Approx(2.0));
            CHECK_EQ(v.z(), doctest::Approx(3.0));
        }
    }

    TEST_CASE("apply_permutation YZX and ZXY are distinct 3-cycles") {
        vector3 base(1.0, 2.0, 3.0);

        vector3 yzx = base;
        yzx.apply_permutation(permutaion_type::YZX);  // (y,z,x)
        CHECK_EQ(yzx.x(), doctest::Approx(2.0));
        CHECK_EQ(yzx.y(), doctest::Approx(3.0));
        CHECK_EQ(yzx.z(), doctest::Approx(1.0));

        vector3 zxy = base;
        zxy.apply_permutation(permutaion_type::ZXY);  // (z,x,y)
        CHECK_EQ(zxy.x(), doctest::Approx(3.0));
        CHECK_EQ(zxy.y(), doctest::Approx(1.0));
        CHECK_EQ(zxy.z(), doctest::Approx(2.0));
    }
}

TEST_SUITE("[vector3] segment intersection") {
    TEST_CASE("parallel non-intersecting segments return nullopt") {
        vector3 A(0.0, 0.0, 0.0);
        vector3 B(2.0, 0.0, 0.0);
        vector3 C(0.0, 1.0, 0.0);
        vector3 D(2.0, 1.0, 0.0);

        auto P = compute_line_line_intersection(A, B, C, D);
        CHECK_FALSE(P.has_value());
    }

    TEST_CASE("orthogonal intersecting segments return intersection point") {
        vector3 A(0.0, 0.0, 0.0);
        vector3 B(2.0, 0.0, 0.0);
        vector3 C(1.0, -1.0, 0.0);
        vector3 D(1.0, 1.0, 0.0);

        auto P = compute_line_line_intersection(A, B, C, D);
        REQUIRE(P.has_value());
        CHECK_EQ(P->x(), doctest::Approx(1.0));
        CHECK_EQ(P->y(), doctest::Approx(0.0));
        CHECK_EQ(P->z(), doctest::Approx(0.0));
    }

    TEST_CASE("skew segments with clear no-intersection") {
        // Segment 1 along x-axis from 0 to 1
        vector3 A(0.0, 0.0, 0.0);
        vector3 B(1.0, 0.0, 0.0);

        // Segment 2 is above and to the right; no overlap with [A,B]
        vector3 C(2.0, 1.0, 0.0);
        vector3 D(3.0, 2.0, 0.0);

        auto P = compute_line_line_intersection(A, B, C, D);
        CHECK_FALSE(P.has_value());
    }
}
