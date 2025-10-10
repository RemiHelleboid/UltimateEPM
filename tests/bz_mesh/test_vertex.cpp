#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <array>

#include "mesh_vertex.hpp"

using bz_mesh::vector3;
using bz_mesh::Vertex;

TEST_SUITE("[Vertex] construction & position") {
    TEST_CASE("default ctor and getters") {
        Vertex v;
        CHECK_EQ(v.get_index(), std::size_t{0});
        CHECK_EQ(v.get_position().x(), 0.0);
        CHECK_EQ(v.get_position().y(), 0.0);
        CHECK_EQ(v.get_position().z(), 0.0);
        CHECK_EQ(v.get_number_bands(), std::size_t{0});
        CHECK(v.get_band_energies().empty());
        CHECK(v.get_electron_phonon_rates_all_bands().empty());
    }

    TEST_CASE("index-only ctor initializes index, position at origin") {
        Vertex v(7);
        CHECK_EQ(v.get_index(), std::size_t{7});
        CHECK_EQ(v.get_position().x(), 0.0);
        CHECK_EQ(v.get_position().y(), 0.0);
        CHECK_EQ(v.get_position().z(), 0.0);
    }

    TEST_CASE("index + vector3 position ctor") {
        vector3 p(1.0, -2.0, 3.5);
        Vertex  v(3, p);
        CHECK_EQ(v.get_index(), std::size_t{3});
        CHECK_EQ(v.get_position().x(), doctest::Approx(1.0));
        CHECK_EQ(v.get_position().y(), doctest::Approx(-2.0));
        CHECK_EQ(v.get_position().z(), doctest::Approx(3.5));
    }

    TEST_CASE("index + xyz position ctor") {
        Vertex v(9, -1.0, 2.0, -3.0);
        CHECK_EQ(v.get_index(), std::size_t{9});
        CHECK_EQ(v.get_position().x(), doctest::Approx(-1.0));
        CHECK_EQ(v.get_position().y(), doctest::Approx(2.0));
        CHECK_EQ(v.get_position().z(), doctest::Approx(-3.0));
    }

    TEST_CASE("set_position & shift_position") {
        Vertex v(1);
        v.set_position(vector3(1.0, 2.0, 3.0));
        CHECK_EQ(v.get_position().x(), doctest::Approx(1.0));
        CHECK_EQ(v.get_position().y(), doctest::Approx(2.0));
        CHECK_EQ(v.get_position().z(), doctest::Approx(3.0));

        v.shift_position(vector3(-2.0, 0.5, 1.0));
        CHECK_EQ(v.get_position().x(), doctest::Approx(-1.0));
        CHECK_EQ(v.get_position().y(), doctest::Approx(2.5));
        CHECK_EQ(v.get_position().z(), doctest::Approx(4.0));
    }
}

TEST_SUITE("[Vertex] band energies") {
    TEST_CASE("add_band_energy_value and get_number_bands / get_energy_at_band") {
        Vertex v;
        v.add_band_energy_value(1.0);
        v.add_band_energy_value(2.5);
        v.add_band_energy_value(-0.75);

        CHECK_EQ(v.get_number_bands(), std::size_t{3});
        CHECK_EQ(v.get_energy_at_band(0), doctest::Approx(1.0));
        CHECK_EQ(v.get_energy_at_band(1), doctest::Approx(2.5));
        CHECK_EQ(v.get_energy_at_band(2), doctest::Approx(-0.75));

        // get_band_energies view matches content
        const auto& all = v.get_band_energies();
        REQUIRE_EQ(all.size(), std::size_t{3});
        CHECK_EQ(all[0], doctest::Approx(1.0));
        CHECK_EQ(all[1], doctest::Approx(2.5));
        CHECK_EQ(all[2], doctest::Approx(-0.75));
    }

    TEST_CASE("set_band_energy modifies existing index") {
        Vertex v;
        v.add_band_energy_value(10.0);
        v.add_band_energy_value(20.0);
        v.add_band_energy_value(30.0);

        v.set_band_energy(1, 42.0);
        CHECK_EQ(v.get_energy_at_band(0), doctest::Approx(10.0));
        CHECK_EQ(v.get_energy_at_band(1), doctest::Approx(42.0));
        CHECK_EQ(v.get_energy_at_band(2), doctest::Approx(30.0));
    }

    TEST_CASE("remove_band_energy erases and shifts following elements") {
        Vertex v;
        v.add_band_energy_value(0.0);
        v.add_band_energy_value(1.0);
        v.add_band_energy_value(2.0);
        v.add_band_energy_value(3.0);

        v.remove_band_energy(1);  // remove the '1.0'
        REQUIRE_EQ(v.get_number_bands(), std::size_t{3});
        CHECK_EQ(v.get_energy_at_band(0), doctest::Approx(0.0));
        CHECK_EQ(v.get_energy_at_band(1), doctest::Approx(2.0));  // shifted
        CHECK_EQ(v.get_energy_at_band(2), doctest::Approx(3.0));
    }

    TEST_CASE("swap_bands exchanges energies at indices") {
        Vertex v;
        v.add_band_energy_value(5.0);  // 0
        v.add_band_energy_value(6.0);  // 1
        v.add_band_energy_value(7.0);  // 2

        v.swap_bands(0, 2);
        CHECK_EQ(v.get_energy_at_band(0), doctest::Approx(7.0));
        CHECK_EQ(v.get_energy_at_band(1), doctest::Approx(6.0));
        CHECK_EQ(v.get_energy_at_band(2), doctest::Approx(5.0));
    }

    TEST_CASE("remove_band_energy throws on OOB index") {
        Vertex v;
        v.add_band_energy_value(1.0);
        CHECK_THROWS_AS(v.remove_band_energy(3), std::invalid_argument);
        CHECK_THROWS_AS(v.remove_band_energy(1), std::invalid_argument);  // size==1, index==1 invalid
    }

    TEST_CASE("set_band_energy throws when index > size()") {
        Vertex v;
        v.add_band_energy_value(9.0);
        // The implementation checks '>' (not '>='); index==size() would be UB.
        // Use size()+1 to ensure it throws without invoking UB.
        const std::size_t bad = v.get_number_bands() + 1;
        CHECK_THROWS_AS(v.set_band_energy(bad, 1.23), std::invalid_argument);
    }

    TEST_CASE("swap_bands throws if either index is OOB") {
        Vertex v;
        v.add_band_energy_value(0.0);
        v.add_band_energy_value(1.0);
        CHECK_THROWS_AS(v.swap_bands(0, 2), std::invalid_argument);
        CHECK_THROWS_AS(v.swap_bands(5, 1), std::invalid_argument);
    }
}

TEST_SUITE("[Vertex] electron-phonon rates") {
    TEST_CASE("add and get per-band rates") {
        Vertex v;

        std::array<double, 8> r0{{0, 1, 2, 3, 4, 5, 6, 7}};
        std::array<double, 8> r1{{7, 6, 5, 4, 3, 2, 1, 0}};

        v.add_electron_phonon_rates(r0);
        v.add_electron_phonon_rates(r1);

        const auto& all = v.get_electron_phonon_rates_all_bands();
        REQUIRE_EQ(all.size(), std::size_t{2});

        const auto& g0 = v.get_electron_phonon_rates(0);
        const auto& g1 = v.get_electron_phonon_rates(1);

        for (std::size_t i = 0; i < 8; ++i) {
            CHECK_EQ(g0[i], doctest::Approx(r0[i]));
            CHECK_EQ(g1[i], doctest::Approx(r1[i]));
        }
    }
}
