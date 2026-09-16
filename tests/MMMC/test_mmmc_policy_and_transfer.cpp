#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>

#include "mmmc_particle_pools.hpp"
#include "mmmc_particle_transfer.hpp"
#include "mmmc_policy.hpp"
#include "unit_conversion.hpp"

TEST_CASE("MMMC bbox policy selects PBMC inside the configured box and ADMC outside") {
    const uepm::MMMC::bbox_transport_policy policy{uepm::mesh::bbox{0.0, 1.0, -1.0, 1.0, 0.0, 0.0}};

    CHECK(policy.method_for_position({0.5, 0.0, 0.0}) == uepm::MMMC::transport_method::pbmc);
    CHECK(policy.method_for_position({1.5, 0.0, 0.0}) == uepm::MMMC::transport_method::admc);
    CHECK(uepm::MMMC::transport_method_name(uepm::MMMC::transport_method::pbmc) == "PBMC");
}

TEST_CASE("MMMC bbox policy smoothly blends transport methods through its outer buffer") {
    const uepm::MMMC::bbox_transport_policy policy{
        .m_pbmc_region_um = uepm::mesh::bbox{0.0, 1.0, -1.0, 1.0, 0.0, 0.0},
        .m_buffer_width_um = 1.0,
    };

    CHECK(policy.admc_probability_at({0.5, 0.0, 0.0}) == doctest::Approx(0.0));
    CHECK(policy.admc_probability_at({1.25, 0.0, 0.0}) == doctest::Approx(0.15625));
    CHECK(policy.admc_probability_at({1.25, 1.75, 0.0}) == doctest::Approx(0.84375));
    CHECK(policy.admc_probability_at({2.0, 0.0, 0.0}) == doctest::Approx(1.0));
    CHECK(policy.admc_probability_at({3.0, 0.0, 0.0}) == doctest::Approx(1.0));
}

TEST_CASE("MMMC bbox buffer gives particles reproducible switching thresholds") {
    const uepm::MMMC::bbox_transport_policy policy{
        .m_pbmc_region_um = uepm::mesh::bbox{0.0, 1.0, -1.0, 1.0, 0.0, 0.0},
        .m_buffer_width_um = 1.0,
    };

    CHECK(policy.particle_threshold(42, 1234) == policy.particle_threshold(42, 1234));
    CHECK(policy.particle_threshold(42, 1234) != policy.particle_threshold(43, 1234));
    CHECK(policy.particle_threshold(42, 1234) != policy.particle_threshold(42, 1235));

    const auto first = policy.method_for_particle({1.5, 0.0, 0.0}, 42, 1234);
    for (std::size_t repetition = 0; repetition < 10; ++repetition) {
        CHECK(policy.method_for_particle({1.5, 0.0, 0.0}, 42, 1234) == first);
    }
}

TEST_CASE("MMMC bbox buffer produces its configured smoothstep population blend") {
    const uepm::MMMC::bbox_transport_policy policy{
        .m_pbmc_region_um = uepm::mesh::bbox{0.0, 1.0, -1.0, 1.0, 0.0, 0.0},
        .m_buffer_width_um = 1.0,
    };
    constexpr std::size_t sample_count = 20000;
    std::size_t           admc_count   = 0;
    for (std::size_t particle_index = 0; particle_index < sample_count; ++particle_index) {
        if (policy.method_for_particle({1.3, 0.0, 0.0}, particle_index, 9876) ==
            uepm::MMMC::transport_method::admc) {
            ++admc_count;
        }
    }

    const double admc_fraction = static_cast<double>(admc_count) / static_cast<double>(sample_count);
    constexpr double expected_probability = 0.3 * 0.3 * (3.0 - 2.0 * 0.3);
    CHECK(admc_fraction == doctest::Approx(expected_probability).epsilon(0.04));
}

TEST_CASE("MMMC bbox policy rejects invalid buffer widths") {
    auto policy = uepm::MMMC::bbox_transport_policy{
        .m_pbmc_region_um = uepm::mesh::bbox{0.0, 1.0, -1.0, 1.0, 0.0, 0.0},
        .m_buffer_width_um = -0.1,
    };
    CHECK_THROWS_AS(policy.validate(), std::invalid_argument);

    policy.m_buffer_width_um = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(policy.validate(), std::invalid_argument);
}

TEST_CASE("MMMC PBMC to ADMC transfer preserves canonical state and converts position units") {
    uepm::PBMC::particle_state state{};
    state.time                      = 2.0e-13;
    state.position                  = {2.0, 3.0, 0.0};
    state.previous_position         = {1.0, 2.0, 0.0};
    state.velocity                  = {10.0, -20.0, 0.0};
    state.electric_field            = {1.0e5, 2.0e5, 0.0};
    state.doping_concentration_cm_3 = 4.0e15;
    state.lattice_temperature_K     = 310.0;

    const uepm::PBMC::pbmc_particle pbmc_particle(7, uepm::PBMC::particle_type::electron, state, 3.0);
    const auto                      admc_particle = uepm::MMMC::convert_pbmc_to_admc(pbmc_particle, 11);

    CHECK(admc_particle.particle.index() == 11);
    CHECK(admc_particle.particle.type() == uepm::ADMC::carrier_type::electron);
    CHECK(admc_particle.weight == doctest::Approx(3.0));
    CHECK(admc_particle.particle.state().time_s == doctest::Approx(state.time));
    CHECK(admc_particle.particle.state().position_m.x() == doctest::Approx(2.0 * uepm::units::micron_to_meter));
    CHECK(admc_particle.particle.state().previous_position_m.y() ==
          doctest::Approx(2.0 * uepm::units::micron_to_meter));
    CHECK(admc_particle.particle.state().total_velocity_m_per_s.x() == doctest::Approx(10.0));
    CHECK(admc_particle.particle.state().electric_field_V_per_m.y() ==
          doctest::Approx(2.0e5 * uepm::units::electric_field_V_per_cm_to_V_per_m));
}

TEST_CASE("MMMC carrier type conversion maps both carrier kinds") {
    CHECK(uepm::MMMC::to_admc_carrier_type(uepm::PBMC::particle_type::electron) == uepm::ADMC::carrier_type::electron);
    CHECK(uepm::MMMC::to_admc_carrier_type(uepm::PBMC::particle_type::hole) == uepm::ADMC::carrier_type::hole);
    CHECK(uepm::MMMC::to_pbmc_particle_type(uepm::ADMC::carrier_type::electron) == uepm::PBMC::particle_type::electron);
    CHECK(uepm::MMMC::to_pbmc_particle_type(uepm::ADMC::carrier_type::hole) == uepm::PBMC::particle_type::hole);
}

TEST_CASE("MMMC ADMC to PBMC transfer preserves canonical state and converts field units") {
    uepm::ADMC::device_admc_particle admc_particle{
        .particle = uepm::ADMC::admc_particle(4, uepm::ADMC::carrier_type::hole, {2.0e-6, 3.0e-6, 0.0}),
        .weight   = 2.5,
    };
    auto& state                     = admc_particle.particle.state();
    state.time_s                    = 4.0e-13;
    state.previous_position_m       = {1.0e-6, 2.0e-6, 0.0};
    state.electric_field_V_per_m    = {1.0e7, -2.0e7, 0.0};
    state.drift_velocity_m_per_s    = {10.0, 20.0, 0.0};
    state.total_velocity_m_per_s    = {30.0, 40.0, 0.0};
    state.doping_concentration_cm_3 = 5.0e16;
    state.lattice_temperature_K     = 325.0;

    uepm::PBMC::pbmc_transport_kernel transport;
    transport.initialize();
    const auto pbmc_particle = uepm::MMMC::convert_admc_to_pbmc(admc_particle, 9, transport);

    CHECK(pbmc_particle.index() == 9);
    CHECK(pbmc_particle.type() == uepm::PBMC::particle_type::hole);
    CHECK(pbmc_particle.weight() == doctest::Approx(2.5));
    CHECK(pbmc_particle.state().time == doctest::Approx(state.time_s));
    CHECK(pbmc_particle.state().position.x() == doctest::Approx(2.0));
    CHECK(pbmc_particle.state().previous_position.y() == doctest::Approx(2.0));
    CHECK(pbmc_particle.state().electric_field.x() == doctest::Approx(1.0e5));
    CHECK(pbmc_particle.state().electric_field.y() == doctest::Approx(-2.0e5));
}

TEST_CASE("MMMC ADMC to PBMC transfer ignores timestep-dependent displacement velocity when drift is zero") {
    uepm::ADMC::device_admc_particle first{
        .particle = uepm::ADMC::admc_particle(0, uepm::ADMC::carrier_type::electron, {0.0, 0.0, 0.0}),
    };
    uepm::ADMC::device_admc_particle second = first;
    first.particle.state().drift_velocity_m_per_s  = {0.0, 0.0, 0.0};
    second.particle.state().drift_velocity_m_per_s = {0.0, 0.0, 0.0};
    first.particle.state().total_velocity_m_per_s  = {1.0e8, 0.0, 0.0};
    second.particle.state().total_velocity_m_per_s = {-1.0e8, 2.0e8, 0.0};

    uepm::PBMC::pbmc_transport_config config;
    config.m_carrier_type = uepm::PBMC::particle_type::electron;
    uepm::PBMC::pbmc_transport_kernel first_transport(config, 1234);
    uepm::PBMC::pbmc_transport_kernel second_transport(config, 1234);
    first_transport.initialize();
    second_transport.initialize();

    const auto first_pbmc  = uepm::MMMC::convert_admc_to_pbmc(first, 0, first_transport);
    const auto second_pbmc = uepm::MMMC::convert_admc_to_pbmc(second, 0, second_transport);

    CHECK(first_pbmc.state().velocity.x() == doctest::Approx(second_pbmc.state().velocity.x()));
    CHECK(first_pbmc.state().velocity.y() == doctest::Approx(second_pbmc.state().velocity.y()));
    CHECK(first_pbmc.state().velocity.z() == doctest::Approx(second_pbmc.state().velocity.z()));
    CHECK(first_pbmc.state().kinetic_energy == doctest::Approx(second_pbmc.state().kinetic_energy));
}

TEST_CASE("MMMC ADMC to PBMC transfer conserves drift velocity statistically") {
    constexpr std::size_t sample_count = 20000;
    const uepm::mesh::vector3 requested_drift{2.0e4, -1.0e4, 0.0};

    uepm::ADMC::device_admc_particle admc_particle{
        .particle = uepm::ADMC::admc_particle(0, uepm::ADMC::carrier_type::electron, {0.0, 0.0, 0.0}),
    };
    auto& admc_state                     = admc_particle.particle.state();
    admc_state.drift_velocity_m_per_s    = requested_drift;
    admc_state.total_velocity_m_per_s    = {9.0e7, -8.0e7, 7.0e7};
    admc_state.lattice_temperature_K     = 300.0;

    uepm::PBMC::pbmc_transport_config config;
    config.m_carrier_type = uepm::PBMC::particle_type::electron;
    uepm::PBMC::pbmc_transport_kernel transport(config, 9876);
    transport.initialize();

    uepm::mesh::vector3 velocity_sum{};
    for (std::size_t i = 0; i < sample_count; ++i) {
        velocity_sum += uepm::MMMC::convert_admc_to_pbmc(admc_particle, i, transport).state().velocity;
    }
    const auto mean_velocity = velocity_sum / static_cast<double>(sample_count);

    CHECK(mean_velocity.x() == doctest::Approx(requested_drift.x()).epsilon(0.12));
    CHECK(mean_velocity.y() == doctest::Approx(requested_drift.y()).epsilon(0.12));
    CHECK(std::abs(mean_velocity.z()) < 2.5e3);
}

TEST_CASE("MMMC particle pools transfer PBMC particles outside the PBMC bbox into ADMC") {
    uepm::MMMC::particle_pools pools;
    uepm::PBMC::particle_state inside_state{};
    inside_state.position = {0.5, 0.0, 0.0};
    uepm::PBMC::particle_state outside_state{};
    outside_state.position = {2.0, 0.0, 0.0};

    pools.m_pbmc_particles.push_back(
        std::make_unique<uepm::PBMC::pbmc_particle>(0, uepm::PBMC::particle_type::electron, inside_state, 1.0));
    pools.m_pbmc_particles.push_back(
        std::make_unique<uepm::PBMC::pbmc_particle>(1, uepm::PBMC::particle_type::hole, outside_state, 2.0));

    uepm::PBMC::pbmc_transport_kernel       electron_transport;
    uepm::PBMC::pbmc_transport_kernel       hole_transport;
    const uepm::MMMC::bbox_transport_policy policy{uepm::mesh::bbox{0.0, 1.0, -1.0, 1.0, 0.0, 0.0}};

    pools.apply_policy(policy, 0, electron_transport, hole_transport);

    CHECK(pools.pbmc_size() == 1);
    CHECK(pools.admc_size() == 1);
    CHECK(pools.m_last_transfer_counters.m_pbmc_to_admc == 1);
    CHECK(pools.m_last_transfer_counters.m_admc_to_pbmc == 0);
    CHECK(pools.m_admc_particles.front().particle.type() == uepm::ADMC::carrier_type::hole);
    CHECK(pools.m_admc_particles.front().weight == doctest::Approx(2.0));
}

TEST_CASE("MMMC particle pools transfer both ways at particle-specific positions in the bbox buffer") {
    constexpr std::uint64_t policy_seed = 2468;
    const uepm::MMMC::bbox_transport_policy policy{
        .m_pbmc_region_um = uepm::mesh::bbox{0.0, 1.0, -1.0, 1.0, 0.0, 0.0},
        .m_buffer_width_um = 1.0,
    };

    constexpr std::size_t pbmc_index = 100;
    constexpr std::size_t admc_index = 101;
    const double pbmc_buffer_coordinate = 0.5 * (1.0 + policy.particle_threshold(pbmc_index, policy_seed));
    const double admc_buffer_coordinate = 0.5 * policy.particle_threshold(admc_index, policy_seed);

    uepm::PBMC::particle_state pbmc_state{};
    pbmc_state.position = {1.0 + pbmc_buffer_coordinate, 0.0, 0.0};

    uepm::MMMC::particle_pools pools;
    pools.m_pbmc_particles.push_back(
        std::make_unique<uepm::PBMC::pbmc_particle>(pbmc_index,
                                                    uepm::PBMC::particle_type::electron,
                                                    pbmc_state,
                                                    1.0));
    pools.m_admc_particles.push_back({
        .particle = uepm::ADMC::admc_particle(
            admc_index,
            uepm::ADMC::carrier_type::hole,
            {(1.0 + admc_buffer_coordinate) * uepm::units::micron_to_meter, 0.0, 0.0}),
        .weight = 2.0,
    });

    uepm::PBMC::pbmc_transport_config electron_config;
    electron_config.m_carrier_type = uepm::PBMC::particle_type::electron;
    uepm::PBMC::pbmc_transport_config hole_config;
    hole_config.m_carrier_type = uepm::PBMC::particle_type::hole;
    uepm::PBMC::pbmc_transport_kernel electron_transport(electron_config, 1);
    uepm::PBMC::pbmc_transport_kernel hole_transport(hole_config, 2);
    electron_transport.initialize();
    hole_transport.initialize();

    pools.apply_policy(policy, policy_seed, electron_transport, hole_transport);

    CHECK(pools.pbmc_size() == 1);
    CHECK(pools.admc_size() == 1);
    CHECK(pools.m_pbmc_particles.front()->index() == admc_index);
    CHECK(pools.m_admc_particles.front().particle.index() == pbmc_index);
    CHECK(pools.m_last_transfer_counters.m_pbmc_to_admc == 1);
    CHECK(pools.m_last_transfer_counters.m_admc_to_pbmc == 1);
}
