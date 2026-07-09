#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>

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

TEST_CASE("MMMC PBMC to ADMC transfer preserves canonical state and converts position units") {
    uepm::PBMC::particle_state state{};
    state.time                         = 2.0e-13;
    state.position                     = {2.0, 3.0, 0.0};
    state.previous_position            = {1.0, 2.0, 0.0};
    state.velocity                     = {10.0, -20.0, 0.0};
    state.electric_field               = {1.0e5, 2.0e5, 0.0};
    state.doping_concentration_cm_3    = 4.0e15;
    state.lattice_temperature_K        = 310.0;

    const uepm::PBMC::pbmc_particle pbmc_particle(7, uepm::PBMC::particle_type::electron, state, 3.0);
    const auto admc_particle = uepm::MMMC::convert_pbmc_to_admc(pbmc_particle, 11);

    CHECK(admc_particle.particle.index() == 11);
    CHECK(admc_particle.particle.type() == uepm::ADMC::carrier_type::electron);
    CHECK(admc_particle.weight == doctest::Approx(3.0));
    CHECK(admc_particle.particle.state().time_s == doctest::Approx(state.time));
    CHECK(admc_particle.particle.state().position_m.x() ==
          doctest::Approx(2.0 * uepm::units::micron_to_meter));
    CHECK(admc_particle.particle.state().previous_position_m.y() ==
          doctest::Approx(2.0 * uepm::units::micron_to_meter));
    CHECK(admc_particle.particle.state().total_velocity_m_per_s.x() == doctest::Approx(10.0));
    CHECK(admc_particle.particle.state().electric_field_V_per_m.y() == doctest::Approx(2.0e5));
}

TEST_CASE("MMMC carrier type conversion maps both carrier kinds") {
    CHECK(uepm::MMMC::to_admc_carrier_type(uepm::PBMC::particle_type::electron) ==
          uepm::ADMC::carrier_type::electron);
    CHECK(uepm::MMMC::to_admc_carrier_type(uepm::PBMC::particle_type::hole) ==
          uepm::ADMC::carrier_type::hole);
    CHECK(uepm::MMMC::to_pbmc_particle_type(uepm::ADMC::carrier_type::electron) ==
          uepm::PBMC::particle_type::electron);
    CHECK(uepm::MMMC::to_pbmc_particle_type(uepm::ADMC::carrier_type::hole) ==
          uepm::PBMC::particle_type::hole);
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

    uepm::PBMC::pbmc_transport_kernel electron_transport;
    uepm::PBMC::pbmc_transport_kernel hole_transport;
    const uepm::MMMC::bbox_transport_policy policy{uepm::mesh::bbox{0.0, 1.0, -1.0, 1.0, 0.0, 0.0}};

    pools.apply_policy(policy, electron_transport, hole_transport);

    CHECK(pools.pbmc_size() == 1);
    CHECK(pools.admc_size() == 1);
    CHECK(pools.m_last_transfer_counters.m_pbmc_to_admc == 1);
    CHECK(pools.m_last_transfer_counters.m_admc_to_pbmc == 0);
    CHECK(pools.m_admc_particles.front().particle.type() == uepm::ADMC::carrier_type::hole);
    CHECK(pools.m_admc_particles.front().weight == doctest::Approx(2.0));
}
