#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "pbmc_transport_kernel.hpp"

namespace {

uepm::PBMC::pbmc_transport_kernel make_electron_transport() {
    uepm::PBMC::pbmc_transport_config config;
    config.m_carrier_type             = uepm::PBMC::particle_type::electron;
    config.m_max_energy_eV            = 2.0;
    config.m_gamma_max_energy_samples = 100;

    uepm::PBMC::pbmc_transport_kernel transport(config, 1234);
    transport.initialize();
    return transport;
}

}  // namespace

TEST_CASE("scattering channels preserve the total rate") {
    auto transport = make_electron_transport();

    uepm::PBMC::pbmc_particle particle(0, uepm::PBMC::particle_type::electron);
    transport.initialize_particle_state(particle);

    const auto channels = transport.build_scattering_channels(particle);
    REQUIRE_FALSE(channels.empty());

    double channel_rate_sum = 0.0;
    for (const auto& channel : channels) {
        channel_rate_sum += channel.rate_s_1;
    }

    CHECK(channel_rate_sum == doctest::Approx(transport.total_scattering_rate(particle)));
}

TEST_CASE("per-valley gamma bound is bounded by the global maximum") {
    auto transport = make_electron_transport();

    uepm::PBMC::pbmc_particle particle(0, uepm::PBMC::particle_type::electron);
    transport.initialize_particle_state(particle);

    CHECK(transport.gamma_max(particle) > 0.0);
    CHECK(transport.gamma_max(particle) <= transport.gamma_max());
}

TEST_CASE("zero-duration scattering does not change the particle") {
    auto transport = make_electron_transport();

    uepm::PBMC::pbmc_particle particle(0, uepm::PBMC::particle_type::electron);
    transport.initialize_particle_state(particle);
    const auto initial_state = particle.state();

    CHECK_FALSE(transport.scatter_particle(particle, 0.0).has_value());
    CHECK(particle.state().kinetic_energy == doctest::Approx(initial_state.kinetic_energy));
    CHECK(particle.state().valley_index == initial_state.valley_index);
}
