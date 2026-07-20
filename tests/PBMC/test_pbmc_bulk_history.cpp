#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>

#include "bulk_pbmc_simulation.hpp"

TEST_CASE("PBMC self-scattering history records each event time once") {
    uepm::PBMC::bulk_pbmc_simulation_config config;
    config.m_number_of_particles      = 1;
    config.m_final_time               = 1.0e-14;
    config.m_record_history           = true;
    config.m_enable_impact_ionization = false;
    config.m_nb_threads               = 1;

    uepm::PBMC::bulk_pbmc_simulation simulation(config);
    simulation.initialize();
    simulation.run_self_scattering_emc();

    REQUIRE(simulation.particles().size() == 1);
    const auto& snapshots = simulation.particles().front().history().snapshots();
    REQUIRE(snapshots.size() >= 2);
    CHECK(snapshots.front().time == doctest::Approx(0.0));
    CHECK(snapshots.back().time == doctest::Approx(config.m_final_time));
    for (std::size_t index = 1; index < snapshots.size(); ++index) {
        CHECK(snapshots[index].time > snapshots[index - 1].time);
    }
}
