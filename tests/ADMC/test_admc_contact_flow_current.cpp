#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <string>

#include "device_admc_simulation.hpp"
#include "element2d.hpp"
#include "physical_constants.hpp"
#include "vertex.hpp"

namespace {

class contact_flow_harness : public uepm::ADMC::device_admc_simulation {
 public:
    using device_admc_simulation::device_admc_simulation;

    void begin_step(double dt_s) { reset_step_contact_flow(dt_s); }

    void collect(const uepm::ADMC::device_admc_particle& particle, std::size_t contact_index) {
        record_contact_collection(particle, contact_index);
    }
    void inject(uepm::ADMC::carrier_type type, double weight, std::size_t contact_index) {
        record_contact_injection(type, weight, contact_index);
    }

    void process_crossing(const uepm::mesh::vector3& previous_um,
                          const uepm::mesh::vector3& current_um,
                          uepm::ADMC::carrier_type  type,
                          double                    weight,
                          uepm::mesh::element*      containing_element) {
        uepm::ADMC::device_admc_particle particle{
            .particle = uepm::ADMC::admc_particle(0, type, previous_um * 1.0e-6),
            .containing_element = containing_element,
            .weight = weight,
            .crossed_contact = false,
        };
        particle.particle.state().previous_position_m = previous_um * 1.0e-6;
        particle.particle.state().position_m = current_um * 1.0e-6;
        update_element_and_check_boundary(particle);
        m_particles.push_back(std::move(particle));
        remove_collected_particles();
    }
};

TEST_CASE("ADMC collected-contact current uses signed particle weight and the actual step duration") {
    uepm::mesh::mesh mesh;
    mesh.set_dimension(2);
    uepm::device::device device{&mesh};
    device.add_contact("source", {-1.0, -1.0, -1.0}, {-0.9, 1.0, 1.0});
    device.add_contact("drain", {0.9, -1.0, -1.0}, {1.0, 1.0, 1.0});

    uepm::ADMC::options_device_ADMC options;
    contact_flow_harness simulation{device, options, 7};
    constexpr double shortened_step_s = 2.5e-16;
    simulation.begin_step(shortened_step_s);

    const uepm::ADMC::device_admc_particle electron{
        .particle = uepm::ADMC::admc_particle(0, uepm::ADMC::carrier_type::electron),
        .weight = 4.0,
    };
    const uepm::ADMC::device_admc_particle hole{
        .particle = uepm::ADMC::admc_particle(1, uepm::ADMC::carrier_type::hole),
        .weight = 1.5,
    };
    simulation.collect(electron, 1);
    simulation.collect(hole, 1);

    const auto electron_currents = simulation.collected_contact_electron_currents_A();
    const auto hole_currents = simulation.collected_contact_hole_currents_A();
    const auto total_currents = simulation.collected_contact_total_currents_A();
    const auto& cumulative_charges = simulation.cumulative_collected_contact_charges_C();

    REQUIRE(total_currents.size() == 2);
    CHECK(electron_currents[1] == doctest::Approx(-4.0 * uepm::constants::q_e / shortened_step_s));
    CHECK(hole_currents[1] == doctest::Approx(1.5 * uepm::constants::q_e / shortened_step_s));
    CHECK(total_currents[1] == doctest::Approx(-2.5 * uepm::constants::q_e / shortened_step_s));
    CHECK(cumulative_charges[1] == doctest::Approx(-2.5 * uepm::constants::q_e));
    simulation.inject(uepm::ADMC::carrier_type::electron, 1.0, 1);
    CHECK(simulation.net_contact_currents_A()[1] ==
          doctest::Approx(-1.5 * uepm::constants::q_e / shortened_step_s));

    simulation.begin_step(1.0e-16);
    CHECK(simulation.collected_contact_total_currents_A()[1] == 0.0);
    CHECK(simulation.cumulative_collected_contact_charges_C()[1] == doctest::Approx(-2.5 * uepm::constants::q_e));
}

TEST_CASE("ADMC boundary processing records a named segment crossing before particle removal") {
    uepm::mesh::mesh mesh;
    mesh.set_dimension(2);
    uepm::device::device device{&mesh};
    device.add_contact("drain", {0.9, -1.0, -1.0}, {1.0, 1.0, 1.0});

    uepm::ADMC::options_device_ADMC options;
    contact_flow_harness simulation{device, options, 11};
    simulation.begin_step(5.0e-16);

    uepm::mesh::vertex v0{0, 0.0, 0.0, 0.0};
    uepm::mesh::vertex v1{1, 0.5, 0.0, 0.0};
    uepm::mesh::vertex v2{2, 0.0, 0.5, 0.0};
    uepm::mesh::element2d element{&v0, &v1, &v2};
    simulation.process_crossing({0.1, 0.1, 0.0},
                                {1.2, 0.1, 0.0},
                                uepm::ADMC::carrier_type::electron,
                                2.0,
                                &element);

    CHECK(simulation.get_number_electrons() == 0);
    CHECK(simulation.collected_contact_total_currents_A()[0] ==
          doctest::Approx(-2.0 * uepm::constants::q_e / 5.0e-16));
}

TEST_CASE("ADMC history exports per-contact collected current and cumulative charge") {
    uepm::ADMC::history_device_ADMC history;
    history.set_contact_flow_names({"drain"});
    history.add(1.0, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, {}, {-2.0}, {0.5}, {-1.5}, {-3.0e-19},
                {-0.5}, {-1.0}, {-1.0e-19}, {-2.0e-19});

    const auto filename = std::filesystem::temp_directory_path() / "ultimate_epm_admc_contact_flow_history.csv";
    history.export_to_csv(filename.string());

    std::ifstream stream(filename);
    REQUIRE(stream.is_open());
    const std::string contents{std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
    CHECK(contents.find("collected_current_electron_drain_A,collected_current_hole_drain_A,") !=
          std::string::npos);
    CHECK(contents.find("collected_current_drain_A,cumulative_collected_charge_drain_C") != std::string::npos);
    CHECK(contents.ends_with(",-2,0.5,-1.5,-2.9999999999999999e-19,-0.5,-1,-9.9999999999999998e-20,"
                             "-2e-19\n"));
}

}  // namespace
