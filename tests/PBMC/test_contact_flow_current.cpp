#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include "device.hpp"
#include "device_pbmc_simulation.hpp"
#include "element2d.hpp"
#include "pbmc_device_history.hpp"
#include "vertex.hpp"

namespace {

class contact_flow_harness : public uepm::PBMC::device_pbmc_simulation {
 public:
    using device_pbmc_simulation::device_pbmc_simulation;

    void use_constant_weighting_field_per_cm(const uepm::mesh::vector3& weighting_field_per_cm) {
        m_state.m_use_constant_RamoUnitaryElectricField = true;
        m_state.m_RamoUnitaryElectricField_Vm_per_cm    = weighting_field_per_cm;
    }

    void collect(const uepm::PBMC::pbmc_particle& particle, std::size_t contact_index) {
        record_contact_collection(particle, contact_index);
    }

    void reset_step() { reset_step_contact_flow(); }
    void inject(uepm::PBMC::particle_type type, double weight, std::size_t contact_index) {
        record_contact_injection(type, weight, contact_index);
    }

    void process_crossing(const uepm::mesh::vector3& previous,
                          const uepm::mesh::vector3& current,
                          uepm::PBMC::particle_type  type,
                          double                     weight,
                          uepm::mesh::element*       containing_element) {
        auto particle                       = std::make_unique<uepm::PBMC::pbmc_particle>(0, type, weight);
        particle->state().previous_position = previous;
        particle->set_position(current);
        particle->set_containing_element(containing_element);
        m_list_particles.push_back(std::move(particle));
        update_element_and_check_boundary();
        remove_collected_particles();
    }
};

TEST_CASE("PBMC Ramo current converts the weighting field from inverse centimeters to inverse meters") {
    uepm::mesh::mesh mesh;
    mesh.set_dimension(2);
    uepm::device::device device{&mesh};

    uepm::PBMC::options_device_PBMC options;
    options.m_time_step = 1.0e-16;
    contact_flow_harness simulation{device, options, 3};
    simulation.use_constant_weighting_field_per_cm({2.0, 0.0, 0.0});

    uepm::PBMC::pbmc_particle electron{0, uepm::PBMC::particle_type::electron, 3.0};
    electron.state().velocity = {4.0, 0.0, 0.0};

    constexpr double elementary_charge_C = 1.602176634e-19;
    const double expected_current_A = -3.0 * elementary_charge_C * 4.0 * 2.0 * 100.0;
    CHECK(simulation.compute_ramo_current_for_particle(electron) == doctest::Approx(expected_current_A));
}

TEST_CASE("device reports the first named contact crossed by a particle segment") {
    uepm::mesh::mesh     mesh;
    uepm::device::device device{&mesh};
    device.add_contact("far", {2.0, -1.0, -1.0}, {2.1, 1.0, 1.0});
    device.add_contact("near", {1.0, -1.0, -1.0}, {1.1, 1.0, 1.0});

    const auto crossing = device.find_first_contact_crossing({0.0, 0.0, 0.0}, {3.0, 0.0, 0.0});

    REQUIRE(crossing.has_value());
    CHECK(crossing->contact_name == "near");
    CHECK(crossing->contact_index == 1);
    CHECK(crossing->segment_fraction == doctest::Approx(1.0 / 3.0));
}

TEST_CASE("PBMC collected-contact current uses signed weighted particle charge per timestep") {
    uepm::mesh::mesh mesh;
    mesh.set_dimension(2);
    uepm::device::device device{&mesh};
    device.add_contact("source", {-1.0, -1.0, -1.0}, {-0.9, 1.0, 1.0});
    device.add_contact("drain", {0.9, -1.0, -1.0}, {1.0, 1.0, 1.0});

    uepm::PBMC::options_device_PBMC options;
    options.m_time_step = 1.0e-16;
    contact_flow_harness simulation{device, options, 3};

    uepm::PBMC::pbmc_particle electron{0, uepm::PBMC::particle_type::electron, 4.0};
    uepm::PBMC::pbmc_particle hole{1, uepm::PBMC::particle_type::hole, 1.5};
    simulation.collect(electron, 1);
    simulation.collect(hole, 1);

    constexpr double elementary_charge_C = 1.602176634e-19;
    const auto       electron_currents   = simulation.collected_contact_electron_currents_A();
    const auto       hole_currents       = simulation.collected_contact_hole_currents_A();
    const auto       total_currents      = simulation.collected_contact_total_currents_A();
    const auto&      cumulative_charges  = simulation.cumulative_collected_contact_charges_C();

    REQUIRE(electron_currents.size() == 2);
    CHECK(electron_currents[0] == 0.0);
    CHECK(electron_currents[1] == doctest::Approx(-4.0 * elementary_charge_C / options.m_time_step));
    CHECK(hole_currents[1] == doctest::Approx(1.5 * elementary_charge_C / options.m_time_step));
    CHECK(total_currents[1] == doctest::Approx(-2.5 * elementary_charge_C / options.m_time_step));
    CHECK(cumulative_charges[1] == doctest::Approx(-2.5 * elementary_charge_C));

    simulation.inject(uepm::PBMC::particle_type::electron, 1.0, 1);
    CHECK(simulation.injected_contact_currents_A()[1] == doctest::Approx(-elementary_charge_C / options.m_time_step));
    CHECK(simulation.net_contact_currents_A()[1] == doctest::Approx(-1.5 * elementary_charge_C / options.m_time_step));
    CHECK(simulation.cumulative_net_contact_charges_C()[1] == doctest::Approx(-1.5 * elementary_charge_C));

    simulation.reset_step();
    CHECK(simulation.collected_contact_total_currents_A()[1] == 0.0);
    CHECK(simulation.cumulative_collected_contact_charges_C()[1] == doctest::Approx(-2.5 * elementary_charge_C));
}

TEST_CASE("PBMC transport attributes a segment crossing to its contact before removing the particle") {
    uepm::mesh::mesh mesh;
    mesh.set_dimension(2);
    uepm::device::device device{&mesh};
    device.add_contact("drain", {0.9, -1.0, -1.0}, {1.0, 1.0, 1.0});

    uepm::PBMC::options_device_PBMC options;
    options.m_time_step = 5.0e-17;
    contact_flow_harness simulation{device, options, 9};

    uepm::mesh::vertex    v0{0, 0.0, 0.0, 0.0};
    uepm::mesh::vertex    v1{1, 0.5, 0.0, 0.0};
    uepm::mesh::vertex    v2{2, 0.0, 0.5, 0.0};
    uepm::mesh::element2d element{&v0, &v1, &v2};

    simulation.process_crossing({0.1, 0.1, 0.0}, {1.2, 0.1, 0.0}, uepm::PBMC::particle_type::electron, 2.0, &element);

    constexpr double elementary_charge_C = 1.602176634e-19;
    CHECK(simulation.get_number_electrons() == 0);
    CHECK(simulation.collected_contact_total_currents_A()[0] ==
          doctest::Approx(-2.0 * elementary_charge_C / options.m_time_step));
}

TEST_CASE("PBMC history exports named contact-flow current and cumulative charge columns") {
    uepm::PBMC::history_device_PBMC history;
    history.set_contact_flow_names({"drain"});
    history.add_data_to_history(1.0,
                                1,
                                0,
                                0.1,
                                0.0,
                                0.1,
                                0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                1.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                {},
                                {-2.0},
                                {0.5},
                                {-1.5},
                                {-3.0e-19},
                                {-0.5},
                                {-1.0},
                                {-1.0e-19},
                                {-2.0e-19});

    const auto filename = std::filesystem::temp_directory_path() / "ultimate_epm_contact_flow_history.csv";
    history.export_to_csv(filename.string());

    std::ifstream stream(filename);
    REQUIRE(stream.is_open());
    const std::string contents{std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
    CHECK(contents.find("collected_current_electron_drain_A,collected_current_hole_drain_A,") != std::string::npos);
    CHECK(contents.find("collected_current_drain_A,cumulative_collected_charge_drain_C") != std::string::npos);
    CHECK(contents.ends_with(",-2,0.5,-1.5,-3e-19,-0.5,-1,-1e-19,-2e-19\n"));
}

TEST_CASE("PBMC history averages every contact-current component over one trailing window") {
    uepm::PBMC::history_device_PBMC history;
    history.set_contact_flow_names({"drain"});
    history.set_contact_current_window_s(2.0);

    const auto add_contact_sample = [&](double time_s,
                                        double collected_electron_A,
                                        double collected_hole_A,
                                        double injected_A,
                                        double cumulative_collected_C,
                                        double cumulative_injected_C) {
        history.add_data_to_history(time_s,
                                    1,
                                    0,
                                    0.1,
                                    0.0,
                                    0.1,
                                    0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    1.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    {},
                                    {collected_electron_A},
                                    {collected_hole_A},
                                    {collected_electron_A + collected_hole_A},
                                    {cumulative_collected_C},
                                    {injected_A},
                                    {collected_electron_A + collected_hole_A - injected_A},
                                    {cumulative_injected_C},
                                    {cumulative_collected_C - cumulative_injected_C});
    };

    add_contact_sample(1.0, 2.0, 1.0, 0.5, 3.0, 0.5);
    add_contact_sample(2.0, 0.0, 0.0, 0.0, 3.0, 0.5);

    CHECK(history.m_list_collected_current_electron_A.back()[0] == doctest::Approx(1.0));
    CHECK(history.m_list_collected_current_hole_A.back()[0] == doctest::Approx(0.5));
    CHECK(history.m_list_collected_current_A.back()[0] == doctest::Approx(1.5));
    CHECK(history.m_list_injected_current_A.back()[0] == doctest::Approx(0.25));
    CHECK(history.m_list_net_contact_current_A.back()[0] == doctest::Approx(1.25));

    add_contact_sample(3.0, 4.0, 0.0, 1.0, 7.0, 1.5);

    CHECK(history.m_list_collected_current_electron_A.back()[0] == doctest::Approx(2.0));
    CHECK(history.m_list_collected_current_hole_A.back()[0] == doctest::Approx(0.0));
    CHECK(history.m_list_collected_current_A.back()[0] == doctest::Approx(2.0));
    CHECK(history.m_list_injected_current_A.back()[0] == doctest::Approx(0.5));
    CHECK(history.m_list_net_contact_current_A.back()[0] == doctest::Approx(1.5));
    CHECK(history.m_list_cumulative_collected_charge_C.back()[0] == doctest::Approx(7.0));
    CHECK(history.m_list_cumulative_injected_charge_C.back()[0] == doctest::Approx(1.5));
    CHECK(history.m_list_cumulative_net_contact_charge_C.back()[0] == doctest::Approx(5.5));
}

}  // namespace
