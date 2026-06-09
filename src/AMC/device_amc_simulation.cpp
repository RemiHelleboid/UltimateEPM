/**
 * @file device_amc_simulation.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-04
 *
 *
 */

#include "device_amc_simulation.hpp"

#include <fmt/chrono.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <fmt/xchar.h>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "physical_constants.hpp"
#include "unit_conversion.hpp"
#include "vtkWriter.hpp"

namespace uepm::amc {

namespace {

void write_vtk_time_collection(const std::string &pvd_filename, std::vector<vtk_time_series_record> records) {
    std::ofstream stream(pvd_filename);

    if (!stream.is_open()) {
        throw std::runtime_error(fmt::format("Could not open VTK collection file '{}'", pvd_filename));
    }

    std::sort(records.begin(), records.end(), [](const vtk_time_series_record &lhs, const vtk_time_series_record &rhs) {
        return lhs.m_time_s < rhs.m_time_s;
    });

    stream << std::setprecision(std::numeric_limits<double>::max_digits10);

    stream << "<?xml version=\"1.0\"?>\n";
    stream << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    stream << "  <Collection>\n";

    for (const auto &record : records) {
        stream << "    <DataSet timestep=\"" << record.m_time_s << "\" group=\"\" part=\"0\" file=\""
               << record.m_filename << "\"/>\n";
    }

    stream << "  </Collection>\n";
    stream << "</VTKFile>\n";
}

namespace {

std::string python_string_literal(const std::string &value) {
    std::string escaped;
    escaped.reserve(value.size() + 2);
    escaped.push_back('"');

    for (const char c : value) {
        switch (c) {
            case '\\':
                escaped += "\\\\";
                break;
            case '"':
                escaped += "\\\"";
                break;
            case '\n':
                escaped += "\\n";
                break;
            case '\r':
                escaped += "\\r";
                break;
            case '\t':
                escaped += "\\t";
                break;
            default:
                escaped.push_back(c);
                break;
        }
    }

    escaped.push_back('"');
    return escaped;
}

void write_paraview_scene_script(const std::filesystem::path &base_directory) {
    std::filesystem::create_directories(base_directory);

    const auto script_path = base_directory / "open_scene.py";
    const auto scene_dir   = std::filesystem::absolute(base_directory).lexically_normal();

    std::ofstream stream(script_path);

    if (!stream.is_open()) {
        throw std::runtime_error(fmt::format("Could not open ParaView scene script '{}'", script_path.string()));
    }

    stream << "from paraview.simple import *\n";
    stream << "import os\n\n";

    stream << "scene_dir = " << python_string_literal(scene_dir.generic_string()) << "\n";
    stream << "mesh_file = os.path.join(scene_dir, \"mesh\", \"mesh.pvd\")\n";
    stream << "particles_file = os.path.join(scene_dir, \"particles\", \"particles.pvd\")\n\n";

    stream << "mesh = OpenDataFile(mesh_file)\n";
    stream << "particles = OpenDataFile(particles_file)\n\n";

    stream << "view = GetActiveViewOrCreate(\"RenderView\")\n\n";

    stream << "mesh_display = Show(mesh, view)\n";
    stream << "mesh_display.Representation = \"Surface\"\n\n";

    stream << "particles_display = Show(particles, view)\n";
    stream << "particles_display.Representation = \"Point Gaussian\"\n";
    stream << "particles_display.PointSize = 6.0\n\n";

    stream << "ColorBy(particles_display, (\"POINTS\", \"energy_eV\"))\n";
    stream << "particles_display.RescaleTransferFunctionToDataRange(True, False)\n\n";

    stream << "ColorBy(mesh_display, (\"POINTS\", \"PoissonSolution_gradient\"))\n";
    stream << "mesh_display.RescaleTransferFunctionToDataRange(True, False)\n\n";

    stream << "animation_scene = GetAnimationScene()\n";
    stream << "animation_scene.UpdateAnimationUsingDataTimeSteps()\n";
    stream << "ResetCamera(view)\n";
    stream << "Render()\n";
}

}  // namespace

}  // namespace

struct impact_ionization_pair_seed {
    mesh::vector3 position;
    double        weight = 1.0;
};

vector3 device_amc_simulation::get_RamoUnitaryElectricField_at_position(const mesh::vector3 &position) const {
    if (m_state.m_use_constant_RamoUnitaryElectricField) {
        return m_state.m_RamoUnitaryElectricField_Vm_per_cm;
    }
    vector3 ramo_unitary_electric_field =
        m_device.interpolate_vector_at_location("RamoUnitaryPotential_gradient", position);
    return ramo_unitary_electric_field;
}

double device_amc_simulation::get_total_electron_weight() const {
    double total = 0.0;
    for (const auto &particle : m_list_particles) {
        if (particle->type() == particle_type::electron) {
            total += particle->weight();
        }
    }
    return total;
}

double device_amc_simulation::get_total_hole_weight() const {
    double total = 0.0;
    for (const auto &particle : m_list_particles) {
        if (particle->type() == particle_type::hole) {
            total += particle->weight();
        }
    }
    return total;
}

amc_transport_config device_amc_simulation::make_transport_config(const options_device_amc &options,
                                                                  particle_type             carrier_type) {
    amc_transport_config cfg;
    cfg.m_carrier_type                     = carrier_type;
    cfg.m_lattice_temperature              = options.m_lattice_temperature;
    cfg.m_max_energy_eV                    = options.m_max_energy_eV;
    cfg.m_self_scattering_safety_factor    = options.m_self_scattering_safety_factor;
    cfg.m_gamma_max_energy_samples         = options.m_gamma_max_energy_samples;
    cfg.m_enable_impact_ionization         = options.m_activate_impact_ionization;
    cfg.m_enable_impurity_scattering       = options.m_enable_impurity_scattering;
    cfg.m_impurity_density_source          = impurity_density_source::particle_local;
    cfg.m_background_impurity_density_cm_3 = 0.0;
    if (options.m_enable_impurity_scattering) {
        cfg.m_impurity_scattering_model = options.m_impurity_scattering_model;
    }
    return cfg;
}

amc_transport_kernel &device_amc_simulation::transport_for(particle_type type) {
    if (type == particle_type::electron) {
        return m_electron_transport;
    }

    return m_hole_transport;
}

const amc_transport_kernel &device_amc_simulation::transport_for(particle_type type) const {
    if (type == particle_type::electron) {
        return m_electron_transport;
    }

    return m_hole_transport;
}

void device_amc_simulation::initialize_particle_transport_state(particle_amc &particle) {
    auto &transport = transport_for(particle.type());

    if (transport.valleys().empty()) {
        throw std::runtime_error("transport kernel has no valleys/bands");
    }

    particle.state().valley_index = particle.index() % transport.valleys().size();

    transport.initialize_particle_state(particle);

    if (m_simulation_options.m_keep_particles_history) {
        particle.record_state();
    }
}

void device_amc_simulation::apply_z_periodicity_to_particles() {
    // No-op by default.
    // Only 2D self-consistent simulations override this.
}

void device_amc_simulation::initialize_scheduled_particle_injection() {
    m_state.m_scheduled_particle_injection_done = !m_simulation_options.m_enable_scheduled_particle_injection;
    if (!m_simulation_options.m_enable_scheduled_particle_injection) {
        return;
    }
    const auto &injection = m_simulation_options.m_scheduled_particle_injection;
    if (injection.m_time_s < 0.0) {
        throw std::invalid_argument("Scheduled particle injection time must be non-negative.");
    }
    if (injection.m_weight <= 0.0) {
        throw std::invalid_argument("Scheduled particle injection weight must be positive.");
    }
}

bool device_amc_simulation::has_pending_scheduled_particle_injection() const {
    if (!m_simulation_options.m_enable_scheduled_particle_injection) {
        return false;
    }
    if (m_state.m_scheduled_particle_injection_done) {
        return false;
    }
    return m_simulation_options.m_scheduled_particle_injection.m_time_s <= m_simulation_options.m_t_max;
}

void device_amc_simulation::inject_scheduled_particle_if_due() {
    if (!has_pending_scheduled_particle_injection()) {
        return;
    }

    const auto &injection = m_simulation_options.m_scheduled_particle_injection;

    const double dt = m_simulation_options.m_time_step;

    // Inject during the timestep that reaches the requested time.
    if (m_state.m_time_s + dt < injection.m_time_s) {
        return;
    }

    if (m_list_particles.size() >= m_simulation_options.m_max_number_particle) {
        fmt::print(stderr,
                   "Warning: scheduled particle injection skipped at t={:.6e} s "
                   "because max particle count was reached.\n",
                   injection.m_time_s);

        m_state.m_scheduled_particle_injection_done                = true;
        m_simulation_options.m_scheduled_particle_injection.m_done = true;
        return;
    }
    add_particle_at_position(injection.m_position_um, injection.m_particle_type, injection.m_weight);
    m_state.m_scheduled_particle_injection_done                = true;
    m_simulation_options.m_scheduled_particle_injection.m_done = true;
    fmt::print("\nScheduled particle injected at t={:.6e} s, "
               "position=({:.6e}, {:.6e}, {:.6e}) um, weight={:.6e}\n",
               injection.m_time_s,
               injection.m_position_um.x(),
               injection.m_position_um.y(),
               injection.m_position_um.z(),
               injection.m_weight);
}

device_amc_simulation::device_amc_simulation(const device::device     &simulation_device,
                                             const options_device_amc &simulation_option,
                                             int                       seed_random_generator)
    : m_device(simulation_device),
      m_electron_transport(make_transport_config(simulation_option, particle_type::electron), seed_random_generator),
      m_hole_transport(make_transport_config(simulation_option, particle_type::hole), seed_random_generator + 1),
      m_dimension(m_device.get_dimension()),
      m_simulation_options(simulation_option) {
    m_simulation_history.m_initial_seed_rng = seed_random_generator;
    m_electron_transport.initialize();
    m_hole_transport.initialize();
    initialize_scheduled_particle_injection();
}

device_amc_simulation::device_amc_simulation(const device::device     &device_simulation,
                                             const options_device_amc &simulation_option,
                                             const mesh::vector3      &starting_position,
                                             std::size_t               number_electrons_start,
                                             std::size_t               number_holes_start,
                                             int                       seed_random_generator)
    : m_device(device_simulation),
      m_electron_transport(make_transport_config(simulation_option, particle_type::electron), seed_random_generator),
      m_hole_transport(make_transport_config(simulation_option, particle_type::hole), seed_random_generator + 1),
      m_dimension(m_device.get_dimension()),
      m_simulation_options(simulation_option) {
    m_electron_transport.initialize();
    m_hole_transport.initialize();
    mesh::element *first_element{nullptr};
    if (m_dimension == 2) {
        first_element = m_device.find_element_at_location(starting_position.to_2d());
    } else {
        first_element = m_device.find_element_at_location(starting_position);
    }
    if (first_element == nullptr) {
        std::cout << "Error : particle can't find its first element. No particle created.    " << starting_position
                  << std::endl;
        return;
    }
    // Creation of electrons and then holes
    m_list_particles.reserve(number_electrons_start + number_holes_start);

    for (std::size_t i = 0; i < number_electrons_start; ++i) {
        const std::size_t particle_index = m_state.m_counter_particles_created++;
        m_list_particles.push_back(std::make_unique<particle_amc>(particle_index, particle_type::electron));
    }

    for (std::size_t i = 0; i < number_holes_start; ++i) {
        const std::size_t particle_index = m_state.m_counter_particles_created++;
        m_list_particles.push_back(std::make_unique<particle_amc>(particle_index, particle_type::hole));
    }

    //  Setup the initial position (unique), random number for RPLA and set containing elements for each particles.
    for (auto &p_particle : m_list_particles) {
        p_particle->set_position(starting_position);
        p_particle->set_weight(1.0);
        p_particle->set_containing_element(first_element);

        initialize_particle_transport_state(*p_particle);
    }

    m_simulation_history.m_initial_seed_rng = seed_random_generator;
    initialize_scheduled_particle_injection();
}

void device_amc_simulation::add_particle_at_position(const mesh::vector3 &location,
                                                     particle_type        type_of_particle,
                                                     double               weight) {
    mesh::element *first_element{nullptr};
    if (m_dimension == 2) {
        first_element = m_device.find_element_at_location(location.to_2d());
    } else {
        first_element = m_device.find_element_at_location(location);
    }
    if (first_element == nullptr) {
        std::cout << "Error : particle can't find its first element. No particle created.    " << location << std::endl;
        return;
    }
    const std::size_t idx_particle = m_state.m_counter_particles_created++;
    particle_state    initial_state{};
    initial_state.position = location;
    if (type_of_particle == particle_type::electron) {
        m_list_particles.push_back(
            std::make_unique<particle_amc>(idx_particle, particle_type::electron, initial_state, weight));
    } else {
        m_list_particles.push_back(
            std::make_unique<particle_amc>(idx_particle, particle_type::hole, initial_state, weight));
    }
    auto &particle = *m_list_particles.back();

    particle.set_containing_element(first_element);
    particle.set_weight(weight);

    initialize_particle_transport_state(particle);
}

void device_amc_simulation::add_particles_at_positions(const std::vector<mesh::vector3> &positions,
                                                       particle_type                     type_of_particle,
                                                       double                            weight) {
    // Reserve memory for all particles upfront
    m_list_particles.reserve(m_list_particles.size() + positions.size());
    for (const auto &location : positions) {
        mesh::element *first_element{nullptr};
        if (m_dimension == 2) {
            first_element = m_device.find_element_at_location(location.to_2d());
        } else {
            first_element = m_device.find_element_at_location(location);
        }

        if (first_element == nullptr) {
            std::cout << "Error : particle can't find its first element. No particle created.    " << location
                      << std::endl;
            continue;
        }

        const std::size_t idx_particle = m_state.m_counter_particles_created++;
        particle_state    initial_state{};
        initial_state.position = location;
        if (type_of_particle == particle_type::electron) {
            m_list_particles.push_back(
                std::make_unique<particle_amc>(idx_particle, particle_type::electron, initial_state, weight));
        } else {
            m_list_particles.push_back(
                std::make_unique<particle_amc>(idx_particle, particle_type::hole, initial_state, weight));
        }
        auto &particle = *m_list_particles.back();

        particle.set_containing_element(first_element);
        particle.set_weight(weight);

        initialize_particle_transport_state(particle);
    }
}

std::size_t device_amc_simulation::get_number_electrons() const {
    std::size_t nb_electron =
        std::accumulate(m_list_particles.begin(),
                        m_list_particles.end(),
                        0,
                        [](const std::size_t nb_part, const auto &p_part) {
                            return nb_part + static_cast<std::size_t>(p_part->type() == particle_type::electron);
                        });
    return nb_electron;
}

std::size_t device_amc_simulation::get_number_holes() const {
    std::size_t nb_hole =
        std::accumulate(m_list_particles.begin(),
                        m_list_particles.end(),
                        0,
                        [](const std::size_t nb_part, const auto &p_part_2) {
                            return nb_part + static_cast<std::size_t>(p_part_2->type() == particle_type::hole);
                        });
    return nb_hole;
}

double device_amc_simulation::ramo_current_scale_factor() const { return 1.0; }

// Compute Ramo current for a given single particle
double device_amc_simulation::compute_ramo_current_for_particle(const particle_amc &particle) const {
    const auto  *mesh_ptr     = m_device.get_p_mesh();
    const double scale_factor = ramo_current_scale_factor();
    auto         position     = particle.state().position;
    if (m_dimension == 2) {
        position.to_2d_inplace();
    }
    const auto weighting_field_m    = get_RamoUnitaryElectricField_at_position(position);
    double     current_contribution = scale_factor * particle.weight() * particle.get_signed_charge() *
                                      particle.state().velocity.dot(weighting_field_m);
    return current_contribution;
}

std::pair<double, double> device_amc_simulation::compute_ramo_current() const {
    double total_electron_current = 0.0;
    double total_hole_current     = 0.0;

    const auto  *mesh_ptr     = m_device.get_p_mesh();
    const double scale_factor = ramo_current_scale_factor();

    for (const auto &particle : m_list_particles) {
        auto position = particle->state().position;

        if (m_dimension == 2) {
            position.to_2d_inplace();
        }

        const auto weighting_field_m = get_RamoUnitaryElectricField_at_position(position);

        if (particle->type() == particle_type::electron) {
            total_electron_current += scale_factor * particle->weight() * particle->get_signed_charge() *
                                      particle->state().velocity.dot(weighting_field_m);
        } else {
            total_hole_current += scale_factor * particle->weight() * particle->get_signed_charge() *
                                  particle->state().velocity.dot(weighting_field_m);
        }
    }

    return std::make_pair(total_electron_current, total_hole_current);
}

void device_amc_simulation::transport_particles_one_time_step() {
    inject_scheduled_particle_if_due();
    const double                             dt = m_simulation_options.m_time_step;
    std::vector<impact_ionization_pair_seed> impact_pair_seeds;
    for (auto &p_particle : m_list_particles) {
        auto &particle = *p_particle;

        particle.state().previous_position = particle.state().position;
        particle.set_data_from_device(m_dimension);
        auto &transport = transport_for(particle.type());
        // Electric field is in V/cm, but we need it in V/m for the transport kernel, so we convert it here.
        constexpr double cm_to_m = 1.0e2;
        transport.drift_particle(particle, particle.state().electric_field * cm_to_m, dt);
    }
    apply_z_periodicity_to_particles();  // In 3D, this does nothing.
    update_element_and_check_boundary();
    remove_collected_particles();

    for (auto &p_particle : m_list_particles) {
        auto &particle = *p_particle;
        // Scattering
        auto      &transport = transport_for(particle.type());
        const auto event     = transport.scatter_particle(particle, dt);

        if (event == scattering_event::impact_ionization) {
            m_simulation_history.m_last_impact_ionization_position = particle.state().position;
            m_simulation_history.m_impact_ionization_positions.push_back(particle.state().position);
            bool enable_part_creation = m_simulation_options.m_particle_creation_activated &&
                                        m_simulation_options.m_scheduled_particle_injection.m_done;
            if (enable_part_creation) {
                const std::size_t queued_particles = 2 * impact_pair_seeds.size();

                if (m_list_particles.size() + queued_particles + 2 <= m_simulation_options.m_max_number_particle) {
                    impact_pair_seeds.push_back(impact_ionization_pair_seed{.position = particle.state().position,
                                                                            .weight   = particle.weight()});
                }
            }
        }

        if (m_simulation_options.m_keep_particles_history) {
            particle.record_state();
        }
    }

    for (const auto &seed : impact_pair_seeds) {
        add_particle_at_position(seed.position, particle_type::electron, seed.weight);
        add_particle_at_position(seed.position, particle_type::hole, seed.weight);
    }

    if (!impact_pair_seeds.empty()) {
        apply_z_periodicity_to_particles();
        update_element_and_check_boundary();
        remove_collected_particles();
    }
}

void device_amc_simulation::advance_particles_one_time_step() {
    transport_particles_one_time_step();
    m_state.m_time_s += m_simulation_options.m_time_step;
    ++m_state.m_iteration;
}

void device_amc_simulation::set_particles_transport_data_from_device() {
    for (auto &p_particle : m_list_particles) {
        p_particle->set_data_from_device(m_dimension);
    }
}

void reflect_particle_to_previous_position(particle_amc &particle) {
    auto &state = particle.state();

    state.position = state.previous_position;

    state.velocity *= -1.0;
    state.local_k *= -1.0;
}

void device_amc_simulation::update_element_and_check_boundary() {
    const bool is_2d = m_dimension == 2;
    for (auto &p_particle : m_list_particles) {
        auto                &particle    = *p_particle;
        const mesh::element *old_element = particle.get_containing_element();
        if (old_element == nullptr) {
            particle.set_crossed_contact(true);
            continue;
        }
        mesh::vector3 current_position = particle.state().position;
        if (is_2d) {
            current_position.to_2d_inplace();
        }
        if (old_element->is_location_inside_element(current_position)) {
            continue;
        }
        if (m_device.check_enters_contact(current_position)) {
            particle.set_crossed_contact(true);
            continue;
        }
        auto *new_element = m_device.find_element_at_location(current_position);
        if (new_element == nullptr) {
            reflect_particle_to_previous_position(particle);
            continue;
        }
        if (m_device.get_material_name_at_element(new_element) != "Silicon") {
            reflect_particle_to_previous_position(particle);
            continue;
        }
        particle.set_containing_element(new_element);
    }
}

void device_amc_simulation::remove_collected_particles() {
    bool remove_particle = false;
    int  nb_part_erased  = 0;
    for (const auto &p_particle : m_list_particles) {
        if (p_particle->state().m_crossed_contact) {
            remove_particle = true;
            nb_part_erased++;
        }
    }
    if (remove_particle) {
        std::erase_if(m_list_particles, [](auto &&p_part) { return p_part->state().m_crossed_contact; });
    }
    std::vector<double> currents = m_device.get_electrode_currents();
    // CHECK SIZE OF CURRENTS VECTOR
    if (currents.size() < 2) {
        throw std::runtime_error("Error: currents vector should have at least 2 elements (anode and cathode currents)");
    }
}

void device_amc_simulation::run() {
    while (m_state.m_time_s < m_simulation_options.m_t_max) {
        if (m_list_particles.empty()) {
            break;
        }
        if (m_simulation_options.m_stop_simu_when_no_electron_remaining && get_number_electrons() == 0) {
            break;
        }
        if (has_reached_avalanche()) {
            break;
        }
        advance_particles_one_time_step();
        const auto nb_electrons         = get_number_electrons();
        const auto nb_holes             = get_number_holes();
        const auto nb_impact_ionization = m_simulation_history.m_impact_ionization_positions.size();

        double dumb_ramo_current_e_h_total = 0.0;
        double dumb_0 = 0.0;
        m_simulation_history.add_data_to_history(m_state.m_time_s,
                                                 nb_electrons,
                                                 nb_holes,
                                                 nb_impact_ionization,
                                                 dumb_ramo_current_e_h_total,
                                                 dumb_ramo_current_e_h_total,
                                                 dumb_ramo_current_e_h_total,
                                                 0.0,
                                                 dumb_0,
                                                 dumb_0,
                                                 dumb_0,
                                                 dumb_0,
                                                 dumb_0,
                                                 dumb_0);

        if (m_simulation_options.m_export_time_step &&
            m_state.m_iteration % static_cast<std::size_t>(m_simulation_options.m_frequency_export_trajectory) == 0) {
            export_current_time_step_as_csv(m_simulation_options.m_prefix_export_filename);
        }
    }
}

std::vector<mesh::vector3> device_amc_simulation::get_all_particles_position() const {
    std::vector<mesh::vector3> all_positions(m_list_particles.size());
    std::transform(m_list_particles.begin(), m_list_particles.end(), all_positions.begin(), [](auto &&p_particle) {
        return p_particle->state().position;
    });
    return all_positions;
}

void device_amc_simulation::export_current_snapshot() const {
    const std::filesystem::path base_directory(m_simulation_options.m_prefix_export_filename);

    const std::filesystem::path mesh_directory      = base_directory / "mesh";
    const std::filesystem::path particles_directory = base_directory / "particles";

    export_current_mesh_as_vtk(mesh_directory.string());
    export_current_particles_as_vtp(particles_directory.string());

    write_paraview_scene_script(base_directory);
}

std::pair<double, double> device_amc_simulation::compute_depletion_region() const {
    double x_min = std::numeric_limits<double>::max();
    double x_max = std::numeric_limits<double>::max();
    // Compute the distance of the particles to the center of the device. xmin is the maximal distance to the center
    // towards the anode and xmax is the maximal distance to the center towards the cathode.
    double center_x = m_device.get_p_mesh()->get_bounding_box().get_x_min() +
                      0.5 * (m_device.get_p_mesh()->get_bounding_box().get_x_max() -
                             m_device.get_p_mesh()->get_bounding_box().get_x_min());
    for (const auto &p_particle : m_list_particles) {
        double distance_to_center = p_particle->state().position.x() - center_x;
        if (distance_to_center < 0) {
            x_min = std::min(std::abs(x_min), std::abs(distance_to_center));
        } else {
            x_max = std::min(x_max, distance_to_center);
        }
    }
    x_min = center_x - std::abs(x_min);
    x_max = center_x + std::abs(x_max);
    return {x_min, x_max};
}

std::string device_amc_simulation::initialize_simulation_history_file()  {
    std::string simulation_name_for_file =
        m_simulation_options.m_simulation_name.empty() ? "simulation" : m_simulation_options.m_simulation_name;
    const std::string history_filename =
        fmt::format("{}/{}_history.csv", m_simulation_options.m_output_directory, simulation_name_for_file);
    // Create output directory if it doesn't exist
    fmt::print("Initializing simulation history file at '{}'\n", history_filename);
    std::filesystem::create_directories(m_simulation_options.m_output_directory);
    std::ofstream stream(history_filename);
    if (!stream.is_open()) {
        throw std::runtime_error(fmt::format("Could not open simulation history CSV file '{}'", history_filename));
    }
    stream.close();
    m_simulation_history.print_header_csv(history_filename);
    return history_filename;
}

void device_amc_simulation::export_current_time_step_as_csv(const std::string &prefix_filename) const {
    const std::string iteration_filename = fmt::format("{}.{:012d}.csv", prefix_filename, m_state.m_iteration);

    std::ofstream stream(iteration_filename);

    if (!stream.is_open()) {
        throw std::runtime_error(fmt::format("Could not open particle CSV file '{}'", iteration_filename));
    }

    stream << std::setprecision(std::numeric_limits<double>::max_digits10);

    stream << "particle_index,"
           << "type,"
           << "time_s,"
           << "x_um,"
           << "y_um,"
           << "z_um,"
           << "local_kx_1_per_m,"
           << "local_ky_1_per_m,"
           << "local_kz_1_per_m,"
           << "vx_m_per_s,"
           << "vy_m_per_s,"
           << "vz_m_per_s,"
           << "energy_eV,"
           << "electric_field_x_V_per_cm,"
           << "electric_field_y_V_per_cm,"
           << "electric_field_z_V_per_cm,"
           << "electric_field_norm_V_per_cm,"
           << "weight,"
           << "valley_index,"
           << "signed_charge_C\n";

    for (const auto &p_particle : m_list_particles) {
        const auto &particle = *p_particle;
        const auto &state    = particle.state();

        const auto &position       = state.position;
        const auto &local_k        = state.local_k;
        const auto &velocity       = state.velocity;
        const auto &electric_field = state.electric_field;

        stream << particle.index() << ',' << static_cast<int>(particle.type()) << ',' << state.time << ','
               << position.x() << ',' << position.y() << ',' << position.z() << ',' << local_k.x() << ',' << local_k.y()
               << ',' << local_k.z() << ',' << velocity.x() << ',' << velocity.y() << ',' << velocity.z() << ','
               << state.kinetic_energy << ',' << electric_field.x() << ',' << electric_field.y() << ','
               << electric_field.z() << ',' << electric_field.norm() << ',' << particle.weight() << ','
               << state.valley_index << ',' << particle.get_signed_charge() << '\n';
    }
}
void device_amc_simulation::export_all_trajectories_as_csv(const std::string &prefix_filename) const {
    for (auto &&p_particle : m_list_particles) {
        // std::cout << "\rExporting trajectory of particle " << p_particle->get_index() << std::flush;
        std::string filename = fmt::format("{}particle_{:06}_trajectory.csv", prefix_filename, p_particle->index());
        p_particle->export_trajectory_as_csv(filename);
    }
    // std::cout << std::endl;
}

void device_amc_simulation::export_current_particles_as_vtp(const std::string &directory) const {
    const std::filesystem::path output_directory(directory);
    std::filesystem::create_directories(output_directory);

    const std::string filename = fmt::format("particles_{:012d}.vtp", m_state.m_iteration);

    const std::filesystem::path vtp_path = output_directory / filename;
    const std::filesystem::path pvd_path = output_directory / "particles.pvd";

    std::ofstream stream(vtp_path);

    if (!stream.is_open()) {
        throw std::runtime_error(fmt::format("Could not open particle VTP file '{}'", vtp_path.string()));
    }

    const std::size_t number_particles = m_list_particles.size();

    stream << std::setprecision(std::numeric_limits<double>::max_digits10);

    stream << "<?xml version=\"1.0\"?>\n";
    stream << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    stream << "  <PolyData>\n";
    stream << "    <Piece NumberOfPoints=\"" << number_particles << "\" NumberOfVerts=\"" << number_particles
           << "\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";

    stream << "      <PointData Scalars=\"energy_eV\" Vectors=\"velocity_m_per_s\">\n";

    stream << "        <DataArray type=\"Int64\" Name=\"particle_index\" format=\"ascii\">\n";
    stream << "          ";
    for (const auto &p_particle : m_list_particles) {
        stream << static_cast<long long>(p_particle->index()) << ' ';
    }
    stream << "\n";
    stream << "        </DataArray>\n";

    stream << "        <DataArray type=\"Int32\" Name=\"particle_type\" format=\"ascii\">\n";
    stream << "          ";
    for (const auto &p_particle : m_list_particles) {
        stream << static_cast<int>(p_particle->type()) << ' ';
    }
    stream << "\n";
    stream << "        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"time_s\" format=\"ascii\">\n";
    stream << "          ";
    for (const auto &p_particle : m_list_particles) {
        stream << p_particle->state().time << ' ';
    }
    stream << "\n";
    stream << "        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"energy_eV\" format=\"ascii\">\n";
    stream << "          ";
    for (const auto &p_particle : m_list_particles) {
        stream << p_particle->state().kinetic_energy << ' ';
    }
    stream << "\n";
    stream << "        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"weight\" format=\"ascii\">\n";
    stream << "          ";
    for (const auto &p_particle : m_list_particles) {
        stream << p_particle->weight() << ' ';
    }
    stream << "\n";
    stream << "        </DataArray>\n";

    stream << "        <DataArray type=\"Int32\" Name=\"valley_index\" format=\"ascii\">\n";
    stream << "          ";
    for (const auto &p_particle : m_list_particles) {
        stream << static_cast<int>(p_particle->state().valley_index) << ' ';
    }
    stream << "\n";
    stream << "        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"signed_charge_C\" format=\"ascii\">\n";
    stream << "          ";
    for (const auto &p_particle : m_list_particles) {
        stream << p_particle->get_signed_charge() << ' ';
    }
    stream << "\n";
    stream << "        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"electric_field_norm_V_per_cm\" format=\"ascii\">\n";
    stream << "          ";
    for (const auto &p_particle : m_list_particles) {
        stream << p_particle->state().electric_field.norm() << ' ';
    }
    stream << "\n";
    stream << "        </DataArray>\n";

    stream
        << "        <DataArray type=\"Float64\" Name=\"velocity_m_per_s\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    stream << "          ";
    for (const auto &p_particle : m_list_particles) {
        const auto &velocity = p_particle->state().velocity;
        stream << velocity.x() << ' ' << velocity.y() << ' ' << velocity.z() << ' ';
    }
    stream << "\n";
    stream << "        </DataArray>\n";

    stream
        << "        <DataArray type=\"Float64\" Name=\"local_k_1_per_m\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    stream << "          ";
    for (const auto &p_particle : m_list_particles) {
        const auto &local_k = p_particle->state().local_k;
        stream << local_k.x() << ' ' << local_k.y() << ' ' << local_k.z() << ' ';
    }
    stream << "\n";
    stream << "        </DataArray>\n";

    stream << "        <DataArray type=\"Float64\" Name=\"electric_field_V_per_cm\" NumberOfComponents=\"3\" "
              "format=\"ascii\">\n";
    stream << "          ";
    for (const auto &p_particle : m_list_particles) {
        const auto &electric_field = p_particle->state().electric_field;
        stream << electric_field.x() << ' ' << electric_field.y() << ' ' << electric_field.z() << ' ';
    }
    stream << "\n";
    stream << "        </DataArray>\n";

    stream << "      </PointData>\n";

    stream << "      <Points>\n";
    stream << "        <DataArray type=\"Float64\" Name=\"position_um\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    stream << "          ";
    for (const auto &p_particle : m_list_particles) {
        const auto &position = p_particle->state().position;
        stream << position.x() << ' ' << position.y() << ' ' << position.z() << ' ';
    }
    stream << "\n";
    stream << "        </DataArray>\n";
    stream << "      </Points>\n";

    stream << "      <Verts>\n";

    stream << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    stream << "          ";
    for (std::size_t i = 0; i < number_particles; ++i) {
        stream << static_cast<long long>(i) << ' ';
    }
    stream << "\n";
    stream << "        </DataArray>\n";

    stream << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    stream << "          ";
    for (std::size_t i = 0; i < number_particles; ++i) {
        stream << static_cast<long long>(i + 1) << ' ';
    }
    stream << "\n";
    stream << "        </DataArray>\n";

    stream << "      </Verts>\n";

    stream << "    </Piece>\n";
    stream << "  </PolyData>\n";
    stream << "</VTKFile>\n";

    const auto already_recorded =
        std::find_if(m_particle_vtp_export_records.begin(),
                     m_particle_vtp_export_records.end(),
                     [&](const vtk_time_series_record &record) { return record.m_filename == filename; });

    if (already_recorded == m_particle_vtp_export_records.end()) {
        m_particle_vtp_export_records.push_back(
            vtk_time_series_record{.m_time_s = m_state.m_time_s, .m_filename = filename});
    }

    write_particle_vtp_time_collection(pvd_path.string());
}

void device_amc_simulation::write_particle_vtp_time_collection(const std::string &pvd_filename) const {
    write_vtk_time_collection(pvd_filename, m_particle_vtp_export_records);
}

void device_amc_simulation::export_current_mesh_as_vtk(const std::string &directory) const {
    const std::filesystem::path output_directory(directory);
    std::filesystem::create_directories(output_directory);

    const std::string filename = fmt::format("mesh_{:012d}.vtu", m_state.m_iteration);

    const std::filesystem::path vtu_path = output_directory / filename;
    const std::filesystem::path pvd_path = output_directory / "mesh.pvd";

    file::export_as_vtu(*(m_device.get_p_mesh()), vtu_path.string(), {}, {}, true);

    const auto already_recorded =
        std::find_if(m_mesh_vtk_export_records.begin(),
                     m_mesh_vtk_export_records.end(),
                     [&](const vtk_time_series_record &record) { return record.m_filename == filename; });

    if (already_recorded == m_mesh_vtk_export_records.end()) {
        m_mesh_vtk_export_records.push_back(
            vtk_time_series_record{.m_time_s = m_state.m_time_s, .m_filename = filename});
    }

    write_mesh_vtk_time_collection(pvd_path.string());
}

void device_amc_simulation::write_mesh_vtk_time_collection(const std::string &pvd_filename) const {
    write_vtk_time_collection(pvd_filename, m_mesh_vtk_export_records);
}

}  // namespace uepm::amc
