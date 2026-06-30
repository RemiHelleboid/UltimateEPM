/**
 * @file device_pbmc_simulation.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-04
 *
 *
 */

#include "device_pbmc_simulation.hpp"

#include <fmt/chrono.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <fmt/xchar.h>
#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "physical_constants.hpp"
#include "unit_conversion.hpp"
#include "vtkWriter.hpp"

namespace uepm::PBMC {

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

    stream << "view = GetActiveViewOrCreate(\"RenderView\")\n";
    stream << "view.InteractionMode = \"2D\"\n";
    stream << "view.Background = [0.0, 0.0, 0.0]\n";
    stream << "view.OrientationAxesVisibility = 1\n\n";

    stream << "mesh_display = Show(mesh, view)\n";
    stream << "mesh_display.Representation = \"Surface\"\n";
    stream << "ColorBy(mesh_display, (\"POINTS\", \"PoissonSolution\"))\n";
    stream << "mesh_display.RescaleTransferFunctionToDataRange(True, False)\n\n";

    stream << "poisson_lut = GetColorTransferFunction(\"PoissonSolution\")\n";
    stream << "poisson_lut.ApplyPreset(\"Jet\", True)\n";
    stream << "poisson_lut.RescaleTransferFunction(-0.18, 1.0)\n";
    stream << "poisson_lut.ColorSpace = \"RGB\"\n\n";

    stream << "poisson_bar = GetScalarBar(poisson_lut, view)\n";
    stream << "poisson_bar.Title = \"PoissonSolution\"\n";
    stream << "poisson_bar.ComponentTitle = \"\"\n";
    stream << "poisson_bar.Visibility = 1\n";
    stream << "poisson_bar.WindowLocation = \"Upper Right Corner\"\n";
    stream << "poisson_bar.TitleColor = [1.0, 1.0, 1.0]\n";
    stream << "poisson_bar.LabelColor = [1.0, 1.0, 1.0]\n";
    stream << "mesh_display.SetScalarBarVisibility(view, True)\n\n";

    stream << "particles_display = Show(particles, view)\n";
    stream << "particles_display.Representation = \"Point Gaussian\"\n";
    stream << "particles_display.GaussianRadius = 0.001\n";
    stream << "particles_display.Opacity = 1.0\n";
    stream << "particles_display.RenderPointsAsSpheres = 1\n";
    stream << "ColorBy(particles_display, (\"POINTS\", \"particle_type\"))\n\n";

    stream << "particle_lut = GetColorTransferFunction(\"particle_type\")\n";
    stream << "particle_lut.InterpretValuesAsCategories = 1\n";
    stream << "particle_lut.AnnotationsInitialized = 1\n";
    stream << "particle_lut.Annotations = [\"0\", \"electron\", \"1\", \"hole\"]\n";
    stream << "particle_lut.IndexedColors = [0.25, 0.0, 1.0, 1.0, 0.0, 0.0]\n\n";

    stream << "particles_display.LookupTable = particle_lut\n";
    stream << "particles_display.SetScalarBarVisibility(view, True)\n\n";

    stream << "particle_bar = GetScalarBar(particle_lut, view)\n";
    stream << "particle_bar.Title = \"particle_type\"\n";
    stream << "particle_bar.ComponentTitle = \"\"\n";
    stream << "particle_bar.Visibility = 1\n";
    stream << "particle_bar.WindowLocation = \"Lower Right Corner\"\n";
    stream << "particle_bar.TitleColor = [1.0, 1.0, 1.0]\n";
    stream << "particle_bar.LabelColor = [1.0, 1.0, 1.0]\n\n";

    stream << "animation_scene = GetAnimationScene()\n";
    stream << "animation_scene.UpdateAnimationUsingDataTimeSteps()\n\n";

    stream << "time_keeper = GetTimeKeeper()\n";
    stream << "AnnotateTimeFilter1 = AnnotateTimeFilter(Input=mesh)\n";
    stream << "AnnotateTimeFilter1.Format = \"time: {time:.3e} s\"\n";
    stream << "time_display = Show(AnnotateTimeFilter1, view)\n";
    stream << "time_display.FontSize = 18\n";
    stream << "time_display.Color = [1.0, 1.0, 1.0]\n";
    stream << "time_display.WindowLocation = \"Upper Center\"\n\n";

    stream << "mesh.UpdatePipeline()\n";
    stream << "particles.UpdatePipeline()\n\n";

    stream << "particles_display.Visibility = 0\n";
    stream << "ResetCamera(view)\n";
    stream << "particles_display.Visibility = 1\n\n";

    stream << "view.CameraPosition = [0.0, 0.0, 1.0]\n";
    stream << "view.CameraFocalPoint = [0.0, 0.0, 0.0]\n";
    stream << "view.CameraViewUp = [0.0, 1.0, 0.0]\n";
    stream << "view.CameraParallelScale *= 0.75\n\n";

    stream << "SetActiveSource(mesh)\n";
    stream << "Render()\n";
}
}  // namespace

}  // namespace

struct impact_ionization_pair_seed {
    mesh::vector3 position;
    double        weight = 1.0;
};

std::vector<std::string> split_csv_line(const std::string &line) {
    std::vector<std::string> fields;
    std::stringstream        stream(line);
    std::string              field;
    while (std::getline(stream, field, ',')) {
        fields.push_back(field);
    }
    if (!line.empty() && line.back() == ',') {
        fields.emplace_back();
    }
    return fields;
}

std::unordered_map<std::string, std::size_t> csv_header_index(const std::string &header_line) {
    std::unordered_map<std::string, std::size_t> indices;
    const auto                                   columns = split_csv_line(header_line);
    for (std::size_t i = 0; i < columns.size(); ++i) {
        indices.emplace(columns[i], i);
    }
    return indices;
}

const std::string &required_csv_field(const std::vector<std::string>                    &fields,
                                      const std::unordered_map<std::string, std::size_t> &indices,
                                      const std::string                                 &name,
                                      std::size_t                                        line_number) {
    const auto it = indices.find(name);
    if (it == indices.end()) {
        throw std::invalid_argument(fmt::format("Particle state CSV is missing required column '{}'.", name));
    }
    if (it->second >= fields.size()) {
        throw std::invalid_argument(fmt::format("Particle state CSV line {} is missing value '{}'.", line_number, name));
    }
    return fields[it->second];
}

double required_csv_double(const std::vector<std::string>                    &fields,
                           const std::unordered_map<std::string, std::size_t> &indices,
                           const std::string                                 &name,
                           std::size_t                                        line_number) {
    return std::stod(required_csv_field(fields, indices, name, line_number));
}

std::size_t required_csv_size(const std::vector<std::string>                    &fields,
                              const std::unordered_map<std::string, std::size_t> &indices,
                              const std::string                                 &name,
                              std::size_t                                        line_number) {
    return static_cast<std::size_t>(std::stoull(required_csv_field(fields, indices, name, line_number)));
}

double optional_csv_double(const std::vector<std::string>                    &fields,
                           const std::unordered_map<std::string, std::size_t> &indices,
                           const std::string                                 &name,
                           double                                             fallback) {
    const auto it = indices.find(name);
    if (it == indices.end() || it->second >= fields.size() || fields[it->second].empty()) {
        return fallback;
    }
    return std::stod(fields[it->second]);
}

std::size_t optional_csv_size(const std::vector<std::string>                    &fields,
                              const std::unordered_map<std::string, std::size_t> &indices,
                              const std::string                                 &name,
                              std::size_t                                        fallback) {
    const auto it = indices.find(name);
    if (it == indices.end() || it->second >= fields.size() || fields[it->second].empty()) {
        return fallback;
    }
    return static_cast<std::size_t>(std::stoull(fields[it->second]));
}

particle_type parse_particle_type_field(const std::string &value) {
    if (value == "0" || value == "electron") {
        return particle_type::electron;
    }
    if (value == "1" || value == "hole") {
        return particle_type::hole;
    }
    throw std::invalid_argument(fmt::format("Invalid particle type '{}'. Expected 0/electron or 1/hole.", value));
}

void options_device_PBMC::validate() const {
    if (m_time_step <= 0.0) {
        throw std::invalid_argument("--dt must be positive.");
    }
    if (m_t_max <= 0.0) {
        throw std::invalid_argument("--time must be positive.");
    }
    if (m_lattice_temperature < 0.0) {
        throw std::invalid_argument("--temperature must be non-negative.");
    }
    if (m_max_energy_eV <= 0.0) {
        throw std::invalid_argument("--max-energy must be positive.");
    }
    if (m_self_scattering_safety_factor <= 0.0) {
        throw std::invalid_argument("--gamma-safety must be positive.");
    }
    if (m_gamma_max_energy_samples < 2) {
        throw std::invalid_argument("--gamma-samples must be at least 2.");
    }
    if (m_max_number_particle == 0) {
        throw std::invalid_argument("--max-particles must be positive.");
    }
    if (m_frequency_export_trajectory <= 0) {
        throw std::invalid_argument("--export-frequency must be positive.");
    }
    if (m_nb_threads <= 0) {
        throw std::invalid_argument("--nthreads must be positive.");
    }
    if (m_enable_scheduled_particle_injection) {
        const auto &injection = m_scheduled_particle_injection;
        if (injection.m_time_s < 0.0) {
            throw std::invalid_argument("--inject-time must be non-negative.");
        }
        if (injection.m_time_s > m_t_max) {
            throw std::invalid_argument("--inject-time cannot be larger than --time.");
        }
        if (injection.m_weight <= 0.0) {
            throw std::invalid_argument("--inject-weight must be positive.");
        }
    }
}

vector3 device_pbmc_simulation::get_RamoUnitaryElectricField_at_position(const mesh::vector3 &position) const {
    return get_RamoUnitaryElectricField_at_position(position, nullptr);
}

vector3 device_pbmc_simulation::get_RamoUnitaryElectricField_at_position(
    const mesh::vector3 &position,
    const mesh::element *containing_element) const {
    if (m_state.m_use_constant_RamoUnitaryElectricField) {
        return m_state.m_RamoUnitaryElectricField_Vm_per_cm;
    }
    if (containing_element != nullptr) {
        return containing_element->interpolate_vector_at_location("RamoUnitaryPotential_gradient", position);
    }
    vector3 ramo_unitary_electric_field = m_device.interpolate_vector_at_location("RamoUnitaryPotential_gradient",
                                                                                  position);
    return ramo_unitary_electric_field;
}

double device_pbmc_simulation::get_total_electron_weight() const {
    double total = 0.0;
    for (const auto &particle : m_list_particles) {
        if (particle->type() == particle_type::electron) {
            total += particle->weight();
        }
    }
    return total;
}

double device_pbmc_simulation::get_total_hole_weight() const {
    double total = 0.0;
    for (const auto &particle : m_list_particles) {
        if (particle->type() == particle_type::hole) {
            total += particle->weight();
        }
    }
    return total;
}

pbmc_transport_config device_pbmc_simulation::make_transport_config(const options_device_PBMC &options,
                                                                    particle_type              carrier_type) {
    pbmc_transport_config cfg;
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
        cfg.m_impurity_screening_model  = options.m_impurity_screening_model;
    }
    return cfg;
}

pbmc_transport_kernel &device_pbmc_simulation::transport_for(particle_type type) {
    if (type == particle_type::electron) {
        return m_electron_transport;
    }

    return m_hole_transport;
}

const pbmc_transport_kernel &device_pbmc_simulation::transport_for(particle_type type) const {
    if (type == particle_type::electron) {
        return m_electron_transport;
    }

    return m_hole_transport;
}

void device_pbmc_simulation::initialize_thread_transports(int seed_random_generator) {
    const auto number_threads = static_cast<std::size_t>(m_simulation_options.m_nb_threads);
    if (number_threads <= 1) {
        return;
    }

    m_thread_electron_transports.assign(number_threads, m_electron_transport);
    m_thread_hole_transports.assign(number_threads, m_hole_transport);

    for (std::size_t thread_index = 0; thread_index < number_threads; ++thread_index) {
        const auto offset = static_cast<std::uint64_t>(7919 * thread_index);
        m_thread_electron_transports[thread_index].seed(static_cast<std::uint64_t>(seed_random_generator) + offset);
        m_thread_hole_transports[thread_index].seed(static_cast<std::uint64_t>(seed_random_generator) + offset + 1);
    }
}

pbmc_transport_kernel &device_pbmc_simulation::transport_for(particle_type type, std::size_t thread_index) {
    if (m_thread_electron_transports.empty()) {
        return transport_for(type);
    }

    if (type == particle_type::electron) {
        return m_thread_electron_transports.at(thread_index);
    }
    return m_thread_hole_transports.at(thread_index);
}

void device_pbmc_simulation::initialize_particle_transport_state(pbmc_particle &particle) {
    auto &transport = transport_for(particle.type());

    if (transport.valleys().empty()) {
        throw std::runtime_error("transport kernel has no valleys/bands");
    }

    particle.state().valley_index = particle.index() % transport.valleys().size();

    particle.set_data_from_device(m_dimension, m_simulation_options.m_enable_impurity_scattering);
    transport.initialize_particle_state(particle, particle.get_lattice_temperature());

    if (m_simulation_options.m_keep_particles_history) {
        particle.record_state();
    }
}

void device_pbmc_simulation::flatten_particle_positions_for_2d() {
    if (m_dimension != 2) {
        return;
    }
    for (auto &p_particle : m_list_particles) {
        p_particle->state().position.set_z(0.0);
        p_particle->state().previous_position.set_z(0.0);
    }
}

void device_pbmc_simulation::initialize_scheduled_particle_injection() {
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

bool device_pbmc_simulation::has_pending_scheduled_particle_injection() const {
    if (!m_simulation_options.m_enable_scheduled_particle_injection) {
        return false;
    }
    if (m_state.m_scheduled_particle_injection_done) {
        return false;
    }
    return m_simulation_options.m_scheduled_particle_injection.m_time_s <= m_simulation_options.m_t_max;
}

void device_pbmc_simulation::validate_time_step_against_scattering_rate() const {
    const double gamma_max_s_1 = std::max(m_electron_transport.gamma_max(), m_hole_transport.gamma_max());
    if (gamma_max_s_1 <= 0.0) {
        throw std::invalid_argument("PBMC gamma_max must be positive before validating the device time step.");
    }

    const double     scattering_step_ratio = m_simulation_options.m_time_step * gamma_max_s_1;
    constexpr double warning_ratio         = 0.05;
    constexpr double error_ratio           = 0.10;

    if (scattering_step_ratio > error_ratio) {
        throw std::invalid_argument(fmt::format(
            "PBMC device time step is too large for fixed-step scattering: dt * gamma_max = {:.6e} "
            "(dt = {:.6e} s, gamma_max = {:.6e} s^-1). Device PBMC applies at most one scattering event per "
            "fixed time step, so use dt <= {:.6e} s or increase gamma/max-energy settings only if gamma_max is "
            "underestimated.",
            scattering_step_ratio,
            m_simulation_options.m_time_step,
            gamma_max_s_1,
            error_ratio / gamma_max_s_1));
    }

    if (scattering_step_ratio > warning_ratio) {
        fmt::print(stderr,
                   "WARNING: PBMC device time step is close to the fixed-step scattering limit: "
                   "dt * gamma_max = {:.6e} (dt = {:.6e} s, gamma_max = {:.6e} s^-1). "
                   "For better accuracy, prefer dt <= {:.6e} s.\n",
                   scattering_step_ratio,
                   m_simulation_options.m_time_step,
                   gamma_max_s_1,
                   warning_ratio / gamma_max_s_1);
    }
}

void device_pbmc_simulation::inject_scheduled_particle_if_due() {
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

device_pbmc_simulation::device_pbmc_simulation(const device::device      &simulation_device,
                                               const options_device_PBMC &simulation_option,
                                               int                        seed_random_generator)
    : m_device(simulation_device),
      m_electron_transport(make_transport_config(simulation_option, particle_type::electron),
                           simulation_option.m_material_model,
                           seed_random_generator),
      m_hole_transport(make_transport_config(simulation_option, particle_type::hole),
                       simulation_option.m_material_model,
                       seed_random_generator + 1),
      m_dimension(m_device.get_dimension()),
      m_simulation_options(simulation_option),
      m_boundary_reflection_rng(static_cast<unsigned int>(seed_random_generator + 2)) {
    m_simulation_history.m_initial_seed_rng = seed_random_generator;
    m_electron_transport.initialize();
    m_hole_transport.initialize();
    validate_time_step_against_scattering_rate();
    initialize_thread_transports(seed_random_generator);
    initialize_scheduled_particle_injection();
}

device_pbmc_simulation::device_pbmc_simulation(const device::device      &device_simulation,
                                               const options_device_PBMC &simulation_option,
                                               const mesh::vector3       &starting_position,
                                               std::size_t                number_electrons_start,
                                               std::size_t                number_holes_start,
                                               int                        seed_random_generator)
    : m_device(device_simulation),
      m_electron_transport(make_transport_config(simulation_option, particle_type::electron),
                           simulation_option.m_material_model,
                           seed_random_generator),
      m_hole_transport(make_transport_config(simulation_option, particle_type::hole),
                       simulation_option.m_material_model,
                       seed_random_generator + 1),
      m_dimension(m_device.get_dimension()),
      m_simulation_options(simulation_option),
      m_boundary_reflection_rng(static_cast<unsigned int>(seed_random_generator + 2)) {
    m_electron_transport.initialize();
    m_hole_transport.initialize();
    validate_time_step_against_scattering_rate();
    initialize_thread_transports(seed_random_generator);
    m_simulation_history.m_initial_seed_rng = seed_random_generator;
    if (number_electrons_start + number_holes_start == 0) {
        initialize_scheduled_particle_injection();
        return;
    }
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
        m_list_particles.push_back(std::make_unique<pbmc_particle>(particle_index, particle_type::electron));
    }

    for (std::size_t i = 0; i < number_holes_start; ++i) {
        const std::size_t particle_index = m_state.m_counter_particles_created++;
        m_list_particles.push_back(std::make_unique<pbmc_particle>(particle_index, particle_type::hole));
    }

    //  Setup the initial position (unique), random number for RPLA and set containing elements for each particles.
    for (auto &p_particle : m_list_particles) {
        auto particle_position = starting_position;
        if (m_dimension == 2) {
            particle_position.set_z(0.0);
        }
        p_particle->set_position(particle_position);
        p_particle->set_weight(1.0);
        p_particle->set_containing_element(first_element);

        initialize_particle_transport_state(*p_particle);
    }

    initialize_scheduled_particle_injection();
}

void device_pbmc_simulation::add_particle_at_position(const mesh::vector3 &location,
                                                      particle_type        type_of_particle,
                                                      double               weight) {
    mesh::vector3 particle_location = location;
    if (m_dimension == 2) {
        particle_location.set_z(0.0);
    }
    mesh::element *first_element{nullptr};
    if (m_dimension == 2) {
        first_element = m_device.find_element_at_location(particle_location);
    } else {
        first_element = m_device.find_element_at_location(particle_location);
    }
    if (first_element == nullptr) {
        std::cout << "Error : particle can't find its first element. No particle created.    " << location << std::endl;
        return;
    }
    const std::size_t idx_particle = m_state.m_counter_particles_created++;
    particle_state    initial_state{};
    initial_state.position          = particle_location;
    initial_state.previous_position = particle_location;
    if (type_of_particle == particle_type::electron) {
        m_list_particles.push_back(
            std::make_unique<pbmc_particle>(idx_particle, particle_type::electron, initial_state, weight));
    } else {
        m_list_particles.push_back(
            std::make_unique<pbmc_particle>(idx_particle, particle_type::hole, initial_state, weight));
    }
    auto &particle = *m_list_particles.back();

    particle.set_containing_element(first_element);
    particle.set_weight(weight);

    initialize_particle_transport_state(particle);
}

void device_pbmc_simulation::add_particles_at_positions(const std::vector<mesh::vector3> &positions,
                                                        particle_type                     type_of_particle,
                                                        double                            weight) {
    // Reserve memory for all particles upfront
    m_list_particles.reserve(m_list_particles.size() + positions.size());
    for (const auto &location : positions) {
        mesh::vector3 particle_location = location;
        if (m_dimension == 2) {
            particle_location.set_z(0.0);
        }
        mesh::element *first_element{nullptr};
        if (m_dimension == 2) {
            first_element = m_device.find_element_at_location(particle_location);
        } else {
            first_element = m_device.find_element_at_location(particle_location);
        }

        if (first_element == nullptr) {
            std::cout << "Error : particle can't find its first element. No particle created.    " << location
                      << std::endl;
            continue;
        }

        const std::size_t idx_particle = m_state.m_counter_particles_created++;
        particle_state    initial_state{};
        initial_state.position          = particle_location;
        initial_state.previous_position = particle_location;
        if (type_of_particle == particle_type::electron) {
            m_list_particles.push_back(
                std::make_unique<pbmc_particle>(idx_particle, particle_type::electron, initial_state, weight));
        } else {
            m_list_particles.push_back(
                std::make_unique<pbmc_particle>(idx_particle, particle_type::hole, initial_state, weight));
        }
        auto &particle = *m_list_particles.back();

        particle.set_containing_element(first_element);
        particle.set_weight(weight);

        initialize_particle_transport_state(particle);
    }
}

std::size_t device_pbmc_simulation::load_particles_from_state_csv(const std::string &filename) {
    std::ifstream stream(filename);
    if (!stream.is_open()) {
        throw std::runtime_error(fmt::format("Could not open initial particle state CSV '{}'.", filename));
    }

    std::string header_line;
    if (!std::getline(stream, header_line)) {
        throw std::invalid_argument(fmt::format("Initial particle state CSV '{}' is empty.", filename));
    }

    const auto indices = csv_header_index(header_line);
    m_list_particles.clear();
    m_state.m_counter_particles_created = 0;

    std::string line;
    std::size_t line_number = 1;
    std::size_t loaded_particles = 0;
    std::size_t next_generated_index = 0;
    std::size_t max_loaded_index_plus_one = 0;

    while (std::getline(stream, line)) {
        ++line_number;
        if (line.empty()) {
            continue;
        }

        const auto fields = split_csv_line(line);

        particle_state state{};
        state.time = 0.0;
        state.position = mesh::vector3{required_csv_double(fields, indices, "x_um", line_number),
                                       required_csv_double(fields, indices, "y_um", line_number),
                                       required_csv_double(fields, indices, "z_um", line_number)};
        state.previous_position = state.position;
        state.local_k = mesh::vector3{required_csv_double(fields, indices, "local_kx_1_per_m", line_number),
                                      required_csv_double(fields, indices, "local_ky_1_per_m", line_number),
                                      required_csv_double(fields, indices, "local_kz_1_per_m", line_number)};
        state.velocity = mesh::vector3{required_csv_double(fields, indices, "vx_m_per_s", line_number),
                                       required_csv_double(fields, indices, "vy_m_per_s", line_number),
                                       required_csv_double(fields, indices, "vz_m_per_s", line_number)};
        state.kinetic_energy = required_csv_double(fields, indices, "energy_eV", line_number);
        state.gamma = optional_csv_double(fields, indices, "gamma_eV", state.kinetic_energy);
        state.valley_index = required_csv_size(fields, indices, "valley_index", line_number);
        if (m_dimension == 2) {
            state.position.set_z(0.0);
            state.previous_position.set_z(0.0);
        }

        const auto type = parse_particle_type_field(required_csv_field(fields, indices, "type", line_number));
        const double weight = required_csv_double(fields, indices, "weight", line_number);
        const std::size_t particle_index =
            optional_csv_size(fields, indices, "particle_index", next_generated_index);
        next_generated_index = std::max(next_generated_index, particle_index + 1);
        max_loaded_index_plus_one = std::max(max_loaded_index_plus_one, particle_index + 1);

        mesh::vector3 lookup_position = state.position;
        if (m_dimension == 2) {
            lookup_position.to_2d_inplace();
        }
        auto *element = m_device.find_element_at_location(lookup_position);
        if (element == nullptr) {
            throw std::invalid_argument(fmt::format(
                "Initial particle state CSV line {} has position outside the device: ({:.6e}, {:.6e}, {:.6e}) um.",
                line_number,
                state.position.x(),
                state.position.y(),
                state.position.z()));
        }

        auto particle = std::make_unique<pbmc_particle>(particle_index, type, state, weight);
        particle->set_containing_element(element);
        particle->set_data_from_device(m_dimension, m_simulation_options.m_enable_impurity_scattering);
        if (m_simulation_options.m_keep_particles_history) {
            particle->record_state();
        }
        m_list_particles.push_back(std::move(particle));
        ++loaded_particles;
    }

    m_state.m_counter_particles_created = std::max(max_loaded_index_plus_one, loaded_particles);
    fmt::print("Loaded {} particles from initial state '{}'.\n", loaded_particles, filename);
    return loaded_particles;
}

std::size_t device_pbmc_simulation::get_number_electrons() const {
    std::size_t nb_electron =
        std::accumulate(m_list_particles.begin(),
                        m_list_particles.end(),
                        0,
                        [](const std::size_t nb_part, const auto &p_part) {
                            return nb_part + static_cast<std::size_t>(p_part->type() == particle_type::electron);
                        });
    return nb_electron;
}

std::size_t device_pbmc_simulation::get_number_holes() const {
    std::size_t nb_hole =
        std::accumulate(m_list_particles.begin(),
                        m_list_particles.end(),
                        0,
                        [](const std::size_t nb_part, const auto &p_part_2) {
                            return nb_part + static_cast<std::size_t>(p_part_2->type() == particle_type::hole);
                        });
    return nb_hole;
}

double device_pbmc_simulation::ramo_current_scale_factor() const { return 1.0; }

// Compute Ramo current for a given single particle
double device_pbmc_simulation::compute_ramo_current_for_particle(const pbmc_particle &particle) const {
    const double scale_factor = ramo_current_scale_factor();
    auto         position     = particle.state().position;
    if (m_dimension == 2) {
        position.to_2d_inplace();
    }
    const auto weighting_field_m    = get_RamoUnitaryElectricField_at_position(position,
                                                                               particle.get_containing_element());
    double     current_contribution = scale_factor * particle.weight() * particle.get_signed_charge() *
                                      particle.state().velocity.dot(weighting_field_m);
    return current_contribution;
}

std::pair<double, double> device_pbmc_simulation::compute_ramo_current() const {
    const auto currents = compute_ramo_currents(true, false);
    return std::make_pair(currents.electron, currents.hole);
}

device_pbmc_simulation::ramo_current_components device_pbmc_simulation::compute_ramo_currents(bool include_full,
                                                                                              bool include_probe) const {
    ramo_current_components currents{};
    const double scale_factor           = ramo_current_scale_factor();
    const bool   current_probe_enabled  = m_simulation_options.m_current_probe.m_enabled;
    const bool   include_probe_currents = include_probe && current_probe_enabled;

    for (const auto &particle : m_list_particles) {
        auto position = particle->state().position;

        if (m_dimension == 2) {
            position.to_2d_inplace();
        }

        const bool inside_probe = current_probe_enabled &&
                                  ((m_dimension == 2)
                                       ? m_simulation_options.m_current_probe.m_box_um.is_inside_2d(position)
                                       : m_simulation_options.m_current_probe.m_box_um.is_inside(position));
        if (current_probe_enabled && !inside_probe) {
            continue;
        }
        if (!include_full && !include_probe_currents) {
            continue;
        }

        const auto weighting_field_m = get_RamoUnitaryElectricField_at_position(position,
                                                                                particle->get_containing_element());
        const auto current = scale_factor * particle->weight() * particle->get_signed_charge() *
                             particle->state().velocity.dot(weighting_field_m);

        if (particle->type() == particle_type::electron) {
            if (include_full) {
                currents.electron += current;
            }
            if (include_probe_currents) {
                currents.probe_electron += current;
            }
        } else {
            if (include_full) {
                currents.hole += current;
            }
            if (include_probe_currents) {
                currents.probe_hole += current;
            }
        }
    }

    return currents;
}

std::pair<double, double> device_pbmc_simulation::compute_probe_ramo_current() const {
    if (!m_simulation_options.m_current_probe.m_enabled) {
        return {0.0, 0.0};
    }
    const auto currents = compute_ramo_currents(false, true);
    return {currents.probe_electron, currents.probe_hole};
}

void device_pbmc_simulation::transport_particles_one_time_step() {
    inject_scheduled_particle_if_due();
    const double                             dt = m_simulation_options.m_time_step;
    std::vector<impact_ionization_pair_seed> impact_pair_seeds;
    const auto                               number_particles = static_cast<std::int64_t>(m_list_particles.size());

#pragma omp parallel for if (m_simulation_options.m_nb_threads > 1) num_threads(m_simulation_options.m_nb_threads)
    for (std::int64_t particle_index = 0; particle_index < number_particles; ++particle_index) {
        auto &particle                     = *m_list_particles[static_cast<std::size_t>(particle_index)];
        particle.state().previous_position = particle.state().position;
        particle.set_data_from_device(m_dimension, m_simulation_options.m_enable_impurity_scattering);
        const auto thread_index = static_cast<std::size_t>(omp_get_thread_num());
        auto      &transport    = transport_for(particle.type(), thread_index);
        // Electric field is in V/cm, but we need it in V/m for the transport kernel, so we convert it here.
        constexpr double cm_to_m = 1.0e2;
        transport.drift_particle(particle, particle.state().electric_field * cm_to_m, dt);
    }
    flatten_particle_positions_for_2d();
    update_element_and_check_boundary();
    remove_collected_particles();

    const auto scattering_particle_count = static_cast<std::int64_t>(m_list_particles.size());
    m_scattering_events_scratch.resize(static_cast<std::size_t>(scattering_particle_count));

#pragma omp parallel for if (m_simulation_options.m_nb_threads > 1) num_threads(m_simulation_options.m_nb_threads)
    for (std::int64_t particle_index = 0; particle_index < scattering_particle_count; ++particle_index) {
        auto &particle = *m_list_particles[static_cast<std::size_t>(particle_index)];
        particle.set_data_from_device(m_dimension, m_simulation_options.m_enable_impurity_scattering);

        // Scattering
        const auto thread_index = static_cast<std::size_t>(omp_get_thread_num());
        auto      &transport    = transport_for(particle.type(), thread_index);
        m_scattering_events_scratch[static_cast<std::size_t>(particle_index)] =
            transport.scatter_particle(particle, dt);

        if (m_simulation_options.m_keep_particles_history) {
            particle.record_state();
        }
    }

    for (std::size_t particle_index = 0; particle_index < m_scattering_events_scratch.size(); ++particle_index) {
        auto &particle = *m_list_particles[particle_index];
        if (m_scattering_events_scratch[particle_index] == scattering_event::impact_ionization) {
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
    }

    for (const auto &seed : impact_pair_seeds) {
        add_particle_at_position(seed.position, particle_type::electron, seed.weight);
        add_particle_at_position(seed.position, particle_type::hole, seed.weight);
    }

    if (!impact_pair_seeds.empty()) {
        flatten_particle_positions_for_2d();
        update_element_and_check_boundary();
        remove_collected_particles();
    }
}

void device_pbmc_simulation::advance_particles_one_time_step() {
    transport_particles_one_time_step();
    m_state.m_time_s += m_simulation_options.m_time_step;
    ++m_state.m_iteration;
}

void device_pbmc_simulation::set_particles_transport_data_from_device() {
    for (auto &p_particle : m_list_particles) {
        p_particle->set_data_from_device(m_dimension, m_simulation_options.m_enable_impurity_scattering);
    }
}

void reflect_particle_to_previous_position(pbmc_particle &particle) {
    auto &state = particle.state();

    state.position = state.previous_position;

    state.velocity *= -1.0;
    state.local_k *= -1.0;
}

void reflect_particle_to_previous_position(pbmc_particle                  &particle,
                                           mesh::boundary_reflection_model reflection_model,
                                           const mesh::element            &old_element,
                                           int                             dimension,
                                           const pbmc_transport_kernel     &transport,
                                           std::minstd_rand               &rng) {
    auto &state = particle.state();

    if (reflection_model == mesh::boundary_reflection_model::reverse) {
        reflect_particle_to_previous_position(particle);
        return;
    }

    const auto hit = mesh::find_boundary_exit_hit(old_element, state.previous_position, state.position, dimension);
    if (!hit.has_value()) {
        reflect_particle_to_previous_position(particle);
        return;
    }

    const mesh::vector3 trial_position         = state.position;
    const mesh::vector3 remaining_displacement = trial_position - hit->position;

    if (reflection_model == mesh::boundary_reflection_model::specular) {
        const mesh::vector3 outgoing_velocity_direction =
            mesh::reflect_vector_specular(state.velocity, hit->inward_normal);
        transport.set_particle_velocity_direction_preserving_energy(particle, outgoing_velocity_direction);
        const mesh::vector3 outgoing_displacement =
            mesh::reflect_vector_specular(remaining_displacement, hit->inward_normal);
        state.position = mesh::place_reflected_position_inside(old_element,
                                                               state.previous_position,
                                                               trial_position,
                                                               *hit,
                                                               outgoing_displacement,
                                                               dimension);
        return;
    }

    const mesh::vector3 outgoing_velocity_direction =
        mesh::draw_diffuse_reflection_vector(state.velocity, hit->inward_normal, dimension, rng);
    transport.set_particle_velocity_direction_preserving_energy(particle, outgoing_velocity_direction);
    const mesh::vector3 outgoing_displacement =
        mesh::align_displacement_with_direction(remaining_displacement, state.velocity);
    state.position = mesh::place_reflected_position_inside(old_element,
                                                           state.previous_position,
                                                           trial_position,
                                                           *hit,
                                                           outgoing_displacement,
                                                           dimension);
}

void device_pbmc_simulation::update_element_and_check_boundary() {
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
            reflect_particle_to_previous_position(particle,
                                                  m_simulation_options.m_boundary_reflection_model,
                                                  *old_element,
                                                  m_dimension,
                                                  transport_for(particle.type()),
                                                  m_boundary_reflection_rng);
            continue;
        }
        if (m_device.get_material_name_at_element(new_element) != "Silicon") {
            reflect_particle_to_previous_position(particle,
                                                  m_simulation_options.m_boundary_reflection_model,
                                                  *old_element,
                                                  m_dimension,
                                                  transport_for(particle.type()),
                                                  m_boundary_reflection_rng);
            continue;
        }
        particle.set_containing_element(new_element);
    }
}

void device_pbmc_simulation::remove_collected_particles() {
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

void device_pbmc_simulation::run() {
    while (m_state.m_time_s < m_simulation_options.m_t_max) {
        if (m_list_particles.empty()) {
            break;
        }
        if (m_simulation_options.m_stop_simu_when_no_electron_remaining && get_number_electrons() == 0) {
            break;
        }
        if (has_reached_particle_limit()) {
            break;
        }
        advance_particles_one_time_step();
        const auto nb_electrons         = get_number_electrons();
        const auto nb_holes             = get_number_holes();
        const auto nb_impact_ionization = m_simulation_history.m_impact_ionization_positions.size();

        double dumb_ramo_current_e_h_total = 0.0;
        double dumb_0                      = 0.0;
        m_simulation_history.add_data_to_history(m_state.m_time_s,
                                                 nb_electrons,
                                                 nb_holes,
                                                 nb_impact_ionization,
                                                 dumb_ramo_current_e_h_total,
                                                 dumb_ramo_current_e_h_total,
                                                 dumb_ramo_current_e_h_total,
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
    export_current_snapshot();
}

std::vector<mesh::vector3> device_pbmc_simulation::get_all_particles_position() const {
    std::vector<mesh::vector3> all_positions(m_list_particles.size());
    std::transform(m_list_particles.begin(), m_list_particles.end(), all_positions.begin(), [](auto &&p_particle) {
        return p_particle->state().position;
    });
    return all_positions;
}

void device_pbmc_simulation::export_current_snapshot() const {
    const std::filesystem::path base_directory(m_simulation_options.m_prefix_export_filename);

    const std::filesystem::path mesh_directory      = base_directory / "mesh";
    const std::filesystem::path particles_directory = base_directory / "particles";

    export_current_mesh_as_vtk(mesh_directory.string());
    export_current_particles_as_vtp(particles_directory.string());

    write_paraview_scene_script(base_directory);
}

std::pair<double, double> device_pbmc_simulation::compute_depletion_region() const {
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

std::string device_pbmc_simulation::initialize_simulation_history_file() {
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

void device_pbmc_simulation::export_current_time_step_as_csv(const std::string &prefix_filename) const {
    const std::string iteration_filename = fmt::format("{}.{:012d}.csv", prefix_filename, m_state.m_iteration);
    export_particle_state_csv(iteration_filename);
}

void device_pbmc_simulation::export_particle_state_csv(const std::string &filename) const {
    std::filesystem::path output_path(filename);
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    std::ofstream stream(filename);

    if (!stream.is_open()) {
        throw std::runtime_error(fmt::format("Could not open particle CSV file '{}'", filename));
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
           << "gamma_eV,"
           << "lattice_temperature_K,"
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
               << state.kinetic_energy << ',' << state.gamma << ',' << state.lattice_temperature_K << ','
               << electric_field.x() << ',' << electric_field.y() << ',' << electric_field.z() << ','
               << electric_field.norm() << ','
               << particle.weight() << ',' << state.valley_index << ',' << particle.get_signed_charge() << '\n';
    }
}
void device_pbmc_simulation::export_all_trajectories_as_csv(const std::string &prefix_filename) const {
    for (auto &&p_particle : m_list_particles) {
        // std::cout << "\rExporting trajectory of particle " << p_particle->get_index() << std::flush;
        std::string filename = fmt::format("{}particle_{:06}_trajectory.csv", prefix_filename, p_particle->index());
        p_particle->export_trajectory_as_csv(filename);
    }
    // std::cout << std::endl;
}

void device_pbmc_simulation::export_current_particles_as_vtp(const std::string &directory) const {
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

    stream << "        <DataArray type=\"Float64\" Name=\"lattice_temperature_K\" format=\"ascii\">\n";
    stream << "          ";
    for (const auto &p_particle : m_list_particles) {
        stream << p_particle->state().lattice_temperature_K << ' ';
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

void device_pbmc_simulation::write_particle_vtp_time_collection(const std::string &pvd_filename) const {
    write_vtk_time_collection(pvd_filename, m_particle_vtp_export_records);
}

double device_pbmc_simulation::current_density_cell_volume_m3(const mesh::element &element) const {
    const double measure_um = std::abs(element.get_measure());
    return measure_um * std::pow(uepm::units::micron_to_meter, m_dimension);
}

void device_pbmc_simulation::publish_mesh_particle_local_average_energy() const {
    auto *mesh = m_device.get_p_mesh();
    if (mesh == nullptr) {
        return;
    }

    const auto                                             list_bulk_elements = mesh->get_list_bulk_element();
    std::unordered_map<const mesh::element *, std::size_t> element_indices;
    element_indices.reserve(list_bulk_elements.size());
    for (std::size_t index = 0; index < list_bulk_elements.size(); ++index) {
        element_indices.emplace(list_bulk_elements[index].get(), index);
    }

    std::vector<double> particle_count(list_bulk_elements.size(), 0.0);
    std::vector<double> energy_sum_eV(list_bulk_elements.size(), 0.0);

    for (const auto &p_particle : m_list_particles) {
        const auto index_it = element_indices.find(p_particle->get_containing_element());
        if (index_it == element_indices.end()) {
            continue;
        }
        const std::size_t index = index_it->second;
        particle_count[index] += 1.0;
        energy_sum_eV[index] += p_particle->state().kinetic_energy;
    }

    std::vector<double> average_energy_eV(list_bulk_elements.size(), 0.0);
    for (std::size_t index = 0; index < list_bulk_elements.size(); ++index) {
        if (particle_count[index] > 0.0) {
            average_energy_eV[index] = energy_sum_eV[index] / particle_count[index];
        }
    }

    const std::string field_name = "particle_local_average_energy_eV";
    if (mesh->scalar_function_exists(field_name)) {
        mesh->remove_scalar_function(field_name);
    }

    for (const auto &region : mesh->get_list_bulk_region()) {
        const auto               region_elements = region.get_list_elements();
        std::vector<double>      region_values;
        std::vector<std::size_t> region_element_indices;
        region_values.reserve(region_elements.size());
        region_element_indices.reserve(region_elements.size());

        for (const auto &element : region_elements) {
            const auto index_it = element_indices.find(element.get());
            region_values.push_back(index_it == element_indices.end() ? 0.0 : average_energy_eV[index_it->second]);
            region_element_indices.push_back(element->get_index());
        }

        auto dataset = std::make_shared<mesh::dataset<double>>(field_name,
                                                               mesh->get_total_number_dataset() + 1,
                                                               region.get_index(),
                                                               region_values,
                                                               region_element_indices,
                                                               mesh::DataType::scalar,
                                                               mesh::DataLocationType::cell,
                                                               1);
        mesh->add_scalar_dataset(dataset);
        mesh->add_scalar_data_to_elements(*dataset);
    }
}

void device_pbmc_simulation::publish_mesh_particle_local_current_density() const {
    auto *mesh = m_device.get_p_mesh();
    if (mesh == nullptr) {
        return;
    }

    const auto                                             list_bulk_elements = mesh->get_list_bulk_element();
    std::unordered_map<const mesh::element *, std::size_t> element_indices;
    element_indices.reserve(list_bulk_elements.size());
    for (std::size_t index = 0; index < list_bulk_elements.size(); ++index) {
        element_indices.emplace(list_bulk_elements[index].get(), index);
    }

    std::vector<mesh::vector3> current_density_A_per_m2(list_bulk_elements.size(), mesh::vector3{0.0, 0.0, 0.0});
    for (const auto &p_particle : m_list_particles) {
        const auto index_it = element_indices.find(p_particle->get_containing_element());
        if (index_it == element_indices.end()) {
            continue;
        }

        const std::size_t index = index_it->second;
        current_density_A_per_m2[index] +=
            p_particle->weight() * p_particle->get_signed_charge() * p_particle->state().velocity;
    }

    for (std::size_t index = 0; index < list_bulk_elements.size(); ++index) {
        const double cell_volume_m3 = current_density_cell_volume_m3(*list_bulk_elements[index]);
        if (cell_volume_m3 > 0.0) {
            current_density_A_per_m2[index] /= cell_volume_m3;
        }
    }

    const std::string field_name = "particle_local_current_density_A_per_m2";
    if (mesh->vector_function_exists(field_name)) {
        mesh->remove_vector_function(field_name);
    }

    for (const auto &region : mesh->get_list_bulk_region()) {
        const auto                 region_elements = region.get_list_elements();
        std::vector<mesh::vector3> region_values;
        std::vector<std::size_t>   region_element_indices;
        region_values.reserve(region_elements.size());
        region_element_indices.reserve(region_elements.size());

        for (const auto &element : region_elements) {
            const auto index_it = element_indices.find(element.get());
            region_values.push_back(index_it == element_indices.end() ? mesh::vector3{0.0, 0.0, 0.0}
                                                                      : current_density_A_per_m2[index_it->second]);
            region_element_indices.push_back(element->get_index());
        }

        auto dataset = std::make_shared<mesh::dataset<mesh::vector3>>(field_name,
                                                                      mesh->get_total_number_dataset() + 1,
                                                                      region.get_index(),
                                                                      region_values,
                                                                      region_element_indices,
                                                                      mesh::DataType::vector,
                                                                      mesh::DataLocationType::cell,
                                                                      m_dimension);
        mesh->add_vector_dataset(dataset);
        mesh->add_vector_data_to_elements(*dataset);
    }
}

void device_pbmc_simulation::export_current_mesh_as_vtk(const std::string &directory) const {
    const std::filesystem::path output_directory(directory);
    std::filesystem::create_directories(output_directory);

    const std::string filename = fmt::format("mesh_{:012d}.vtu", m_state.m_iteration);

    const std::filesystem::path vtu_path = output_directory / filename;
    const std::filesystem::path pvd_path = output_directory / "mesh.pvd";

    publish_mesh_particle_local_average_energy();
    publish_mesh_particle_local_current_density();
    file::export_as_vtu(*(m_device.get_p_mesh()), vtu_path.string());

    const auto already_recorded =
        std::find_if(m_mesh_vtk_export_records.begin(),
                     m_mesh_vtk_export_records.end(),
                     [&](const vtk_time_series_record &record) { return record.m_filename == filename; });

    if (already_recorded == m_mesh_vtk_export_records.end()) {
        m_mesh_vtk_export_records.push_back(
            vtk_time_series_record{.m_time_s = m_state.m_time_s, .m_filename = filename});
    }

    write_mesh_vtk_time_collection(pvd_path.string());

    const bool export_x_cut_enabled = true;
    if (export_x_cut_enabled) {
        const auto                  device_box = m_device.get_p_mesh()->get_bounding_box();
        const auto                  y_middle   = 0.5 * (device_box.get_y_min() + device_box.get_y_max());
        const auto                  z_middle   = 0.5 * (device_box.get_z_min() + device_box.get_z_max());
        const auto                  dx         = 1e-3;  // 1 nm
        const std::filesystem::path x_cut_path =
            output_directory / fmt::format("mesh_x_cut_{:012d}.csv", m_state.m_iteration);
        m_device.get_p_mesh()->export_x_cut(x_cut_path.string(), y_middle, z_middle, dx);
    }
}

void device_pbmc_simulation::write_mesh_vtk_time_collection(const std::string &pvd_filename) const {
    write_vtk_time_collection(pvd_filename, m_mesh_vtk_export_records);
}

}  // namespace uepm::PBMC
