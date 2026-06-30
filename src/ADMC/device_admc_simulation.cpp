/**
 * @file device_admc_simulation.cpp
 * @brief Device-level ADMC drift-diffusion particle simulation.
 */

#include "device_admc_simulation.hpp"

#include <fmt/core.h>
#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <stdexcept>
#include <unordered_map>

#include "physical_constants.hpp"
#include "unit_conversion.hpp"
#include "vtkWriter.hpp"

namespace uepm::ADMC {
namespace {

double finite_or_zero(double value) { return std::isfinite(value) ? value : 0.0; }

double local_net_doping_cm_3(const mesh::element& element, const mesh::vector3& position_um) {
    return finite_or_zero(element.interpolate_doping_at_location(position_um));
}

double signed_particle_charge_C(const device_admc_particle& particle) {
    return particle.weight * carrier_charge_sign(particle.particle.type()) * uepm::constants::q_e;
}

double average_or_zero(double sum, double count) { return count > 0.0 ? sum / count : 0.0; }

void write_vtk_time_collection(const std::string& pvd_filename, std::vector<admc_vtk_time_series_record> records) {
    std::ofstream stream(pvd_filename);
    if (!stream.is_open()) {
        throw std::runtime_error("Could not open ADMC VTK collection file '" + pvd_filename + "'.");
    }

    std::sort(records.begin(), records.end(), [](const auto& lhs, const auto& rhs) { return lhs.time_s < rhs.time_s; });

    stream << std::setprecision(std::numeric_limits<double>::max_digits10);
    stream << "<?xml version=\"1.0\"?>\n";
    stream << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    stream << "  <Collection>\n";
    for (const auto& record : records) {
        stream << "    <DataSet timestep=\"" << record.time_s << "\" group=\"\" part=\"0\" file=\"" << record.filename
               << "\"/>\n";
    }
    stream << "  </Collection>\n";
    stream << "</VTKFile>\n";
}

std::string python_string_literal(const std::string& value) {
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

void write_paraview_scene_script(const std::filesystem::path& base_directory) {
    std::filesystem::create_directories(base_directory);
    const auto    script_path = base_directory / "open_scene.py";
    const auto    scene_dir   = std::filesystem::absolute(base_directory).lexically_normal();
    std::ofstream stream(script_path);
    if (!stream.is_open()) {
        throw std::runtime_error("Could not open ADMC ParaView scene script '" + script_path.string() + "'.");
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

void options_device_ADMC::validate() const {
    if (!std::isfinite(m_lattice_temperature_K) || m_lattice_temperature_K <= 0.0) {
        throw std::invalid_argument("ADMC device temperature must be positive and finite.");
    }
    if (!std::isfinite(m_time_step_s) || m_time_step_s <= 0.0) {
        throw std::invalid_argument("ADMC device time step must be positive and finite.");
    }
    if (!std::isfinite(m_final_time_s) || m_final_time_s < 0.0) {
        throw std::invalid_argument("ADMC device final time must be finite and non-negative.");
    }
    if (m_max_number_particles == 0) {
        throw std::invalid_argument("ADMC device max particle count must be positive.");
    }
    if (m_frequency_export <= 0) {
        throw std::invalid_argument("ADMC export frequency must be positive.");
    }
    if (m_enable_scheduled_particle_injection) {
        if (!std::isfinite(m_scheduled_injection_time_s) || m_scheduled_injection_time_s < 0.0) {
            throw std::invalid_argument("ADMC scheduled injection time must be non-negative and finite.");
        }
        if (m_scheduled_injection_time_s > m_final_time_s) {
            throw std::invalid_argument("ADMC scheduled injection time cannot exceed final time.");
        }
        if (!std::isfinite(m_scheduled_injection_weight) || m_scheduled_injection_weight <= 0.0) {
            throw std::invalid_argument("ADMC scheduled injection weight must be positive and finite.");
        }
    }
}

void history_device_ADMC::add(double      time_s,
                              std::size_t electrons,
                              std::size_t holes,
                              double      electron_current_A,
                              double      hole_current_A,
                              double      total_current_A,
                              double      probe_electron_current_A,
                              double      probe_hole_current_A,
                              double      probe_total_current_A,
                              double      max_field_V_per_m) {
    times_s.push_back(time_s);
    nb_electrons.push_back(electrons);
    nb_holes.push_back(holes);
    ramo_current_electron_A.push_back(electron_current_A);
    ramo_current_hole_A.push_back(hole_current_A);
    ramo_current_A.push_back(total_current_A);
    probe_ramo_current_electron_A.push_back(probe_electron_current_A);
    probe_ramo_current_hole_A.push_back(probe_hole_current_A);
    probe_ramo_current_A.push_back(probe_total_current_A);
    max_electric_field_V_per_m.push_back(max_field_V_per_m);
}

void history_device_ADMC::print_header_csv(const std::string& filename) const {
    std::ofstream stream(filename);
    if (!stream.is_open()) {
        throw std::runtime_error("Could not open ADMC history CSV file '" + filename + "'.");
    }
    stream << "time,nb_electrons,nb_holes,nb_impact_ionization,ramo_current_electron,ramo_current_hole,ramo_current,"
              "probe_ramo_current_electron,probe_ramo_current_hole,probe_ramo_current,"
              "max_electric_field,ramo_electrode_voltage_V,reference_electrode_voltage_V,quench_bias_voltage_V,"
              "quench_device_current_A,"
              "quench_resistor_current_A,quench_voltage_drop_V\n";
}

void history_device_ADMC::append_last_iter_to_csv(std::fstream& file) const {
    if (times_s.empty()) {
        return;
    }
    const std::size_t i = times_s.size() - 1;
    constexpr double electric_field_V_per_m_to_V_per_cm = 0.01;
    file << times_s[i] << ',' << nb_electrons[i] << ',' << nb_holes[i] << ',' << 0 << ','
         << ramo_current_electron_A[i] << ',' << ramo_current_hole_A[i] << ',' << ramo_current_A[i] << ','
         << probe_ramo_current_electron_A[i] << ',' << probe_ramo_current_hole_A[i] << ','
         << probe_ramo_current_A[i] << ',' << max_electric_field_V_per_m[i] * electric_field_V_per_m_to_V_per_cm
         << ',' << 0.0 << ',' << 0.0 << ',' << 0.0 << ',' << 0.0 << ',' << 0.0 << ',' << 0.0 << '\n';
}

void history_device_ADMC::export_to_csv(const std::string& filename) const {
    std::ofstream stream(filename);
    if (!stream.is_open()) {
        throw std::runtime_error("Could not open ADMC history CSV file '" + filename + "'.");
    }
    stream.close();
    print_header_csv(filename);
    stream.open(filename, std::ios::app);
    if (!stream.is_open()) {
        throw std::runtime_error("Could not open ADMC history CSV file '" + filename + "'.");
    }
    stream << std::setprecision(std::numeric_limits<double>::max_digits10);
    constexpr double electric_field_V_per_m_to_V_per_cm = 0.01;
    for (std::size_t i = 0; i < times_s.size(); ++i) {
        stream << times_s[i] << ',' << nb_electrons[i] << ',' << nb_holes[i] << ',' << 0 << ','
               << ramo_current_electron_A[i] << ',' << ramo_current_hole_A[i] << ',' << ramo_current_A[i] << ','
               << probe_ramo_current_electron_A[i] << ',' << probe_ramo_current_hole_A[i] << ','
               << probe_ramo_current_A[i] << ',' << max_electric_field_V_per_m[i] * electric_field_V_per_m_to_V_per_cm
               << ',' << 0.0 << ',' << 0.0 << ',' << 0.0 << ',' << 0.0 << ',' << 0.0 << ',' << 0.0 << '\n';
    }
}

device_admc_simulation::device_admc_simulation(const device::device&      simulation_device,
                                               const options_device_ADMC& options,
                                               std::uint64_t              random_seed)
    : m_device(simulation_device),
      m_options(options),
      m_dimension(m_device.get_dimension()),
      m_random_generator(random_seed) {
    m_options.validate();
    initialize_scheduled_particle_injection();
}

device_admc_simulation::device_admc_simulation(const device::device&      simulation_device,
                                               const options_device_ADMC& options,
                                               const mesh::vector3&       starting_position_um,
                                               std::size_t                number_electrons_start,
                                               std::size_t                number_holes_start,
                                               std::uint64_t              random_seed)
    : device_admc_simulation(simulation_device, options, random_seed) {
    if (number_holes_start > std::numeric_limits<std::size_t>::max() - number_electrons_start) {
        throw std::invalid_argument("ADMC device particle count overflows size_t.");
    }
    m_particles.reserve(number_electrons_start + number_holes_start);
    for (std::size_t i = 0; i < number_electrons_start; ++i) {
        add_particle_at_position(starting_position_um, carrier_type::electron);
    }
    for (std::size_t i = 0; i < number_holes_start; ++i) {
        add_particle_at_position(starting_position_um, carrier_type::hole);
    }
}

mesh::vector3 device_admc_simulation::normalize_mesh_position_for_dimension(mesh::vector3 position_um) const {
    if (m_dimension == 2) {
        position_um.to_2d_inplace();
    }
    return position_um;
}

mesh::vector3 device_admc_simulation::to_mesh_position_um(const vector3& position_m) const {
    return normalize_mesh_position_for_dimension(position_m * uepm::units::meter_to_micron);
}

vector3 device_admc_simulation::to_admc_position_m(const mesh::vector3& position_um) const {
    return normalize_mesh_position_for_dimension(position_um) * uepm::units::micron_to_meter;
}

bool device_admc_simulation::is_transport_material_element(mesh::element& element) {
    const std::string material_name = m_device.get_material_name_at_element(&element);
    return material_name == "Si" || material_name == "Silicon";
}

vector3 device_admc_simulation::draw_standard_normal() {
    return {m_standard_normal(m_random_generator),
            m_standard_normal(m_random_generator),
            m_standard_normal(m_random_generator)};
}

vector3 device_admc_simulation::get_RamoUnitaryElectricField_at_position(const mesh::vector3& position) const {
    if (m_state.m_use_constant_RamoUnitaryElectricField) {
        return m_state.m_RamoUnitaryElectricField_Vm_per_cm;
    }
    return m_device.interpolate_vector_at_location("RamoUnitaryPotential_gradient", position);
}

void device_admc_simulation::initialize_scheduled_particle_injection() {
    m_state.m_scheduled_particle_injection_done = !m_options.m_enable_scheduled_particle_injection;
}

bool device_admc_simulation::has_pending_scheduled_particle_injection() const {
    return m_options.m_enable_scheduled_particle_injection && !m_state.m_scheduled_particle_injection_done;
}

void device_admc_simulation::inject_scheduled_particle_if_due() {
    if (!has_pending_scheduled_particle_injection()) {
        return;
    }
    if (m_state.m_time_s + m_options.m_time_step_s < m_options.m_scheduled_injection_time_s) {
        return;
    }
    add_particle_at_position(m_options.m_scheduled_injection_position_um,
                             m_options.m_scheduled_injection_type,
                             m_options.m_scheduled_injection_weight);
    m_state.m_scheduled_particle_injection_done = true;
    fmt::print("ADMC scheduled particle injected at t={:.6e} s\n", m_options.m_scheduled_injection_time_s);
}

void device_admc_simulation::add_particle_at_position(const mesh::vector3& position_um,
                                                      carrier_type         type,
                                                      double               weight) {
    if (m_particles.size() >= m_options.m_max_number_particles) {
        return;
    }
    if (!std::isfinite(weight) || weight <= 0.0) {
        throw std::invalid_argument("ADMC device particle weight must be positive and finite.");
    }

    const mesh::vector3 mesh_position_um = normalize_mesh_position_for_dimension(position_um);
    auto*               element          = m_device.find_element_at_location(mesh_position_um);
    if (element == nullptr) {
        return;
    }
    if (!is_transport_material_element(*element)) {
        return;
    }

    const std::size_t index = m_state.m_counter_particles_created++;
    m_particles.push_back(device_admc_particle{
        .particle           = admc_particle(index, type, to_admc_position_m(mesh_position_um)),
        .containing_element = element,
        .weight             = weight,
        .crossed_contact    = false,
    });
    initialize_particle_device_state(m_particles.back());
}

admc_local_environment device_admc_simulation::local_environment(const device_admc_particle& particle) const {
    if (particle.containing_element == nullptr) {
        throw std::runtime_error("ADMC device particle has no containing element.");
    }

    const mesh::vector3 position_um = to_mesh_position_um(particle.particle.state().position_m);
    const mesh::vector3 electric_field_V_per_cm =
        particle.containing_element->interpolate_electric_field_at_location(position_um);
    const double net_doping_cm_3 = local_net_doping_cm_3(*particle.containing_element, position_um);

    return {
        .electric_field_V_per_m    = electric_field_V_per_cm * uepm::units::electric_field_V_per_cm_to_V_per_m,
        .doping_concentration_cm_3 = std::abs(net_doping_cm_3),
        .lattice_temperature_K     = m_options.m_lattice_temperature_K,
    };
}

void device_admc_simulation::initialize_particle_device_state(device_admc_particle& particle) {
    const auto environment          = local_environment(particle);
    auto&      state                = particle.particle.state();
    state.electric_field_V_per_m    = environment.electric_field_V_per_m;
    state.doping_concentration_cm_3 = environment.doping_concentration_cm_3;
    state.lattice_temperature_K     = environment.lattice_temperature_K;
}

void device_admc_simulation::update_element_and_check_boundary(device_admc_particle& particle) {
    if (particle.containing_element == nullptr) {
        particle.crossed_contact = true;
        return;
    }

    const mesh::vector3 current_position_um = to_mesh_position_um(particle.particle.state().position_m);
    if (particle.containing_element->is_location_inside_element(current_position_um)) {
        return;
    }
    if (m_device.check_enters_contact(current_position_um)) {
        particle.crossed_contact = true;
        return;
    }

    auto* new_element = m_device.find_element_at_location(current_position_um);
    if (new_element == nullptr || !is_transport_material_element(*new_element)) {
        auto&              state                = particle.particle.state();
        const mesh::vector3 previous_position_um = to_mesh_position_um(state.previous_position_m);
        const auto          hit = mesh::find_boundary_exit_hit(*particle.containing_element,
                                                              previous_position_um,
                                                              current_position_um,
                                                              m_dimension);
        if (m_options.m_boundary_reflection_model == mesh::boundary_reflection_model::reverse || !hit.has_value()) {
            state.position_m = state.previous_position_m;
            state.total_velocity_m_per_s *= -1.0;
            state.drift_velocity_m_per_s *= -1.0;
            return;
        }

        const mesh::vector3 remaining_displacement_um = current_position_um - hit->position;
        mesh::vector3       outgoing_displacement_um;

        if (m_options.m_boundary_reflection_model == mesh::boundary_reflection_model::specular) {
            state.total_velocity_m_per_s =
                mesh::reflect_vector_specular(state.total_velocity_m_per_s, hit->inward_normal);
            state.drift_velocity_m_per_s =
                mesh::reflect_vector_specular(state.drift_velocity_m_per_s, hit->inward_normal);
            outgoing_displacement_um =
                mesh::reflect_vector_specular(remaining_displacement_um, hit->inward_normal);
        } else {
            state.total_velocity_m_per_s = mesh::draw_diffuse_reflection_vector(state.total_velocity_m_per_s,
                                                                                hit->inward_normal,
                                                                                m_dimension,
                                                                                m_random_generator);
            state.drift_velocity_m_per_s = mesh::draw_diffuse_reflection_vector(state.drift_velocity_m_per_s,
                                                                                hit->inward_normal,
                                                                                m_dimension,
                                                                                m_random_generator);
            outgoing_displacement_um =
                mesh::align_displacement_with_direction(remaining_displacement_um, state.total_velocity_m_per_s);
        }
        state.position_m = to_admc_position_m(mesh::place_reflected_position_inside(*particle.containing_element,
                                                                                    previous_position_um,
                                                                                    current_position_um,
                                                                                    *hit,
                                                                                    outgoing_displacement_um,
                                                                                    m_dimension));
        return;
    }
    particle.containing_element = new_element;
}

void device_admc_simulation::remove_collected_particles() {
    std::erase_if(m_particles, [](const auto& particle) { return particle.crossed_contact; });
}

void device_admc_simulation::advance_particles_one_time_step() {
    inject_scheduled_particle_if_due();
    const double dt_s = std::min(m_options.m_time_step_s, m_options.m_final_time_s - m_state.m_time_s);
    if (dt_s <= 0.0) {
        return;
    }

    for (auto& particle : m_particles) {
        const auto environment = local_environment(particle);
        m_transport.step(particle.particle, environment, dt_s, draw_standard_normal());
        if (m_dimension == 2) {
            particle.particle.state().position_m.set_z(0.0);
        }
        update_element_and_check_boundary(particle);
    }
    remove_collected_particles();
    m_state.m_time_s += dt_s;
    ++m_state.m_iteration;
    const auto [electron_current_A, hole_current_A] = compute_ramo_current();
    m_state.m_last_ramo_current_electron_A          = electron_current_A;
    m_state.m_last_ramo_current_hole_A              = hole_current_A;
    if (m_options.m_current_probe.m_enabled) {
        m_state.m_last_probe_ramo_current_electron_A = electron_current_A;
        m_state.m_last_probe_ramo_current_hole_A     = hole_current_A;
    } else {
        m_state.m_last_probe_ramo_current_electron_A = 0.0;
        m_state.m_last_probe_ramo_current_hole_A     = 0.0;
    }
    record_history(electron_current_A, hole_current_A);
}

void device_admc_simulation::run() {
    while (m_state.m_time_s < m_options.m_final_time_s) {
        if (m_particles.empty() && !has_pending_scheduled_particle_injection()) {
            break;
        }
        if (m_options.m_stop_when_no_electrons && get_number_electrons() == 0 &&
            !has_pending_scheduled_particle_injection()) {
            break;
        }
        advance_particles_one_time_step();
    }
    export_current_snapshot();
}

void device_admc_simulation::record_history(double electron_current_A, double hole_current_A) {
    m_history.add(m_state.m_time_s,
                  get_number_electrons(),
                  get_number_holes(),
                  electron_current_A,
                  hole_current_A,
                  electron_current_A + hole_current_A,
                  m_state.m_last_probe_ramo_current_electron_A,
                  m_state.m_last_probe_ramo_current_hole_A,
                  m_state.m_last_probe_ramo_current_electron_A + m_state.m_last_probe_ramo_current_hole_A,
                  max_particle_electric_field_V_per_m());
}

std::string device_admc_simulation::initialize_simulation_history_file() {
    const std::string simulation_name_for_file =
        m_options.m_simulation_name.empty() ? "simulation" : m_options.m_simulation_name;
    const std::string history_filename =
        fmt::format("{}/{}_history.csv", m_options.m_output_directory, simulation_name_for_file);
    fmt::print("Initializing ADMC simulation history file at '{}'\n", history_filename);
    std::filesystem::create_directories(m_options.m_output_directory);
    m_history.print_header_csv(history_filename);
    return history_filename;
}

std::pair<double, double> device_admc_simulation::last_ramo_current() const {
    return {m_state.m_last_ramo_current_electron_A, m_state.m_last_ramo_current_hole_A};
}

std::pair<double, double> device_admc_simulation::last_probe_ramo_current() const {
    return {m_state.m_last_probe_ramo_current_electron_A, m_state.m_last_probe_ramo_current_hole_A};
}

std::size_t device_admc_simulation::get_number_electrons() const {
    return static_cast<std::size_t>(std::count_if(m_particles.begin(), m_particles.end(), [](const auto& particle) {
        return particle.particle.type() == carrier_type::electron;
    }));
}

std::size_t device_admc_simulation::get_number_holes() const {
    return static_cast<std::size_t>(std::count_if(m_particles.begin(), m_particles.end(), [](const auto& particle) {
        return particle.particle.type() == carrier_type::hole;
    }));
}

double device_admc_simulation::get_total_electron_weight() const {
    double total = 0.0;
    for (const auto& particle : m_particles) {
        if (particle.particle.type() == carrier_type::electron) {
            total += particle.weight;
        }
    }
    return total;
}

double device_admc_simulation::get_total_hole_weight() const {
    double total = 0.0;
    for (const auto& particle : m_particles) {
        if (particle.particle.type() == carrier_type::hole) {
            total += particle.weight;
        }
    }
    return total;
}

std::pair<double, double> device_admc_simulation::compute_ramo_current() const {
    double electron_current_A = 0.0;
    double hole_current_A     = 0.0;

    for (const auto& particle : m_particles) {
        mesh::vector3 position_um = to_mesh_position_um(particle.particle.state().position_m);
        if (m_dimension == 2) {
            position_um.to_2d_inplace();
        }
        if (m_options.m_current_probe.m_enabled) {
            const bool inside_probe = (m_dimension == 2)
                                          ? m_options.m_current_probe.m_box_um.is_inside_2d(position_um)
                                          : m_options.m_current_probe.m_box_um.is_inside(position_um);
            if (!inside_probe) {
                continue;
            }
        }
        const vector3       weighting_field = get_RamoUnitaryElectricField_at_position(position_um);
        const double        current_A       = particle.weight * carrier_charge_sign(particle.particle.type()) *
                                 uepm::constants::q_e *
                                 particle.particle.state().total_velocity_m_per_s.dot(weighting_field);
        if (particle.particle.type() == carrier_type::electron) {
            electron_current_A += current_A;
        } else {
            hole_current_A += current_A;
        }
    }
    return {electron_current_A, hole_current_A};
}

std::pair<double, double> device_admc_simulation::compute_probe_ramo_current() const {
    if (!m_options.m_current_probe.m_enabled) {
        return {0.0, 0.0};
    }
    return compute_ramo_current();
}

double device_admc_simulation::max_particle_electric_field_V_per_m() const {
    double max_field = 0.0;
    for (const auto& particle : m_particles) {
        max_field = std::max(max_field, particle.particle.state().electric_field_V_per_m.norm());
    }
    return max_field;
}

void device_admc_simulation::export_current_time_step_as_csv(const std::string& prefix_filename) const {
    const std::filesystem::path prefix(prefix_filename);
    const std::filesystem::path parent_directory =
        prefix.parent_path().empty() ? std::filesystem::path(".") : prefix.parent_path();
    std::filesystem::create_directories(parent_directory);
    const std::string filename = fmt::format("{}.{:012d}.csv", prefix_filename, m_state.m_iteration);
    std::ofstream     stream(filename);
    if (!stream.is_open()) {
        throw std::runtime_error("Could not open ADMC particle CSV file '" + filename + "'.");
    }
    stream << std::setprecision(std::numeric_limits<double>::max_digits10);
    stream << "particle_index,type,time_s,x_um,y_um,z_um,vx_m_per_s,vy_m_per_s,vz_m_per_s,"
              "drift_vx_m_per_s,drift_vy_m_per_s,drift_vz_m_per_s,mobility_m2_per_V_s,diffusion_m2_per_s,"
              "lattice_temperature_K,doping_concentration_cm_3,electric_field_x_V_per_m,"
              "electric_field_y_V_per_m,electric_field_z_V_per_m,electric_field_norm_V_per_m,"
              "electric_field_x_V_per_cm,electric_field_y_V_per_cm,electric_field_z_V_per_cm,"
              "electric_field_norm_V_per_cm,weight,signed_charge_C\n";
    for (const auto& device_particle : m_particles) {
        const auto& particle                = device_particle.particle;
        const auto& state                   = particle.state();
        const auto  position_um             = to_mesh_position_um(state.position_m);
        const auto  electric_field_V_per_cm = 0.01 * state.electric_field_V_per_m;
        stream << particle.index() << ',' << static_cast<int>(particle.type()) << ',' << state.time_s << ','
               << position_um.x() << ',' << position_um.y() << ',' << position_um.z() << ','
               << state.total_velocity_m_per_s.x() << ',' << state.total_velocity_m_per_s.y() << ','
               << state.total_velocity_m_per_s.z() << ',' << state.drift_velocity_m_per_s.x() << ','
               << state.drift_velocity_m_per_s.y() << ',' << state.drift_velocity_m_per_s.z() << ','
               << state.mobility_m2_per_V_s << ',' << state.diffusion_m2_per_s << ',' << state.lattice_temperature_K
               << ',' << state.doping_concentration_cm_3 << ',' << state.electric_field_V_per_m.x() << ','
               << state.electric_field_V_per_m.y() << ',' << state.electric_field_V_per_m.z() << ','
               << state.electric_field_V_per_m.norm() << ',' << electric_field_V_per_cm.x() << ','
               << electric_field_V_per_cm.y() << ',' << electric_field_V_per_cm.z() << ','
               << electric_field_V_per_cm.norm() << ',' << device_particle.weight << ','
               << signed_particle_charge_C(device_particle) << '\n';
    }
}

void device_admc_simulation::export_current_particles_as_vtp(const std::string& directory) const {
    const std::filesystem::path output_directory(directory);
    std::filesystem::create_directories(output_directory);
    const std::filesystem::path vtp_path = output_directory / fmt::format("particles_{:012d}.vtp", m_state.m_iteration);
    const std::string           filename = vtp_path.filename().string();
    const std::filesystem::path pvd_path = output_directory / "particles.pvd";
    std::ofstream               stream(vtp_path);
    if (!stream.is_open()) {
        throw std::runtime_error("Could not open ADMC particle VTP file '" + vtp_path.string() + "'.");
    }
    const std::size_t n = m_particles.size();
    stream << std::setprecision(std::numeric_limits<double>::max_digits10);
    stream << "<?xml version=\"1.0\"?>\n";
    stream << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    stream << "  <PolyData>\n";
    stream << "    <Piece NumberOfPoints=\"" << n << "\" NumberOfVerts=\"" << n
           << "\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
    stream << "      <PointData Scalars=\"mobility_m2_per_V_s\" Vectors=\"velocity_m_per_s\">\n";
    stream << "        <DataArray type=\"Int64\" Name=\"particle_index\" format=\"ascii\">\n          ";
    for (const auto& p : m_particles) {
        stream << static_cast<long long>(p.particle.index()) << ' ';
    }
    stream << "\n        </DataArray>\n";
    stream << "        <DataArray type=\"Int32\" Name=\"particle_type\" format=\"ascii\">\n          ";
    for (const auto& p : m_particles) {
        stream << static_cast<int>(p.particle.type()) << ' ';
    }
    stream << "\n        </DataArray>\n";
    stream << "        <DataArray type=\"Float64\" Name=\"time_s\" format=\"ascii\">\n          ";
    for (const auto& p : m_particles) {
        stream << p.particle.state().time_s << ' ';
    }
    stream << "\n        </DataArray>\n";
    stream << "        <DataArray type=\"Float64\" Name=\"mobility_m2_per_V_s\" format=\"ascii\">\n          ";
    for (const auto& p : m_particles) {
        stream << p.particle.state().mobility_m2_per_V_s << ' ';
    }
    stream << "\n        </DataArray>\n";
    stream << "        <DataArray type=\"Float64\" Name=\"weight\" format=\"ascii\">\n          ";
    for (const auto& p : m_particles) {
        stream << p.weight << ' ';
    }
    stream << "\n        </DataArray>\n";
    stream << "        <DataArray type=\"Float64\" Name=\"signed_charge_C\" format=\"ascii\">\n          ";
    for (const auto& p : m_particles) {
        stream << signed_particle_charge_C(p) << ' ';
    }
    stream << "\n        </DataArray>\n";
    stream << "        <DataArray type=\"Float64\" Name=\"electric_field_norm_V_per_m\" format=\"ascii\">\n          ";
    for (const auto& p : m_particles) {
        stream << p.particle.state().electric_field_V_per_m.norm() << ' ';
    }
    stream << "\n        </DataArray>\n";
    stream << "        <DataArray type=\"Float64\" Name=\"electric_field_norm_V_per_cm\" format=\"ascii\">\n          ";
    for (const auto& p : m_particles) {
        stream << 0.01 * p.particle.state().electric_field_V_per_m.norm() << ' ';
    }
    stream << "\n        </DataArray>\n";
    stream << "        <DataArray type=\"Float64\" Name=\"velocity_m_per_s\" NumberOfComponents=\"3\" "
              "format=\"ascii\">\n          ";
    for (const auto& p : m_particles) {
        const auto& v = p.particle.state().total_velocity_m_per_s;
        stream << v.x() << ' ' << v.y() << ' ' << v.z() << ' ';
    }
    stream << "\n        </DataArray>\n";
    stream << "        <DataArray type=\"Float64\" Name=\"electric_field_V_per_cm\" NumberOfComponents=\"3\" "
              "format=\"ascii\">\n          ";
    for (const auto& p : m_particles) {
        const auto electric_field_V_per_cm = 0.01 * p.particle.state().electric_field_V_per_m;
        stream << electric_field_V_per_cm.x() << ' ' << electric_field_V_per_cm.y() << ' '
               << electric_field_V_per_cm.z() << ' ';
    }
    stream << "\n        </DataArray>\n";
    stream << "      </PointData>\n";
    stream << "      <Points>\n";
    stream << "        <DataArray type=\"Float64\" Name=\"position_um\" NumberOfComponents=\"3\" format=\"ascii\">\n   "
              "       ";
    for (const auto& p : m_particles) {
        const auto position_um = to_mesh_position_um(p.particle.state().position_m);
        stream << position_um.x() << ' ' << position_um.y() << ' ' << position_um.z() << ' ';
    }
    stream << "\n        </DataArray>\n";
    stream << "      </Points>\n";
    stream << "      <Verts>\n";
    stream << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n          ";
    for (std::size_t i = 0; i < n; ++i) {
        stream << static_cast<long long>(i) << ' ';
    }
    stream << "\n        </DataArray>\n";
    stream << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n          ";
    for (std::size_t i = 0; i < n; ++i) {
        stream << static_cast<long long>(i + 1) << ' ';
    }
    stream << "\n        </DataArray>\n";
    stream << "      </Verts>\n";
    stream << "    </Piece>\n";
    stream << "  </PolyData>\n";
    stream << "</VTKFile>\n";

    const auto already_recorded = std::find_if(m_particle_vtp_export_records.begin(),
                                               m_particle_vtp_export_records.end(),
                                               [&](const auto& record) { return record.filename == filename; });
    if (already_recorded == m_particle_vtp_export_records.end()) {
        m_particle_vtp_export_records.push_back({.time_s = m_state.m_time_s, .filename = filename});
    }
    write_particle_vtp_time_collection(pvd_path.string());
}

void device_admc_simulation::publish_mesh_particle_local_averages() const {
    auto* mesh = m_device.get_p_mesh();
    if (mesh == nullptr) {
        return;
    }

    const auto                                            list_bulk_elements = mesh->get_list_bulk_element();
    std::unordered_map<const mesh::element*, std::size_t> element_indices;
    element_indices.reserve(list_bulk_elements.size());
    for (std::size_t index = 0; index < list_bulk_elements.size(); ++index) {
        element_indices.emplace(list_bulk_elements[index].get(), index);
    }

    std::vector<double> particle_count(list_bulk_elements.size(), 0.0);
    std::vector<double> total_weight(list_bulk_elements.size(), 0.0);
    std::vector<double> mobility_sum(list_bulk_elements.size(), 0.0);
    std::vector<double> diffusion_sum(list_bulk_elements.size(), 0.0);
    std::vector<double> electric_field_norm_sum(list_bulk_elements.size(), 0.0);
    std::vector<double> speed_sum(list_bulk_elements.size(), 0.0);

    for (const auto& particle : m_particles) {
        const auto index_it = element_indices.find(particle.containing_element);
        if (index_it == element_indices.end()) {
            continue;
        }
        const std::size_t index = index_it->second;
        const auto&       state = particle.particle.state();

        particle_count[index] += 1.0;
        total_weight[index] += particle.weight;
        mobility_sum[index] += state.mobility_m2_per_V_s;
        diffusion_sum[index] += state.diffusion_m2_per_s;
        electric_field_norm_sum[index] += state.electric_field_V_per_m.norm();
        speed_sum[index] += state.total_velocity_m_per_s.norm();
    }

    std::vector<double> average_mobility(list_bulk_elements.size(), 0.0);
    std::vector<double> average_diffusion(list_bulk_elements.size(), 0.0);
    std::vector<double> average_electric_field_norm(list_bulk_elements.size(), 0.0);
    std::vector<double> average_speed(list_bulk_elements.size(), 0.0);
    for (std::size_t index = 0; index < list_bulk_elements.size(); ++index) {
        average_mobility[index]            = average_or_zero(mobility_sum[index], particle_count[index]);
        average_diffusion[index]           = average_or_zero(diffusion_sum[index], particle_count[index]);
        average_electric_field_norm[index] = average_or_zero(electric_field_norm_sum[index], particle_count[index]);
        average_speed[index]               = average_or_zero(speed_sum[index], particle_count[index]);
    }

    const auto publish_scalar_cell_field = [mesh, &element_indices](const std::string&         field_name,
                                                                    const std::vector<double>& values) {
        if (mesh->scalar_function_exists(field_name)) {
            mesh->remove_scalar_function(field_name);
        }

        for (const auto& region : mesh->get_list_bulk_region()) {
            const auto               region_elements = region.get_list_elements();
            std::vector<double>      region_values;
            std::vector<std::size_t> region_element_indices;
            region_values.reserve(region_elements.size());
            region_element_indices.reserve(region_elements.size());

            for (const auto& element : region_elements) {
                const auto index_it = element_indices.find(element.get());
                region_values.push_back(index_it == element_indices.end() ? 0.0 : values[index_it->second]);
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
    };
    publish_scalar_cell_field("particle_local_count", particle_count);
    publish_scalar_cell_field("particle_local_total_weight", total_weight);
    publish_scalar_cell_field("particle_local_average_mobility_m2_per_V_s", average_mobility);
    publish_scalar_cell_field("particle_local_average_diffusion_m2_per_s", average_diffusion);
    publish_scalar_cell_field("particle_local_average_electric_field_norm_V_per_m", average_electric_field_norm);
    publish_scalar_cell_field("particle_local_average_speed_m_per_s", average_speed);
}

void device_admc_simulation::export_current_mesh_as_vtu(const std::string& directory, bool export_x_cut_enabled) const {
    const std::filesystem::path output_directory(directory);
    std::filesystem::create_directories(output_directory);
    const std::filesystem::path vtu_path = output_directory / fmt::format("mesh_{:012d}.vtu", m_state.m_iteration);
    const std::string           filename = vtu_path.filename().string();
    const std::filesystem::path pvd_path = output_directory / "mesh.pvd";
    if (m_options.m_export_mesh_particle_local_averages) {
        publish_mesh_particle_local_averages();
    }
    file::export_as_vtu(*(m_device.get_p_mesh()), vtu_path.string());

    const auto already_recorded = std::find_if(m_mesh_vtu_export_records.begin(),
                                               m_mesh_vtu_export_records.end(),
                                               [&](const auto& record) { return record.filename == filename; });
    if (already_recorded == m_mesh_vtu_export_records.end()) {
        m_mesh_vtu_export_records.push_back({.time_s = m_state.m_time_s, .filename = filename});
    }
    write_mesh_vtu_time_collection(pvd_path.string());

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

void device_admc_simulation::export_current_snapshot() const {
    const std::filesystem::path base_directory(m_options.m_prefix_export_filename);
    export_current_mesh_as_vtu((base_directory / "mesh").string(), true);
    export_current_particles_as_vtp((base_directory / "particles").string());
    write_paraview_scene_script(base_directory);
}

void device_admc_simulation::write_particle_vtp_time_collection(const std::string& filename) const {
    write_vtk_time_collection(filename, m_particle_vtp_export_records);
}

void device_admc_simulation::write_mesh_vtu_time_collection(const std::string& filename) const {
    write_vtk_time_collection(filename, m_mesh_vtu_export_records);
}

}  // namespace uepm::ADMC
