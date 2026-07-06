/**
 * @file device_admc_simulation.hpp
 * @brief Device-level ADMC drift-diffusion particle simulation.
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <map>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "admc_transport.hpp"
#include "boundary_reflection.hpp"
#include "bbox.hpp"
#include "device.hpp"

namespace uepm::ADMC {

struct current_probe_options {
    bool       m_enabled = false;
    mesh::bbox m_box_um{};
};

struct options_device_ADMC {
    std::string m_simulation_name;
    std::string m_output_directory = "./";

    double      m_lattice_temperature_K               = 300.0;
    double      m_time_step_s                         = 1.0e-15;
    double      m_final_time_s                        = 1.0e-12;
    std::size_t m_max_number_particles                = 100000000;
    bool        m_stop_when_no_electrons              = true;
    bool        m_export_time_step                    = false;
    int         m_frequency_export                    = 10;
    bool        m_export_mesh_particle_local_averages = true;
    std::string m_prefix_export_filename              = "trajectory";
    mesh::boundary_reflection_model m_boundary_reflection_model = mesh::boundary_reflection_model::reverse;
    current_probe_options m_current_probe{};

    bool          m_enable_scheduled_particle_injection = false;
    double        m_scheduled_injection_time_s          = 0.0;
    mesh::vector3 m_scheduled_injection_position_um{0.0, 0.0, 0.0};
    carrier_type  m_scheduled_injection_type   = carrier_type::electron;
    double        m_scheduled_injection_weight = 1.0;

    void validate() const;
};

struct device_admc_particle {
    admc_particle  particle;
    mesh::element* containing_element = nullptr;
    double         weight             = 1.0;
    bool           crossed_contact    = false;
};

struct state_device_ADMC {
    double      m_time_s                                = 0.0;
    std::size_t m_iteration                             = 0;
    std::size_t m_counter_particles_created             = 0;
    bool        m_scheduled_particle_injection_done     = true;
    bool        m_use_constant_RamoUnitaryElectricField = false;
    vector3     m_RamoUnitaryElectricField_Vm_per_cm{0.0, 0.0, 0.0};
    double      m_last_ramo_current_electron_A = 0.0;
    double      m_last_ramo_current_hole_A     = 0.0;
    double      m_last_probe_ramo_current_electron_A = 0.0;
    double      m_last_probe_ramo_current_hole_A     = 0.0;
};

struct admc_vtk_time_series_record {
    double      time_s = 0.0;
    std::string filename;
};

struct history_device_ADMC {
    std::vector<double>      times_s;
    std::vector<std::size_t> nb_electrons;
    std::vector<std::size_t> nb_holes;
    std::vector<double>      ramo_current_electron_A;
    std::vector<double>      ramo_current_hole_A;
    std::vector<double>      ramo_current_A;
    std::vector<double>      probe_ramo_current_electron_A;
    std::vector<double>      probe_ramo_current_hole_A;
    std::vector<double>      probe_ramo_current_A;
    std::vector<double>      max_electric_field_V_per_m;
    std::vector<std::string> contact_voltage_names;
    std::vector<std::vector<double>> contact_voltages_V;

    void add(double      time_s,
             std::size_t electrons,
             std::size_t holes,
             double      electron_current_A,
             double      hole_current_A,
             double      total_current_A,
             double      probe_electron_current_A,
             double      probe_hole_current_A,
             double      probe_total_current_A,
             double      max_field_V_per_m,
             const std::vector<double>& active_contact_voltages_V = {});
    void set_contact_voltage_names(const std::vector<std::string>& contact_names);
    std::vector<double> contact_voltage_values_from_map(const std::map<std::string, double>& active_contact_voltages_V)
        const;
    void set_last_contact_voltages(const std::vector<double>& active_contact_voltages_V);
    void print_header_csv(const std::string& filename) const;
    void append_last_iter_to_csv(std::fstream& file) const;
    void export_to_csv(const std::string& filename) const;
};

class device_admc_simulation {
 public:
    device_admc_simulation(const device::device&      simulation_device,
                           const options_device_ADMC& options,
                           std::uint64_t              random_seed = 5489u);

    device_admc_simulation(const device::device&      simulation_device,
                           const options_device_ADMC& options,
                           const mesh::vector3&       starting_position_um,
                           std::size_t                number_electrons_start,
                           std::size_t                number_holes_start,
                           std::uint64_t              random_seed = 5489u);

    void add_particle_at_position(const mesh::vector3& position_um, carrier_type type, double weight = 1.0);
    void advance_particles_one_time_step();
    void run();

    std::size_t get_number_electrons() const;
    std::size_t get_number_holes() const;
    double      get_total_electron_weight() const;
    double      get_total_hole_weight() const;
    double      current_time_s() const noexcept { return m_state.m_time_s; }

    std::pair<double, double> compute_ramo_current() const;
    std::pair<double, double> compute_probe_ramo_current() const;
    double                    max_particle_electric_field_V_per_m() const;

    const std::vector<device_admc_particle>& particles() const noexcept { return m_particles; }
    const history_device_ADMC&               history() const noexcept { return m_history; }
    void export_history_to_csv(const std::string& filename) const { m_history.export_to_csv(filename); }
    void export_current_time_step_as_csv(const std::string& prefix_filename) const;
    void export_current_snapshot() const;
    void set_prefix_export_trajectory_filename(const std::string& new_prefix) {
        m_options.m_prefix_export_filename = new_prefix;
    }

 protected:
    mesh::vector3 to_mesh_position_um(const vector3& position_m) const;
    vector3       to_admc_position_m(const mesh::vector3& position_um) const;
    mesh::vector3 normalize_mesh_position_for_dimension(mesh::vector3 position_um) const;
    bool          is_transport_material_element(mesh::element& element);
    vector3       draw_standard_normal();

    admc_local_environment    local_environment(const device_admc_particle& particle) const;
    void                      initialize_particle_device_state(device_admc_particle& particle);
    void                      update_element_and_check_boundary(device_admc_particle& particle);
    void                      remove_collected_particles();
    void                      initialize_scheduled_particle_injection();
    bool                      has_pending_scheduled_particle_injection() const;
    void                      inject_scheduled_particle_if_due();
    void                      record_history(double electron_current_A, double hole_current_A);
    std::string               initialize_simulation_history_file();
    std::pair<double, double> last_ramo_current() const;
    std::pair<double, double> last_probe_ramo_current() const;
    vector3                   get_RamoUnitaryElectricField_at_position(const mesh::vector3& position) const;
    void                      export_current_particles_as_vtp(const std::string& directory) const;
    void export_current_mesh_as_vtu(const std::string& directory, bool export_x_cut_enabled = false) const;
    void publish_mesh_particle_local_averages() const;
    void write_particle_vtp_time_collection(const std::string& filename) const;
    void write_mesh_vtu_time_collection(const std::string& filename) const;
    bool has_reached_particle_limit() const { return m_particles.size() >= m_options.m_max_number_particles; }

    state_device_ADMC                                m_state;
    device::device                                   m_device;
    options_device_ADMC                              m_options;
    admc_transport_kernel                            m_transport;
    int                                              m_dimension = 3;
    std::vector<device_admc_particle>                m_particles;
    history_device_ADMC                              m_history;
    mutable std::vector<admc_vtk_time_series_record> m_particle_vtp_export_records;
    mutable std::vector<admc_vtk_time_series_record> m_mesh_vtu_export_records;
    std::mt19937_64                                  m_random_generator;
    std::normal_distribution<double>                 m_standard_normal{0.0, 1.0};
};

}  // namespace uepm::ADMC
