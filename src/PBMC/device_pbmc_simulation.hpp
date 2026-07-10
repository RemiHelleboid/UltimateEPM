/**
 * @file device_pbmc_simulation.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-05-04
 *
 *
 */

#pragma once

#include <fstream>
#include <memory>
#include <random>
#include <vector>

#include "bbox.hpp"
#include "boundary_reflection.hpp"
#include "device.hpp"
#include "pbmc_device_history.hpp"
#include "pbmc_particle.hpp"
#include "pbmc_transport_kernel.hpp"
#include "vector.hpp"

namespace uepm::PBMC {

struct scheduled_particle_injection {
    double        m_time_s = 0.0;
    mesh::vector3 m_position_um{0.0, 0.0, 0.0};
    particle_type m_particle_type = particle_type::electron;
    double        m_weight        = 1.0;
    bool          m_done          = false;
};

struct current_probe_options {
    bool       m_enabled = false;
    mesh::bbox m_box_um{};
};

/**
 * @brief Struct containing the options of the admc device simulation.
 *
 */
struct options_device_PBMC {
    pbmc_material_model m_material_model = make_silicon_pbmc_material_model();

    std::string m_simulation_name  = "";
    std::string m_output_directory = "./";

    double      m_lattice_temperature                  = 300.0;
    double      m_max_energy_eV                        = 2.0;
    double      m_self_scattering_safety_factor        = 1.2;
    std::size_t m_gamma_max_energy_samples             = 1000;
    double      m_time_step                            = 1e-15;      // s
    double      m_t_max                                = 1e-9;       // s
    std::size_t m_max_number_particle                  = 100000000;  // Hard limit on nb of particles in the simulation.
    bool        m_particle_creation_activated          = true;
    bool        m_stop_simu_when_no_electron_remaining = true;
    bool        m_keep_particles_history               = false;
    bool        m_export_time_step                     = false;
    int         m_frequency_export_trajectory          = 10;
    int         m_nb_threads                           = 1;
    std::string m_prefix_export_filename               = "trajectory/time_step.csv";

    bool                         m_enable_scheduled_particle_injection = false;
    scheduled_particle_injection m_scheduled_particle_injection{};
    current_probe_options        m_current_probe{};

    bool                            m_activate_impact_ionization = true;
    bool                            m_enable_impurity_scattering = false;
    impurity_scattering_model       m_impurity_scattering_model  = impurity_scattering_model::mobility_empirical;
    impurity_screening_model        m_impurity_screening_model   = impurity_screening_model::debye_analytic;
    mesh::boundary_reflection_model m_boundary_reflection_model  = mesh::boundary_reflection_model::reverse;

    options_device_PBMC() = default;

    options_device_PBMC(double      t_max,
                        double      time_step,
                        std::size_t max_number_particle,
                        bool        activate_impact_ionization,
                        bool        particle_creation_activated,
                        bool        stop_simu_when_no_electron_remaining,
                        bool        keep_particles_history,
                        bool        export_time_step,
                        int         frequency_export_trajectory,
                        int         nb_threads = 1)
        : m_time_step(time_step),
          m_t_max(t_max),
          m_max_number_particle(max_number_particle),
          m_particle_creation_activated(particle_creation_activated),
          m_stop_simu_when_no_electron_remaining(stop_simu_when_no_electron_remaining),
          m_keep_particles_history(keep_particles_history),
          m_export_time_step(export_time_step),
          m_frequency_export_trajectory(frequency_export_trajectory),
          m_nb_threads(nb_threads),
          m_activate_impact_ionization(activate_impact_ionization) {}

    void validate() const;

    void print_options() const {
        std::cout << "Simulation options: " << std::endl;
        std::cout << "Time step: " << m_time_step << std::endl;
        std::cout << "Final time: " << m_t_max << std::endl;
        std::cout << "Max number of particles: " << m_max_number_particle << std::endl;
        std::cout << "Activate impact ionization: " << m_activate_impact_ionization << std::endl;
        std::cout << "Activate particle creation: " << m_particle_creation_activated << std::endl;
        std::cout << "Stop simulation when no electron remaining: " << m_stop_simu_when_no_electron_remaining
                  << std::endl;
        std::cout << "Keep particles history: " << m_keep_particles_history << std::endl;
        std::cout << "Export time step: " << m_export_time_step << std::endl;
        std::cout << "Frequency export trajectory: " << m_frequency_export_trajectory << std::endl;
        std::cout << "Number of threads: " << m_nb_threads << std::endl;
    }
};

struct state_device_pbmc_simulation {
    double      m_time_s                    = 0.0;
    std::size_t m_iteration                 = 0;
    std::size_t m_counter_particles_created = 0;

    bool    m_scheduled_particle_injection_done     = true;
    bool    m_use_constant_RamoUnitaryElectricField = false;
    vector3 m_RamoUnitaryElectricField_Vm_per_cm{0.0, 0.0, 0.0};
};

struct vtk_time_series_record {
    double      m_time_s = 0.0;
    std::string m_filename;
};

/**
 * @brief Device PBMC simulation class. This class contains the main loop of the simulation and the list of particles.
 *
 */
class device_pbmc_simulation {
    struct ramo_current_components {
        double electron       = 0.0;
        double hole           = 0.0;
        double probe_electron = 0.0;
        double probe_hole     = 0.0;
    };

 protected:
    state_device_pbmc_simulation       m_state;
    device::device                     m_device;
    pbmc_transport_kernel              m_electron_transport;
    pbmc_transport_kernel              m_hole_transport;
    std::vector<pbmc_transport_kernel> m_thread_electron_transports;
    std::vector<pbmc_transport_kernel> m_thread_hole_transports;
    int                                m_dimension;
    options_device_PBMC                m_simulation_options;
    history_device_PBMC                m_simulation_history{};
    std::minstd_rand                   m_boundary_reflection_rng;

    std::vector<std::unique_ptr<pbmc_particle>>  m_list_particles;
    std::vector<std::optional<scattering_event>> m_scattering_events_scratch;

    void                         initialize_scheduled_particle_injection();
    bool                         has_pending_scheduled_particle_injection() const;
    void                         inject_scheduled_particle_if_due();
    void                         validate_time_step_against_scattering_rate() const;
    static pbmc_transport_config make_transport_config(const options_device_PBMC &options, particle_type carrier_type);
    void                         initialize_thread_transports(int seed_random_generator);
    pbmc_transport_kernel       &transport_for(particle_type type);
    const pbmc_transport_kernel &transport_for(particle_type type) const;
    pbmc_transport_kernel       &transport_for(particle_type type, std::size_t thread_index);
    void                         initialize_particle_transport_state(pbmc_particle &particle);
    std::string                  initialize_simulation_history_file();
    void                         flatten_particle_positions_for_2d();

    // Export functions

    mutable std::vector<vtk_time_series_record> m_particle_vtp_export_records;
    mutable std::vector<vtk_time_series_record> m_mesh_vtk_export_records;

    void export_current_particles_as_vtp(const std::string &directory) const;
    void write_particle_vtp_time_collection(const std::string &pvd_filename) const;
    void export_current_mesh_as_vtk(const std::string &directory) const;
    void publish_mesh_particle_local_average_energy() const;
    void publish_mesh_particle_local_current_density() const;
    void write_mesh_vtk_time_collection(const std::string &pvd_filename) const;
    void export_current_snapshot() const;

    vector3        get_RamoUnitaryElectricField_at_position(const mesh::vector3 &position) const;
    vector3        get_RamoUnitaryElectricField_at_position(const mesh::vector3 &position,
                                                            const mesh::element *containing_element) const;
    virtual double current_density_cell_volume_m3(const mesh::element &element) const;

 public:
    /**
     * @brief Construct an ampty simulation object.
     *
     * @param simulation_device
     * @param simulation_option
     * @param seed_random_generator
     */
    device_pbmc_simulation(const device::device      &simulation_device,
                           const options_device_PBMC &simulation_option,
                           int                        seed_random_generator = 0);

    /**
     * @brief Construct a new device admc simulation object
     *
     * @param simulation_device
     * @param simulation_option
     * @param starting_position
     * @param number_electrons_start
     * @param number_holes_start
     */
    device_pbmc_simulation(const device::device      &simulation_device,
                           const options_device_PBMC &simulation_option,
                           const mesh::vector3       &starting_position,
                           std::size_t                number_electrons_start,
                           std::size_t                number_holes_start,
                           int                        seed_random_generator = 0);

    void add_particle_at_position(const mesh::vector3 &location, particle_type type_of_particle, double weight = 1.0);
    void add_particles_at_positions(const std::vector<mesh::vector3> &positions,
                                    particle_type                     type_of_particle,
                                    double                            weight = 1.0);
    std::size_t load_particles_from_state_csv(const std::string &filename);

    void transport_particles_one_time_step();
    void advance_particles_one_time_step();

    void set_particles_transport_data_from_device();
    void update_element_and_check_boundary();

    void                      remove_collected_particles();
    ramo_current_components   compute_ramo_currents(bool include_full, bool include_probe) const;
    std::pair<double, double> compute_ramo_current() const;
    std::pair<double, double> compute_probe_ramo_current() const;
    double                    compute_ramo_current_for_particle(const pbmc_particle &particle) const;
    virtual double            ramo_current_scale_factor() const;

    void run();

    void set_max_number_particles(std::size_t new_value) { m_simulation_options.m_max_number_particle = new_value; }
    void set_keep_particles_history(bool new_value) { m_simulation_options.m_keep_particles_history = new_value; }
    void set_exporting_iterations(bool new_value) { m_simulation_options.m_export_time_step = new_value; }
    void set_exporting_frequency(int new_value) { m_simulation_options.m_frequency_export_trajectory = new_value; }
    void set_prefix_export_trajectory_filename(const std::string &new_prefix) {
        m_simulation_options.m_prefix_export_filename = new_prefix;
    }
    void set_stop_simulation_without_electron(bool new_value) {
        m_simulation_options.m_stop_simu_when_no_electron_remaining = new_value;
    }
    bool has_reached_particle_limit() const {
        return m_list_particles.size() >= m_simulation_options.m_max_number_particle;
    }

    const history_device_PBMC &get_simulation_history() const { return m_simulation_history; }

    std::size_t                get_number_electrons() const;
    std::size_t                get_number_holes() const;
    double                     get_total_electron_weight() const;
    double                     get_total_hole_weight() const;
    std::optional<double>      get_current_time() const { return m_state.m_time_s; }
    std::vector<mesh::vector3> get_all_particles_position() const;

    /**
     * @brief Compute the depletion region of the device (x_min, x_max).
     * Works only for device with 1D symmetry.
     *
     * @return std::pair<double, double>
     */
    std::pair<double, double> compute_depletion_region() const;

    void export_history_to_csv(const std::string &filename, std::size_t frequency = 1) {
        m_simulation_history.export_to_csv(filename, frequency);
    }

    /**
     * @brief Export current state of the simulation in a single csv file with all particle positions and data.
     *
     * @param prefix_filename
     */
    void export_current_time_step_as_csv(const std::string &prefix_filename) const;
    void export_particle_state_csv(const std::string &filename) const;

    /**
     * @brief Export all trajectories in a single csv files.
     * Works only if the particles saved their history.
     *
     * @param prefix_filename
     */
    void export_all_trajectories_as_csv(const std::string &prefix_filename) const;
};

}  // namespace uepm::PBMC
