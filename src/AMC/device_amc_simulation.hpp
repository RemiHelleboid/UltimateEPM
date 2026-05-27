/**
 * @file device_amc_simulation.hpp
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
#include <vector>

#include "amc_device_history.hpp"
#include "amc_transport_kernel.hpp"
#include "device.hpp"
#include "particle_amc.hpp"
#include "vector.hpp"

namespace uepm::amc {

/**
 * @brief Struct containing the options of the admc device simulation.
 *
 */
struct options_device_amc {
    double      m_lattice_temperature                  = 300.0;
    double      m_max_energy_eV                        = 2.0;
    double      m_self_scattering_safety_factor        = 1.2;
    std::size_t m_gamma_max_energy_samples             = 1000;
    double      m_time_step                            = 1e-15;    // s
    double      m_t_max                                = 1e-9;     // s
    std::size_t m_max_number_particle                  = 1000000;  // Hard limit on nb of particles in the simulation.
    std::size_t m_avalanche_threshold                  = 1000;
    bool        m_activate_impact_ionization           = true;
    bool        m_particle_creation_activated          = true;
    bool        m_stop_simu_when_no_electron_remaining = true;
    bool        m_keep_particles_history               = false;
    bool        m_export_time_step                     = false;
    int         m_frequency_export_trajectory          = 10;
    int         m_nb_threads                           = 1;

    options_device_amc() = default;

    options_device_amc(double      t_max,
                       double      time_step,
                       std::size_t max_number_particle,
                       std::size_t avalanche_threshold,
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
          m_avalanche_threshold(avalanche_threshold),
          m_activate_impact_ionization(activate_impact_ionization),
          m_particle_creation_activated(particle_creation_activated),
          m_stop_simu_when_no_electron_remaining(stop_simu_when_no_electron_remaining),
          m_keep_particles_history(keep_particles_history),
          m_export_time_step(export_time_step),
          m_frequency_export_trajectory(frequency_export_trajectory),
          m_nb_threads(nb_threads) {}

    void print_options() const {
        std::cout << "Simulation options: " << std::endl;
        std::cout << "Time step: " << m_time_step << std::endl;
        std::cout << "Final time: " << m_t_max << std::endl;
        std::cout << "Max number of particles: " << m_max_number_particle << std::endl;
        std::cout << "Avalanche threshold: " << m_avalanche_threshold << std::endl;
        std::cout << "Activate impact ionization: " << m_activate_impact_ionization << std::endl;
        std::cout << "Activate particle creation: " << m_particle_creation_activated << std::endl;
        std::cout << "Stop simulation when no electron remaining: " << m_stop_simu_when_no_electron_remaining << std::endl;
        std::cout << "Keep particles history: " << m_keep_particles_history << std::endl;
        std::cout << "Export time step: " << m_export_time_step << std::endl;
        std::cout << "Frequency export trajectory: " << m_frequency_export_trajectory << std::endl;
        std::cout << "Number of threads: " << m_nb_threads << std::endl;
    }
};

/**
 * @brief Device amc simulation class. This class contains the main loop of the simulation and the list of particles.
 * 
 */
class device_amc_simulation {
 protected:
    device::device       m_device;
    amc_transport_kernel m_electron_transport;
    amc_transport_kernel m_hole_transport;
    std::string          m_simulation_name = "";
    int                  m_dimension;
    options_device_amc   m_simulation_options;
    history_device_amc   m_simulation_history{};
    double               m_time = 0.0;
    std::size_t          m_iteration;

    std::vector<std::unique_ptr<particle_amc>> m_list_particles;

    double m_anode_current{0.0};
    double m_cathode_current{0.0};

    std::string m_prefix_export_filename = "trajectory/time_step.csv";

    static amc_transport_config make_transport_config(const options_device_amc &options, particle_type carrier_type);
    amc_transport_kernel       &transport_for(particle_type type);
    const amc_transport_kernel &transport_for(particle_type type) const;
    void                        initialize_particle_transport_state(particle_amc &particle);

    virtual void apply_z_periodicity_to_particles();

 public:
    /**
     * @brief Construct an ampty simulation object.
     *
     * @param simulation_device
     * @param simulation_option
     */
    device_amc_simulation(const device::device     &simulation_device,
                          const options_device_amc &simulation_option,
                          const std::string        &simulation_name       = "",
                          int                       seed_random_generator = 0);

    /**
     * @brief Construct a new device admc simulation object
     *
     * @param simulation_device
     * @param simulation_option
     * @param starting_position
     * @param number_electrons_start
     * @param number_holes_start
     */
    device_amc_simulation(const device::device     &simulation_device,
                          const options_device_amc &simulation_option,
                          const std::string        &simulation_name,
                          const mesh::vector3      &starting_position,
                          std::size_t               number_electrons_start,
                          std::size_t               number_holes_start,
                          int                       seed_random_generator = 0);

    void               set_simulation_name(const std::string &new_name) { m_simulation_name = new_name; }
    const std::string &get_simulation_name() const { return m_simulation_name; }


    void add_particle_at_position(const mesh::vector3 &location, particle_type type_of_particle, double weight = 1.0);
    void add_particles_at_positions(const std::vector<mesh::vector3> &positions, particle_type type_of_particle, double weight = 1.0);

    void transport_particles_one_time_step();
    void advance_particles_one_time_step();
    
    void set_particles_transport_data_from_device();
    void update_element_and_check_boundary();
    
    void remove_collected_particles();
    double compute_ramo_current() const;
    
    void run();

    void set_max_number_particles(std::size_t new_value) { m_simulation_options.m_max_number_particle = new_value; }
    void set_keep_particles_history(bool new_value) { m_simulation_options.m_keep_particles_history = new_value; }
    void set_exporting_iterations(bool new_value) { m_simulation_options.m_export_time_step = new_value; }
    void set_exporting_frequency(int new_value) { m_simulation_options.m_frequency_export_trajectory = new_value; }
    void set_prefix_export_trajectory_filename(const std::string &new_prefix) { m_prefix_export_filename = new_prefix; }
    void set_stop_simulation_without_electron(bool new_value) { m_simulation_options.m_stop_simu_when_no_electron_remaining = new_value; }
    bool has_reached_avalanche() const { return m_list_particles.size() >= m_simulation_options.m_max_number_particle; }

    /**
     * @brief Return a const reference to the simulation history.
     *
     * @return const history_device_amc&
     */
    const history_device_amc &get_simulation_history() const { return m_simulation_history; }

    std::size_t get_number_electrons() const;
    std::size_t get_number_holes() const;
    std::optional<double> get_current_time() const { return m_time; }

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

    /**
     * @brief Export all trajectories in a single csv files.
     * Works only if the particles saved their history.
     *
     * @param prefix_filename
     */
    void export_all_trajectories_as_csv(const std::string &prefix_filename) const;
};

}  // namespace uepm::amc