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
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include "intervalley_phonon.hpp"
#include "particle_amc.hpp"
#include "scattering_channels.hpp"
#include "valley_model.hpp"
#include "vector.hpp"
#include "device.hpp"
#include "physical_constants.hpp"
#include "materials.hpp"

namespace uepm::amc {

/**
 * @brief Struct containing the options of the admc device simulation.
 *
 */
struct options_device_admc {
    /**
     * @brief Time step used for the simulation.
     *
     */
    double m_time_step;

    /**
     * @brief Final time of the simulation.
     *
     */
    double m_t_max;

    /**
     * @brief Maximum number of active particles. If this number is reached, the simulation stops.
     *
     */
    std::size_t m_max_number_particle{0};

    /**
     * @brief Threshold number of particle above wich we consider to have an avalanche.
     *
     */
    std::size_t m_avalanche_threshold{0};

    /**
     * @brief If True, in the simulation, impact ionization will be computed, and generations will occurr if m_particle_creation_activated
     * is true.
     *
     */
    bool m_activate_impact_ionization = true;

    /**
     * @brief If m_activate_impact_ionization is true and m_activate_particle is false, impact ionization process will be computed BUT
     * no e-h pairs will be created.
     *
     */
    bool m_particle_creation_activated = true;

    /**
     * @brief If this parameter is set to true, the simulation stops when there is no electron remaining in the device.
     * It is usefull to improve simulation time.
     *
     */
    bool m_stop_simu_when_no_electron_remaining = true;

    /**
     * @brief If true, all particles will keep their history of position, impact-ionization, time, etc.
     *
     */
    bool m_keep_particles_history = false;

    /**
     * @brief If true, all iterations data will be exported in a csv file.
     *
     */
    bool m_export_time_step = false;

    /**
     * @brief Frequency in time step at which the trajectory is exported.
     *
     */
    int m_frequency_export_trajectory = 10;

    /**
     * @brief Number of threads used for the simulation. (OpenMP)
     *
     */
    int m_nb_threads = 1;

    options_device_admc(double      t_max,
                        double      time_step,
                        std::size_t max_number_particle,
                        std::size_t avalanche_threshod,
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
          m_avalanche_threshold(avalanche_threshod),
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
 * @brief Structure contatining the history of the device_simulation states.
 *
 */
struct history_device_admc {
    mesh::vector3              m_last_impact_ionization_position{};
    std::vector<double>        m_list_times{};
    std::vector<std::size_t>   m_list_nb_electrons{};
    std::vector<std::size_t>   m_list_nb_holes{};
    std::vector<std::size_t>   m_list_nb_impact_ionization{};
    std::vector<double>        m_list_anode_current{};
    std::vector<double>        m_list_cathode_current{};
    std::vector<mesh::vector3> m_impact_ionization_positions{};
    std::vector<double>        m_max_electric_field{};
    std::size_t                m_initial_seed_rng{0};

    history_device_admc() = default;

    void reserve_memory(std::size_t size) {
        m_list_times.reserve(size);
        m_list_nb_electrons.reserve(size);
        m_list_nb_holes.reserve(size);
        m_list_nb_impact_ionization.reserve(size);
        m_list_anode_current.reserve(size);
        m_list_cathode_current.reserve(size);
        m_max_electric_field.reserve(size);
    }

    void add_data_to_history(double      time,
                             std::size_t nb_electrons,
                             std::size_t nb_holes,
                             std::size_t nb_impact_ionization,
                             double      anode_current,
                             double      cathode_current,
                             double      max_electric_field) {
        m_list_times.push_back(time);
        m_list_nb_electrons.push_back(nb_electrons);
        m_list_nb_holes.push_back(nb_holes);
        m_list_nb_impact_ionization.push_back(nb_impact_ionization);
        m_list_anode_current.push_back(anode_current);
        m_list_cathode_current.push_back(cathode_current);
    }

    void print_header_csv(const std::string &filename) {
        std::ofstream file(filename);
        file << "time,nb_electrons,nb_holes,nb_impact_ionization,anode_current,cathode_current,max_electric_field\n";
        file.close();
    }

    void append_last_iter_to_csv(std::fstream &file) {
        file << m_list_times.back() << ',' << m_list_nb_electrons.back() << ',' << m_list_nb_holes.back() << ','
             << m_list_nb_impact_ionization.back() << ',' << m_list_anode_current.back() << ',' << m_list_cathode_current.back() << ','
             << m_max_electric_field.back() << '\n';
    }

    void export_to_csv(const std::string &filename, std::size_t frequency = 1) {
        std::vector<double> double_list_time;
        std::vector<double> double_list_nb_electrons;
        std::vector<double> double_list_nb_hole;
        std::vector<double> double_list_nb_impact_ionization;
        std::vector<double> double_list_anode_current;
        std::vector<double> double_list_cathode_current;
        std::vector<double> double_list_max_electric_field;
        for (std::size_t iter_nb = 0; iter_nb < m_list_times.size() - 1; iter_nb++) {
            double_list_time.push_back(m_list_times[iter_nb]);
            double_list_nb_electrons.push_back(m_list_nb_electrons[iter_nb]);
            double_list_nb_hole.push_back(m_list_nb_holes[iter_nb]);
            double_list_nb_impact_ionization.push_back(m_list_nb_impact_ionization[iter_nb]);
            double_list_anode_current.push_back(m_list_anode_current[iter_nb]);
            double_list_cathode_current.push_back(m_list_cathode_current[iter_nb]);
            double_list_max_electric_field.push_back(m_max_electric_field[iter_nb]);
        }
        // Always add the last iteration
        double_list_time.push_back(m_list_times[m_list_nb_electrons.size() - 1]);
        double_list_nb_electrons.push_back(m_list_nb_electrons[m_list_nb_electrons.size() - 1]);
        double_list_nb_hole.push_back(m_list_nb_holes[m_list_nb_electrons.size() - 1]);
        double_list_nb_impact_ionization.push_back(m_list_nb_impact_ionization[m_list_nb_electrons.size() - 1]);
        double_list_anode_current.push_back(m_list_anode_current[m_list_nb_electrons.size() - 1]);
        double_list_cathode_current.push_back(m_list_cathode_current[m_list_nb_electrons.size() - 1]);
        double_list_max_electric_field.push_back(m_max_electric_field[m_list_nb_electrons.size() - 1]);

        std::vector<std::string> header_csv =
            {"time", "nb_electrons", "nb_holes", "nb_impact_ionization", "anode_current", "cathode_current", "max_electric_field"};
        utils::export_multiple_vector_to_csv(filename,
                                             header_csv,
                                             {double_list_time,
                                              double_list_nb_electrons,
                                              double_list_nb_hole,
                                              double_list_nb_impact_ionization,
                                              double_list_anode_current,
                                              double_list_cathode_current,
                                              double_list_max_electric_field});
    }

    std::vector<std::size_t> get_history_total_nb_particles() const {
        std::vector<std::size_t> history_total_number_particles(m_list_nb_electrons.size());
        for (std::size_t index_iteration = 0; index_iteration < m_list_nb_electrons.size(); ++index_iteration) {
            history_total_number_particles[index_iteration] = m_list_nb_electrons[index_iteration] + m_list_nb_holes[index_iteration];
        }
        return history_total_number_particles;
    }

    std::optional<double> get_time_at_number_particle(std::size_t nb_particle) const;
};

/**
 * @brief Bulk Advection-Diffusion Monte Carlo Simulation Class
 *
 * This class handle bulk simulation with ADMC Monte Carlo.
 * It means that we suppose that electrons and holes are in the full R^3 space, without border etc.
 * The electric field and the doping profile are specified through "analytical" function that are given as input.
 *
 */
class device_admc_simulation {
 protected:
    /**
     * @brief Device in which the simulation will be performed.
     *
     */
    device::device m_device;

    /**
     * @brief Name of the simulation.
     *
     */
    std::string m_simulation_name = "";

    /**
     * @brief Dimension of the simulation.
     *
     */
    int m_dimension;

    /**
     * @brief Structure containing the options of the simulation
     *
     */
    options_device_admc m_simulation_options;

    /**
     * @brief History of the simulation.
     *
     */
    history_device_admc m_simulation_history{};

    /**
     * @brief Time of the simulation.
     *
     */
    double m_time = 0.0;

    /**
     * @brief Number of iterations performed.
     *
     */
    std::size_t m_iteration;

    /**
     * @brief Time at which the number of particle for avalanche threshold is reached.
     *
     */
    std::optional<double> m_time_to_avalanche{};

    /**
     * @brief Time at which the number of particle for avalanche threshold is reached, but for a secondary avalanche (AfterPulsing).
     *
     */
    std::optional<double> m_time_to_secondary_avalanche{};

    /**
     * @brief Vector of unique_ptr on the particles of the simulation.
     *
     */
    std::vector<std::unique_ptr<particle_amc>> m_list_particles;

    /**
     * @brief Total number of impact ionization.
     *
     */
    std::size_t m_number_impact_ionizations{0};

    /**
     * @brief Current at the anode of the device.
     *
     */
    double m_anode_current{0.0};

    /**
     * @brief Current at the anode of the device.
     *
     */
    double m_cathode_current{0.0};

    /**
     * @brief Maximum electric field at this iteration.
     *
     */
    double m_max_global_electric_field{0.0};

    /**
     * @brief Vector of unique_ptr on the particles of the simulation.
     *
     */
    std::vector<std::unique_ptr<particle_amc>> m_list_deleted_particles;

    /**
     * @brief Mersenne random number generator for the Random Path Length Algorithm and the brownian motion transport.
     *
     */
    std::minstd_rand m_random_generator;

    /**
     * @brief Random normal distribution for the brownian motion simulation.
     *
     */
    std::normal_distribution<double> m_normal_distribution;

    /**
     * @brief Uniform in [0, 1] distribution for the Random Path Length Algorithm
     *
     */
    std::uniform_real_distribution<double> m_uniform_distribution;

    /**
     * @brief Prefix path of location where the iterations are saved.
     *
     */
    std::string m_prefix_export_filename = "trajectory/time_step.csv";

 public:
    /**
     * @brief Construct an ampty simulation object.
     *
     * @param simulation_device
     * @param simulation_option
     */
    device_admc_simulation(const device::device      &simulation_device,
                           const options_device_admc &simulation_option,
                           const std::string         &simulation_name       = "",
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
    device_admc_simulation(const device::device      &simulation_device,
                           const options_device_admc &simulation_option,
                           const std::string         &simulation_name,
                           const mesh::vector3       &starting_position,
                           std::size_t                number_electrons_start,
                           std::size_t                number_holes_start,
                           int                        seed_random_generator = 0);

    void seed_random_generator(int new_seed) {
        m_random_generator.seed(new_seed);
        m_simulation_history.m_initial_seed_rng = new_seed;
    }

    void               set_simulation_name(const std::string &new_name) { m_simulation_name = new_name; }
    const std::string &get_simulation_name() const { return m_simulation_name; }


    /**
     * @brief Add a particle at a given position.
     *
     * @param location
     * @param type_of_particle
     */
    void add_particle_at_position(const mesh::vector3 &location, particle_type type_of_particle, double weight = 1.0);

    /**
     * @brief Add multiple particles at given positions.
     *
     * @param positions
     * @param type_of_particle
     */
    void add_particles_at_positions(const std::vector<mesh::vector3> &positions, particle_type type_of_particle, double weight = 1.0);

    /**
     * @brief Set transport data of the particles from the device (drift velocity and diffusion coefficient).
     * This function is called at the beginning of each iteration to update the transport data of the particles.
     *
     * @param device
     */
    void set_particles_transport_data_from_device();

    void perform_drift_diffusion_transport_without_data_setting();

    /**
     * @brief Move all particles according to the ADMC algorithm.
     * Only transport is performed, nothing about impact ionization or so.
     *
     */
    void perform_drift_diffusion_transport();

    /**
     * @brief For each particle, perform the step of the random path length algorithm to compute potential impact ionization process.
     *
     */
    void perform_impact_ionization(double weight_poisson_new_particle = 1.0);

    /**
     * @brief Run the main loop of transport simulation, without any impact ionization process.
     *
     */
    void run_transport_simulation();

    /**
     * @brief Run the main loop of the transport simulation with impact ionization process.
     *
     */
    void run_transport_ionization_simulation();

    /**
     * @brief Run the main loop of the transport simulation for a single particle only. Returns true if the particle reaches the contact
     * or a region with a doping concentration that exceeds a set doping threshold (note that we might have problems with signs in the
     * future if we want to investigate -ive doping concentration ie PP). Returns false if the simulation exceeds max time.
     *
     */
    bool run_transport_to_doping_threshold(const double doping_threshold);

    /**
     * @brief Update the m_containing_elements of the particles and check if a particle
     * change of region or went outside of the mesh.
     *
     */
    void update_element_and_check_boundary();

    /**
     * @brief Return true if the segment [point_A, point_B] is crossing one of the device contact.
     *
     * @param point_A
     * @param point_B
     * @return true
     * @return false
     */
    bool check_crossing_contact(const mesh::vector3 &point_A, const mesh::vector3 &point_B) const;

    /**
     * @brief Remove particle that was spotted at crossing the contacts.
     *
     */
    void remove_collected_particles();

    /**
     * @brief Compute current through Ramu-like formula.
     *
     * @return double
     */
    double compute_ramo_current() const;

    /**
     * @brief Set the maximum number of particles in the simulation.
     * This is also the avalanche threshold.
     *
     * @param new_value
     */
    void set_max_number_particles(std::size_t new_value) { m_simulation_options.m_max_number_particle = new_value; }

    /**
     * @brief Set the keep particles history object
     *
     * @param new_value
     */
    void set_keep_particles_history(bool new_value) { m_simulation_options.m_keep_particles_history = new_value; }

    /**
     * @brief Set the status of exporting each iteration. Set to true to save the particles location and data after each step.
     *
     * @param new_value
     */
    void set_exporting_iterations(bool new_value) { m_simulation_options.m_export_time_step = new_value; }

    /**
     * @brief Set the exporting frequency.
     * If different to 1, and if m_export_time_step is true, then the export occurs only each frequency iteration step.
     *
     * @param new_value
     */
    void set_exporting_frequency(int new_value) { m_simulation_options.m_frequency_export_trajectory = new_value; }

    /**
     * @brief Set the prefix export trajectory filename object
     *
     * @param new_prefix
     */
    void set_prefix_export_trajectory_filename(const std::string &new_prefix) { m_prefix_export_filename = new_prefix; }

    /**
     * @brief If set with true, the simulation is stopped when their is no electron remaining in the device.
     *
     * @param new_value
     */
    void set_stop_simulation_without_electron(bool new_value) { m_simulation_options.m_stop_simu_when_no_electron_remaining = new_value; }

    /**
     * @brief Return true is the avalanche threshold was reached.
     *
     * @return true
     * @return false
     */
    bool has_reached_avalanche() const { return m_list_particles.size() >= m_simulation_options.m_max_number_particle; }

    /**
     * @brief Return a const reference to the simulation history.
     *
     * @return const history_device_admc&
     */
    const history_device_admc &get_simulation_history() const { return m_simulation_history; }

    /**
     * @brief Return a vector with all the current velocities of the particles.
     *
     * @return std::vector<mesh::vector3>
     */
    std::vector<mesh::vector3> get_all_global_velocities() const;

    /**
     * @brief Return a vector with all the number of impact ionization for each particle.
     *
     * @return std::vector<std::size_t>
     */
    std::vector<std::size_t> get_all_number_impact_ionization() const;

    /**
     * @brief Return a vector with all the positions of impact ionizations of the simulation.
     *
     * @return std::vector<std::size_t>
     */
    std::vector<mesh::vector3> get_all_positions_impact_ionization() const;

    /**
     * @brief Return a vector with the impaxct ionization coef of each particle.
     *
     * @return std::vector<double>
     */
    std::vector<double> get_all_coefficient_impact_ionization() const;

    /**
     * @brief Return a vector with the average dead space for each particle.
     *
     * @return std::vector<double>
     */
    std::vector<double> get_all_mean_dead_spaces() const;

    /**
     * @brief Return a vector with the position of the first impact ionzation for each particle.
     * If a particle has no II yet, no value is returned.
     *
     * @return std::vector<mesh::vector3>
     */
    std::vector<mesh::vector3> get_all_first_impact_ionization_position() const;

    /**
     * @brief Get the number electrons.
     *
     * @return std::size_t
     */
    std::size_t get_number_electrons() const;

    /**
     * @brief Get the number holes.
     *
     * @return std::size_t
     */
    std::size_t get_number_holes() const;

    /**
     * @brief Return the time at wich the avalanche threshold was reached.
     *
     * @return std::optional<double>
     */
    std::optional<double> get_time_to_avalanche() const { return m_time_to_avalanche; }

    /**
     * @brief Return the current time of the simulation.
     *
     * @return std::optional<double>
     */
    std::optional<double> get_current_time() const { return m_time; }

    /**
     * @brief Return the vector with all the particles positions.
     *
     * @return std::vector<mesh::vector3>
     */
    std::vector<mesh::vector3> get_all_particles_position() const;

    /**
     * @brief Get the average particles position.
     *
     * @return mesh::vector3
     */
    mesh::vector3 get_average_particles_position() const { return utils::compute_vector_mean(get_all_particles_position()); }

    /**
     * @brief Return the position of the last impact ionization of the simulation.
     *
     * @return mesh::vector3
     */
    mesh::vector3 get_last_impact_ionization_position() const { return m_simulation_history.m_last_impact_ionization_position; }

    /**
     * @brief Return the time at wich the avalanche threshold was reached.
     *
     * @return std::optional<double>
     */
    std::optional<double> get_time_to_secondary_avalanche() const { return m_time_to_secondary_avalanche; }

    /**
     * @brief Compute the depletion region of the device (x_min, x_max).
     * Works only for device with 1D symmetry.
     *
     * @return std::pair<double, double>
     */
    std::pair<double, double> compute_depletion_region() const;

    // /**
    //  * @brief Ad current simlation state to its history.
    //  *
    //  */
    // void add_data_to_history(bool populate_particle_history = false) {
    //     m_simulation_history.add_data_to_history(m_time,
    //                                              get_number_electrons(),
    //                                              get_number_holes(),
    //                                              m_number_impact_ionizations,
    //                                              m_anode_current,
    //                                              m_cathode_current,
    //                                              m_max_global_electric_field);

    //     if (populate_particle_history) {
    //         for (const auto &particle : m_list_particles) {
    //             particle->add_step_to_history();
    //         }
    //     }
    // }

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