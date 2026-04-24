/**
 * @file impact_ionization.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-11-08
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once

#include <cmath>
#include <functional>

#include "export_vector_to_csv.hpp"
#include "physical_constants.hpp"
#include "physical_model.hpp"
#include "utils_containers.hpp"

namespace uepm {

namespace physic {
namespace model {

/**
 * @brief Compute the coefficient gamma that represents the temperature dependence of many impact ionization models for silicon.
 *
 *
 * @param temperature
 * @return double
 */
inline static double compute_gamma_temperature_dependence(double temperature) {
    constexpr double ref_temperature = 300.0;
    // reduced planck constant x w_optical_phonon
    constexpr double h_bar_omega = 0.063;
    double           nominator   = tanh(h_bar_omega / (2 * physic::constant::boltzmann_constant_SI * ref_temperature));
    const double     denominator = tanh(h_bar_omega / (2 * physic::constant::boltzmann_constant_SI * temperature));
    const double     gamma       = nominator / denominator;
    return gamma;
}

/**
 * @brief Return the local dead-space depending on electric field and energy threshold for impact ionization.
 *
 * @param electric_field
 * @param impact_ionization_energy_threshold
 * @return double
 */
inline static double compute_dead_space(double electric_field, double impact_ionization_energy_threshold) {
    return (impact_ionization_energy_threshold / electric_field);
    // return (0.0);
}

/**
 * @brief Transform the local coefficient of impact ionization from local model to dead-space model.
 * See 10.1109/TED.2015.2422789 for details.
 *
 * @param local_coeff
 * @return double
 */
inline static double local_to_dead_space_coefficients(double local_coeff, double dead_space) {
    const double dead_space_coefficient = 1.0 / ((1.0 / local_coeff) - 1.0 * dead_space);
    return dead_space_coefficient;
}

/**
 * @brief Impact ionization coefficient alpha for electron according to VanOverstraten and DeMan model in Silicon.
 *
 * @param F_ava
 * @param gamma_temperature_dependence
 * @return double
 */
inline double impact_ionization_rate_electron_van_overstraten_silicon(double F_ava, double gamma_temperature_dependence) {
    constexpr double E_threshold = 0.0;
    if (F_ava <= E_threshold) {
        return 0.0;
    } else {
        constexpr double a_e                  = 7.03e5;
        constexpr double b_e                  = 1.231e6;
        const double     imapact_ionization_e = gamma_temperature_dependence * a_e * exp(-gamma_temperature_dependence * b_e / F_ava);
        return imapact_ionization_e;
    }
}

/**
 * @brief Impact ionization coefficient alpha for hole according to VanOverstraten and DeMan model in Silicon.
 *
 * @param F_ava
 * @param gamma_temperature_dependence
 * @return double
 */
inline double impact_ionization_rate_hole_van_overstraten_silicon(double F_ava, double gamma_temperature_dependence) {
    constexpr double E_threshold = 0.0;
    constexpr double E_0         = 4.0e5;
    if (F_ava <= E_threshold) {
        return 0.0;
    } else if (F_ava <= E_0) {
        constexpr double a_h                  = 1.582e6;
        constexpr double b_h                  = 2.036e6;
        double           imapact_ionization_h = gamma_temperature_dependence * a_h * exp(-gamma_temperature_dependence * b_h / F_ava);
        return imapact_ionization_h;
    } else {
        constexpr double a_h                  = 6.71e5;
        constexpr double b_h                  = 1.693e6;
        double           imapact_ionization_h = gamma_temperature_dependence * a_h * exp(-gamma_temperature_dependence * b_h / F_ava);
        return imapact_ionization_h;
    }
}

/**
 * @brief This class can represent any impact ionization model that take as argument the
 * field strength and the temperature of the lattice.
 *
 */
class impact_ionization_model : public physical_model {
 private:
    std::function<double(double, double)> m_function_impact_ionization_electron_from_field_and_gamma;
    std::function<double(double, double)> m_function_impact_ionization_hole_from_field_and_gamma;
    double                                m_gamma_temperature_dependence;
    bool                                  m_all_argument_constant;
    double                                m_precomputed_impact_ionization_electron = 0.0;
    double                                m_precomputed_impact_ionization_hole     = 0.0;
    bool                                  m_use_dead_space                         = false;
    double                                m_soften_dead_space_factor               = 0.0;
    double                                m_dead_space                             = 0.0;
    double                                m_e_impact_ionization_energy_threshold   = 0.0;
    double                                m_h_impact_ionization_energy_threshold   = 0.0;

 public:
    impact_ionization_model() : m_gamma_temperature_dependence(1.0), m_all_argument_constant(false){};
    // Assignement operator
    // impact_ionization_model& operator=(const impact_ionization_model& other) {
    //     m_function_impact_ionization_electron_from_field_and_gamma = other.m_function_impact_ionization_electron_from_field_and_gamma;
    //     m_function_impact_ionization_hole_from_field_and_gamma     = other.m_function_impact_ionization_hole_from_field_and_gamma;
    //     m_gamma_temperature_dependence                             = other.m_gamma_temperature_dependence;
    //     m_all_argument_constant                                    = other.m_all_argument_constant;
    //     m_precomputed_impact_ionization_electron                   = other.m_precomputed_impact_ionization_electron;
    //     m_precomputed_impact_ionization_hole                       = other.m_precomputed_impact_ionization_hole;
    //     m_use_dead_space                                           = other.m_use_dead_space;
    //     m_soften_dead_space_factor                                 = other.m_soften_dead_space_factor;
    //     m_dead_space                                               = other.m_dead_space;
    //     m_e_impact_ionization_energy_threshold                     = other.m_e_impact_ionization_energy_threshold;
    //     m_h_impact_ionization_energy_threshold                     = other.m_h_impact_ionization_energy_threshold;
    //     return *this;
    // }

    /**
     * @brief Construct a new impact ionization model.
     *
     * @param model_name
     * @param model_literature_reference
     * @param electron_impact_ionization_function
     * @param hole_impact_ionization_function
     * @param lattice_temperature
     */
    impact_ionization_model(const std::string&                          model_name,
                            const std::function<double(double, double)> electron_impact_ionization_function,
                            const std::function<double(double, double)> hole_impact_ionization_function,
                            double                                      lattice_temperature)
        : physical_model(model_name, ModelType::impact_ionization),
          m_function_impact_ionization_electron_from_field_and_gamma{electron_impact_ionization_function},
          m_function_impact_ionization_hole_from_field_and_gamma{hole_impact_ionization_function},
          m_gamma_temperature_dependence{compute_gamma_temperature_dependence(lattice_temperature)},
          m_all_argument_constant{false} {}

    /**
     * @brief Construct a new impact ionization model with constant electric_field.
     *
     * @param model_name
     * @param model_literature_reference
     * @param electron_impact_ionization_function
     * @param hole_impact_ionization_function
     * @param lattice_temperature
     * @param constant_electric_field
     */
    impact_ionization_model(const std::string&                          model_name,
                            const std::function<double(double, double)> electron_impact_ionization_function,
                            const std::function<double(double, double)> hole_impact_ionization_function,
                            double                                      lattice_temperature,
                            double                                      constant_electric_field)
        : physical_model(model_name, ModelType::impact_ionization),
          m_function_impact_ionization_electron_from_field_and_gamma{electron_impact_ionization_function},
          m_function_impact_ionization_hole_from_field_and_gamma{hole_impact_ionization_function},
          m_gamma_temperature_dependence{compute_gamma_temperature_dependence(lattice_temperature)},
          m_all_argument_constant(true),
          m_precomputed_impact_ionization_electron{
              electron_impact_ionization_function(constant_electric_field, m_gamma_temperature_dependence)},
          m_precomputed_impact_ionization_hole{hole_impact_ionization_function(constant_electric_field, m_gamma_temperature_dependence)} {}

    /**
     * @brief Get the electron impact ionization for a given field strength.
     *
     * @param field_strength
     * @return double
     */
    double get_electron_impact_ionization(double field_strength) const {
        double local_e_impact_ionization = 0.0;
        if (m_all_argument_constant) {
            local_e_impact_ionization = m_precomputed_impact_ionization_electron;
        } else {
            local_e_impact_ionization =
                m_function_impact_ionization_electron_from_field_and_gamma(field_strength, m_gamma_temperature_dependence);
        }
        if (m_use_dead_space) {
            const double dead_space = compute_dead_space(field_strength, m_e_impact_ionization_energy_threshold);
            return local_to_dead_space_coefficients(local_e_impact_ionization, dead_space);
        } else {
            return local_e_impact_ionization;
        }
    }

    /**
     * @brief Get the hole impact ionization for a given field strength.
     *
     * @param field_strength
     * @return double
     */
    double get_hole_impact_ionization(double field_strength) const {
        double local_h_impact_ionization = 0.0;
        if (m_all_argument_constant) {
            local_h_impact_ionization = m_precomputed_impact_ionization_hole;
        } else {
            local_h_impact_ionization =
                m_function_impact_ionization_hole_from_field_and_gamma(field_strength, m_gamma_temperature_dependence);
        }
        if (m_use_dead_space) {
            const double dead_space = compute_dead_space(field_strength, m_h_impact_ionization_energy_threshold);
            return local_to_dead_space_coefficients(local_h_impact_ionization, dead_space);
        } else {
            return local_h_impact_ionization;
        }
    }

    void set_use_dead_space(bool use_dead_space) { m_use_dead_space = use_dead_space; }
    bool get_dead_space_activated() const { return m_use_dead_space; }

    double get_dead_space_soft_factor() const { return m_soften_dead_space_factor; }
    void   set_dead_space_soft_factor(double new_value) { m_soften_dead_space_factor = new_value; }

    void make_constant(double constant_electric_field) {
        m_all_argument_constant = true;
        m_precomputed_impact_ionization_electron =
            m_function_impact_ionization_electron_from_field_and_gamma(constant_electric_field, m_gamma_temperature_dependence);
        m_precomputed_impact_ionization_hole =
            m_function_impact_ionization_hole_from_field_and_gamma(constant_electric_field, m_gamma_temperature_dependence);
    }

    void export_impact_ionization_as_csv_function_of_field(const std::string& filename,
                                                           double             min_field,
                                                           double             max_field,
                                                           int                NbPoints) const {
        std::vector<double> electron_impact_ionization(NbPoints);
        std::vector<double> hole_impact_ionization(NbPoints);
        std::vector<double> electric_field_values = utils::linspace(min_field, max_field, NbPoints);
        std::transform(electric_field_values.begin(), electric_field_values.end(), electron_impact_ionization.begin(),
                       [&](auto&& field) { return get_electron_impact_ionization(field); });
        std::transform(electric_field_values.begin(), electric_field_values.end(), hole_impact_ionization.begin(),
                       [&](auto&& field) { return get_hole_impact_ionization(field); });
        utils::export_multiple_vector_to_csv(filename, {"ElectricField", "e_impact_ionization", "h_impact_ionization"},
                                             {electric_field_values, electron_impact_ionization, hole_impact_ionization});
    }
};

static impact_ionization_model VanOverstratenDeManSilicon300K("VanOverstraten&DeMan",
                                                              impact_ionization_rate_electron_van_overstraten_silicon,
                                                              impact_ionization_rate_hole_van_overstraten_silicon,
                                                              300.0);

}  // namespace model

}  // namespace physic

}  // namespace uepm