/**
 * @file mobility_model.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-11-09
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

#include "export_vector_to_csv.hpp"
#include "physical_model.hpp"
#include "numerical_helper.hpp"
#include "utils_containers.hpp"

namespace uepm {

namespace physic {

namespace model {

// Constants for Arora model
constexpr double REF_TEMPERATURE = 300.0;

// Electron mobility constants
constexpr double A_min_electron     = 88.0;
constexpr double alpha_min_electron = -0.57;
constexpr double A_dop_electron     = 1252.0;
constexpr double alpha_dop_electron = -2.33;
constexpr double A_N_electron       = 1.25e17;
constexpr double alpha_N_electron   = 2.4;
constexpr double A_a_electron       = 0.88;
constexpr double alpha_a_electron   = -0.146;

// Hole mobility constants
constexpr double A_min_hole     = 54.3;
constexpr double alpha_min_hole = -0.57;
constexpr double A_dop_hole     = 407.0;
constexpr double alpha_dop_hole = -2.23;
constexpr double A_N_hole       = 2.35e17;
constexpr double alpha_N_hole   = 2.4;
constexpr double A_a_hole       = 0.88;
constexpr double alpha_a_hole   = -0.146;

// Saturation velocity constants
constexpr double v_sat_min_electron      = 1.07e7;
constexpr double v_sat_exponent_electron = 0.87;
constexpr double v_sat_min_hole          = 8.37e6;
constexpr double v_sat_exponent_hole     = 0.52;

// Canali model constants
constexpr double beta_0_electron   = 1.109;
constexpr double beta_exp_electron = 0.66;
constexpr double beta_0_hole       = 1.213;
constexpr double beta_exp_hole     = 0.17;
constexpr double alpha             = 0.0;

/**
 * @brief Compute the electron mobility using Arora model for Silicon.
 *
 * @param doping_concentration
 * @param Temperature
 * @return double
 */
inline double electron_mobility_Arora(double doping_concentration, double Temperature) {
    const double relative_temperature = Temperature / REF_TEMPERATURE;
    const double mu_min               = A_min_electron * std::pow(relative_temperature, alpha_min_electron);
    const double mu_d                 = A_dop_electron * std::pow(relative_temperature, alpha_dop_electron);
    const double N_0                  = A_N_electron * std::pow(relative_temperature, alpha_N_electron);
    const double A_star               = A_a_electron * std::pow(relative_temperature, alpha_a_electron);
    const double doping_coefficient   = std::pow((std::fabs(doping_concentration) / N_0), A_star);
    const double mu_dop               = mu_d / (1.0 + doping_coefficient);
    return mu_min + mu_dop;
}

/**
 * @brief Compute the electron mobility using Arora model for Silicon at 300K.
 * Exists to fasten the computation, by avoiding the temperature dependency.
 *
 * @param doping_concentration
 * @param Temperature
 * @return double
 */
inline double electron_mobility_Arora_300K(double doping_concentration) {
    const double relative_temperature = 1.0;
    const double mu_min               = A_min_electron;
    const double mu_d                 = A_dop_electron;
    const double N_0                  = A_N_electron;
    const double A_star               = A_a_electron;
    const double doping_coefficient   = std::pow((std::fabs(doping_concentration) / N_0), A_star);
    const double mu_dop               = mu_d / (1.0 + doping_coefficient);
    return mu_min + mu_dop;
}

/**
 * @brief Compute the hole mobility using Arora model for Silicon.
 *
 * @param doping_concentration
 * @param Temperature
 * @return double
 */
inline double hole_mobility_Arora(double doping_concentration, double Temperature) {
    const double relative_temperature = Temperature / REF_TEMPERATURE;
    const double mu_min               = A_min_hole * std::pow(relative_temperature, alpha_min_hole);
    const double mu_d                 = A_dop_hole * std::pow(relative_temperature, alpha_dop_hole);
    const double N_0                  = A_N_hole * std::pow(relative_temperature, alpha_N_hole);
    const double A_star               = A_a_hole * std::pow(relative_temperature, alpha_a_hole);
    const double doping_coefficient   = std::pow((std::fabs(doping_concentration) / N_0), A_star);
    const double mu_dop               = mu_d / (1.0 + doping_coefficient);
    return mu_min + mu_dop;
}

/**
 * @brief Compute the hole mobility using Arora model for Silicon at 300K.
 * Exists to fasten the computation, by avoiding the temperature dependency.
 *
 * @param doping_concentration
 * @param Temperature
 * @return double
 */
inline double hole_mobility_Arora_300K(double doping_concentration) {
    const double relative_temperature = 1.0;
    const double mu_min               = A_min_hole;
    const double mu_d                 = A_dop_hole;
    const double N_0                  = A_N_hole;
    const double A_star               = A_a_hole;
    const double doping_coefficient   = std::pow((std::fabs(doping_concentration) / N_0), A_star);
    const double mu_dop               = mu_d / (1.0 + doping_coefficient);
    return mu_min + mu_dop;
}

/**
 * @brief Compute the electron saturation velocity.
 *
 * @param temperature
 * @return double
 */
inline double electron_saturation_velocity(double temperature) {
    return v_sat_min_electron * std::pow((temperature / REF_TEMPERATURE), v_sat_exponent_electron);
}

/**
 * @brief Compute the electron saturation velocity at 300K.
 *  Exists to fasten the computation, by avoiding the temperature dependency.
 * @return double
 */
inline double electron_saturation_velocity_300K() { return v_sat_min_electron; }

/**
 * @brief Compute the hole saturation velocity.
 *
 * @param temperature
 * @return double
 */
inline double hole_saturation_velocity(double temperature) {
    return v_sat_min_hole * std::pow((temperature / REF_TEMPERATURE), v_sat_exponent_hole);
}

/**
 * @brief Compute the hole saturation velocity at 300K.
 *  Exists to fasten the computation, by avoiding the temperature dependency.
 * @return double
 */
inline double hole_saturation_velocity_300K() { return v_sat_min_hole; }

/**
 * @brief Mobility for electron using Canali model for high fields for Silicon.
 *
 * @param doping_concentration
 * @param ElectricField
 * @param Temperature
 * @return double
 */
inline double electron_mobility_Arora_Canali(double doping_concentration, double ElectricField, double Temperature) {
    const double beta               = beta_0_electron * std::pow((Temperature / REF_TEMPERATURE), beta_exp_electron);
    const double mobility_low_field = electron_mobility_Arora(doping_concentration, Temperature);
    const double v_sat              = electron_saturation_velocity(Temperature);
    const double numerator          = (alpha + 1.0) * mobility_low_field;
    const double denominator =
        alpha + std::pow(1.0 + std::pow((((alpha + 1.0) * mobility_low_field * ElectricField) / v_sat), beta), (1.0 / beta));
    return numerator / denominator;
}

/**
 * @brief Mobility for electron using Canali model for high fields for Silicon at 300K.
 * Exists to fasten the computation, by avoiding the temperature dependency.
 *
 * @param doping_concentration
 * @param ElectricField
 * @return double
 */
inline double electron_mobility_Arora_Canali_300K(double doping_concentration, double ElectricField, double T_null) {
    const double beta               = beta_0_electron;
    const double mobility_low_field = electron_mobility_Arora_300K(doping_concentration);
    const double v_sat              = electron_saturation_velocity_300K();
    const double numerator          = (alpha + 1.0) * mobility_low_field;
    const double denominator =
        alpha + std::pow(1.0 + std::pow((((alpha + 1.0) * mobility_low_field * ElectricField) / v_sat), beta), (1.0 / beta));
    return numerator / denominator;
}

/**
 * @brief Mobility for hole using Canali model for high fields for Silicon.
 *
 * @param doping_concentration
 * @param ElectricField
 * @param Temperature
 * @return double
 */
inline double hole_mobility_Arora_Canali(double doping_concentration, double ElectricField, double Temperature) {
    const double beta               = beta_0_hole * std::pow((Temperature / REF_TEMPERATURE), beta_exp_hole);
    const double mobility_low_field = hole_mobility_Arora(doping_concentration, Temperature);
    const double v_sat              = hole_saturation_velocity(Temperature);
    const double numerator          = (alpha + 1.0) * mobility_low_field;
    const double denominator =
        alpha + std::pow(1.0 + std::pow((((alpha + 1.0) * mobility_low_field * ElectricField) / v_sat), beta), (1.0 / beta));
    return numerator / denominator;
}

/**
 * @brief Mobility for hole using Canali model for high fields for Silicon at 300K.
 * Exists to fasten the computation, by avoiding the temperature dependency.
 *
 * @param doping_concentration
 * @param ElectricField
 * @return double
 */
inline double hole_mobility_Arora_Canali_300K(double doping_concentration, double ElectricField, double T_null) {
    const double beta               = beta_0_hole;
    const double mobility_low_field = hole_mobility_Arora_300K(doping_concentration);
    const double v_sat              = hole_saturation_velocity_300K();
    const double numerator          = (alpha + 1.0) * mobility_low_field;
    const double denominator =
        alpha + std::pow(1.0 + std::pow((((alpha + 1.0) * mobility_low_field * ElectricField) / v_sat), beta), (1.0 / beta));
    return numerator / denominator;
}

/**
 * @brief Mobility model class.
 *
 */
class mobility_model : public physical_model {
 private:
    std::function<double(double, double, double)> m_electron_mobility_temperature_doping_electric_field;
    std::function<double(double, double, double)> m_hole_mobility_temperature_doping_electric_field;
    bool                                          m_all_argument_constant   = false;
    double                                        m_pre_computed_e_mobility = 0.0;
    double                                        m_pre_computed_h_mobility = 0.0;

 public:
    mobility_model() = default;

    mobility_model(const std::string                            &model_name,
                   std::function<double(double, double, double)> electron_mobility_function,
                   std::function<double(double, double, double)> hole_mobility_function)
        : physical_model(model_name, ModelType::mobility),
          m_electron_mobility_temperature_doping_electric_field{std::move(electron_mobility_function)},
          m_hole_mobility_temperature_doping_electric_field{std::move(hole_mobility_function)} {}

    mobility_model(const std::string                            &model_name,
                   std::function<double(double, double, double)> electron_mobility_function,
                   std::function<double(double, double, double)> hole_mobility_function,
                   double                                        constant_electric_field,
                   double                                        constant_doping,
                   double                                        constant_temperature)
        : physical_model(model_name, ModelType::mobility),
          m_electron_mobility_temperature_doping_electric_field{std::move(electron_mobility_function)},
          m_hole_mobility_temperature_doping_electric_field{std::move(hole_mobility_function)},
          m_all_argument_constant(true),
          m_pre_computed_e_mobility{
              m_electron_mobility_temperature_doping_electric_field(constant_doping, constant_electric_field, constant_temperature)},
          m_pre_computed_h_mobility{
              m_hole_mobility_temperature_doping_electric_field(constant_doping, constant_electric_field, constant_temperature)} {}

    double get_electron_mobility(double temperature, double doping_concentration, double electric_field) const {
        return m_electron_mobility_temperature_doping_electric_field(doping_concentration, electric_field, temperature);
    }

    double get_constant_electron_mobility() const {
        if (!m_all_argument_constant) {
            throw std::runtime_error("The mobility model is not constant, you cannot call get_constant_electron_mobility");
        }
        return m_pre_computed_e_mobility;
    }

    double get_hole_mobility(double temperature, double doping_concentration, double electric_field) const {
        return m_hole_mobility_temperature_doping_electric_field(doping_concentration, electric_field, temperature);
    }

    double get_constant_hole_mobility() const {
        if (!m_all_argument_constant) {
            throw std::runtime_error("The mobility model is not constant, you cannot call get_constant_hole_mobility");
        }
        return m_pre_computed_h_mobility;
    }

    void make_constant(double constant_temperature, double constant_doping, double constant_electric_field) {
        m_all_argument_constant = true;
        m_pre_computed_e_mobility =
            m_electron_mobility_temperature_doping_electric_field(constant_doping, constant_electric_field, constant_temperature);
        m_pre_computed_h_mobility =
            m_hole_mobility_temperature_doping_electric_field(constant_doping, constant_electric_field, constant_temperature);
    }

    void export_velocity_as_csv_function_of_field(const std::string &filename,
                                                  double             min_field,
                                                  double             max_field,
                                                  int                NbPoints,
                                                  double             doping,
                                                  double             temperature) const {
        std::vector<double> e_velocity;
        std::vector<double> h_velocity;
        std::vector<double> electric_field_values = utils::linspace(min_field, max_field, NbPoints);

        e_velocity.reserve(NbPoints);
        h_velocity.reserve(NbPoints);

        for (const auto &field : electric_field_values) {
            e_velocity.push_back(field * get_electron_mobility(temperature, doping, field));
            h_velocity.push_back(field * get_hole_mobility(temperature, doping, field));
        }

        utils::export_multiple_vector_to_csv(filename, {"ElectricField", "e_Velocity", "h_Velocity"},
                                             {electric_field_values, e_velocity, h_velocity});
    }
};

// Predefined mobility model for Silicon
static mobility_model Arora_Canali_Mobility_Silicon("Arora&Canali_Silicon", electron_mobility_Arora_Canali, hole_mobility_Arora_Canali);
static mobility_model Arora_Canali_Mobility_Silicon300K("Arora&Canali_Silicon",
                                                        electron_mobility_Arora_Canali_300K,
                                                        hole_mobility_Arora_Canali_300K);

}  // namespace model

}  // namespace physic

} // namespace uepm