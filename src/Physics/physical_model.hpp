/**
 * @file physical_model.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-11-09
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <algorithm>
#include <fstream>
#include <string>

namespace uepm {

namespace physic {

namespace model {

enum class ModelDataType { analytic, data_file };

enum class ModelType { mobility, impact_ionization, band_gap, ionization_threshold, dielectric_constant, unknown };

class physical_model {
 protected:
    std::string m_model_name;
    ModelType   m_model_type;
    std::string m_model_literature_reference;

 public:
    physical_model() : m_model_name(""), m_model_type(ModelType::unknown), m_model_literature_reference("") {}
    physical_model(const std::string& model_name, ModelType type_of_model) : m_model_name(model_name), m_model_type(type_of_model) {}
    physical_model(const std::string& model_name, ModelType type_of_model, const std::string& model_literature_reference)
        : m_model_name{model_name},
          m_model_type{type_of_model},
          m_model_literature_reference(model_literature_reference) {}
    std::string get_model_name() const { return m_model_name; }
    ModelType   get_model_type() const { return m_model_type; }
    std::string get_literature_reference() const { return m_model_literature_reference; }
};

}  // namespace model

}  // namespace physic

}  // namespace uepm