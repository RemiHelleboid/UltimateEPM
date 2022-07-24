#pragma once

#include <map>
#include <string>

#include "yaml-cpp/yaml.h"


#include "Pseudopotential.h"

namespace EmpiricalPseudopotential {

class Material {
 public:
    Material() : m_lattice_constant(0) {}

    Material(const std::string& Name, double a, double V3S, double V8S, double V11S, double V3A = 0, double V4A = 0, double V11A = 0);

    std::string name;

    double m_lattice_constant;

    Pseudopotential pseudopotential;
};

class Materials {
 public:
    Materials();
    void load_material_parameters(const std::string& filename);

    std::map<std::string, Material> materials;

    void print_materials_list()const;
    void print_material_parameters(const std::string& material_name)const;
    void print_material_parameters()const ;

};

}  // namespace EmpiricalPseudopotential