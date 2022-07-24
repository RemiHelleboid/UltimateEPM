#include "Material.h"

#include "yaml-cpp/yaml.h"

namespace EmpiricalPseudopotential {

const double Bohr = 0.52917721092;  // in Angstroms

Material::Material(const std::string& Name, double a, double V3S, double V8S, double V11S, double V3A, double V4A, double V11A)
    : name(Name),
      m_lattice_constant(a / Bohr),  // convert it to Bohrs, the value is given in Angstroms
      pseudopotential(V3S / 2.,
                      V8S / 2.,
                      V11S / 2.,
                      V3A / 2.,
                      V4A / 2.,
                      V11A / 2.)  // the values are in Rydbergs, convert them to Hartree, one Rydberg is half a Hartree
{}

Materials::Materials() {
    // values are given in Angstroms and Rydbergs, the constructor converts them to atomic units, see above
    // materials["Si"]   = Material("Si", 5.43, -0.2241, 0.052, 0.0724);
    // materials["Ge"]   = Material("Ge", 5.66, -0.23, 0.01, 0.06);
    // materials["Sn"]   = Material("Sn", 6.49, -0.20, 0.00, 0.04);
    // materials["GaP"]  = Material("GaP", 5.44, -0.22, 0.03, 0.07, 0.12, 0.07, 0.02);
    // materials["GaAs"] = Material("GaAs", 5.64, -0.23, 0.01, 0.06, 0.07, 0.05, 0.01);
    // materials["AlSb"] = Material("AlSb", 6.13, -0.21, 0.02, 0.06, 0.06, 0.04, 0.02);
    // materials["InP"]  = Material("InP", 5.86, -0.23, 0.01, 0.06, 0.07, 0.05, 0.01);
    // materials["GaSb"] = Material("GaSb", 6.12, -0.22, 0.00, 0.05, 0.06, 0.05, 0.01);
    // materials["InAs"] = Material("InAs", 6.04, -0.22, 0.00, 0.05, 0.08, 0.05, 0.03);
    // materials["InSb"] = Material("InSb", 6.48, -0.20, 0.00, 0.04, 0.06, 0.05, 0.01);
    // materials["ZnS"]  = Material("ZnS", 5.41, -0.22, 0.03, 0.07, 0.24, 0.14, 0.04);
    // materials["ZnSe"] = Material("ZnSe", 5.65, -0.23, 0.01, 0.06, 0.18, 0.12, 0.03);
    // materials["ZnTe"] = Material("ZnTe", 6.07, -0.22, 0.00, 0.05, 0.13, 0.10, 0.01);
    // materials["CdTe"] = Material("CdTe", 6.41, -0.20, 0.00, 0.04, 0.15, 0.09, 0.04);

    // // SiGe alloy, x = 0.5
    // materials["SiGe"] = Material("SiGe", 5.495, -0.22, 0.025, 0.07, 0., 0., 0.);

    // // Added AlAs from here: https://www.ece.nus.edu.sg/stfpage/eleadj/pseudopotential.htm
    // // lattice constant from here: https://sector7.xray.aps.anl.gov/calculators/crystal_lattice_parameters.html
    // materials["AlAs"] = Material("AlAs", 5.6605, -0.221, 0.025, 0.07, 0.08, 0.05, -0.004);
}

void Materials::load_material_parameters(const std::string& filename) {
    YAML::Node config         = YAML::LoadFile(filename);
    auto       list_materials = config["materials"];
    for (const auto& material : list_materials) {
        std::string name   = material["name"].as<std::string>();
        std::string symbol = material["symbol"].as<std::string>();
        double      a      = material["lattice-constant"].as<double>();
        double      V3S    = material["V3S"].as<double>();
        double      V8S    = material["V8S"].as<double>();
        double      V11S   = material["V11S"].as<double>();
        double      V3A    = material["V3A"].as<double>();
        double      V4A    = material["V4A"].as<double>();
        double      V11A   = material["V11A"].as<double>();
        materials[name]    = Material(symbol, a, V3S, V8S, V11S, V3A, V4A, V11A);
    }
}

void Materials::print_materials_list() const {
    for (const auto& material : materials) {
        std::cout << material.first << std::endl;
    }
}

void Materials::print_material_parameters(const std::string& name) const {
    if (materials.find(name) == materials.end()) {
        std::cout << "Material " << name << " not found" << std::endl;
        return;
    }
    const Material& material = materials.at(name);
    std::cout << "Material: " << name << std::endl;
    std::cout << "Lattice constant: " << material.m_lattice_constant << " Bohr" << std::endl;
    std::cout << "Empirical Pseudopotential Parameters:" << std::endl;
    material.pseudopotential.print_parameters();
    std::cout << "-------------------------------------" << std::endl;
}

void Materials::print_material_parameters() const {
    for (const auto& material : materials) {
        print_material_parameters(material.first);
    }
}

}  // namespace EmpiricalPseudopotential
