#include <fmt/format.h>
#include <tclap/CmdLine.h>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "BandStructure.h"
#include "Vector3D.h"
#include "effective_mass.hpp"
#include "epm_material.hpp"
#include "materials.hpp"

namespace {

Vector3D<double> parse_kpoint(const std::string& value) {
    std::stringstream   stream(value);
    std::string         token;
    std::vector<double> components;
    while (std::getline(stream, token, ',')) {
        components.push_back(std::stod(token));
    }
    if (components.size() != 3) {
        throw std::invalid_argument("k-point must have the form kx,ky,kz");
    }
    return Vector3D<double>(components[0], components[1], components[2]);
}

std::vector<Vector3D<double>> make_stencil(const Vector3D<double>& k0, double radius, int shells) {
    if (!(radius > 0.0)) {
        throw std::invalid_argument("radius must be positive");
    }
    if (shells <= 0) {
        throw std::invalid_argument("shell count must be positive");
    }

    std::vector<Vector3D<double>> kpoints;
    kpoints.push_back(k0);
    for (int shell = 1; shell <= shells; ++shell) {
        const double h = radius * static_cast<double>(shell) / static_cast<double>(shells);
        for (int ix = -1; ix <= 1; ++ix) {
            for (int iy = -1; iy <= 1; ++iy) {
                for (int iz = -1; iz <= 1; ++iz) {
                    if (ix == 0 && iy == 0 && iz == 0) {
                        continue;
                    }
                    kpoints.emplace_back(k0.X + h * static_cast<double>(ix),
                                         k0.Y + h * static_cast<double>(iy),
                                         k0.Z + h * static_cast<double>(iz));
                }
            }
        }
    }
    return kpoints;
}

std::vector<uepm::pseudopotential::valley_fit_sample> compute_samples(
    const uepm::pseudopotential::epm_material& material,
    const std::vector<Vector3D<double>>&       kpoints,
    int                                        band_index,
    int                                        nearest_neighbors,
    bool                                       enable_nonlocal,
    bool                                       enable_soc,
    int                                        nb_threads) {
    uepm::pseudopotential::BandStructure band_structure;
    band_structure.Initialize(material,
                              static_cast<std::size_t>(band_index + 1),
                              kpoints,
                              static_cast<unsigned int>(nearest_neighbors),
                              enable_nonlocal,
                              enable_soc);
    band_structure.Compute_parallel(nb_threads);

    std::vector<uepm::pseudopotential::valley_fit_sample> samples;
    samples.reserve(kpoints.size());
    const auto& energies = band_structure.get_band_energies();
    for (std::size_t i = 0; i < kpoints.size(); ++i) {
        samples.push_back({kpoints[i], energies.at(i).at(static_cast<std::size_t>(band_index))});
    }
    return samples;
}

void write_csv(const std::filesystem::path& path, const uepm::pseudopotential::effective_mass_fit_result& result) {
    std::ofstream file(path);
    if (!file) {
        throw std::runtime_error("cannot open output file '" + path.string() + "'");
    }

    file << "quantity,value\n";
    file << "edge_kind," << uepm::pseudopotential::to_string(result.edge_kind) << "\n";
    file << "k0_x_reduced," << result.k0_reduced.X << "\n";
    file << "k0_y_reduced," << result.k0_reduced.Y << "\n";
    file << "k0_z_reduced," << result.k0_reduced.Z << "\n";
    file << "edge_energy_eV," << result.edge_energy_eV << "\n";
    file << "m1_m0," << result.principal_masses_m0[0] << "\n";
    file << "m2_m0," << result.principal_masses_m0[1] << "\n";
    file << "m3_m0," << result.principal_masses_m0[2] << "\n";
    file << "non_parabolicity_eV_inv," << result.non_parabolicity_eV_inv << "\n";
    file << "mass_rms_error_meV," << result.mass_rms_error_meV << "\n";
    file << "alpha_rms_error_meV," << result.alpha_rms_error_meV << "\n";
    file << "sample_count," << result.sample_count << "\n";
    file << "mass_sample_count," << result.mass_sample_count << "\n";
    file << "alpha_sample_count," << result.alpha_sample_count << "\n";
}

void print_result(const uepm::pseudopotential::effective_mass_fit_result& result) {
    fmt::print("Valley fit ({})\n", uepm::pseudopotential::to_string(result.edge_kind));
    fmt::print("  k0 reduced: ({:.8g}, {:.8g}, {:.8g})\n",
               result.k0_reduced.X,
               result.k0_reduced.Y,
               result.k0_reduced.Z);
    fmt::print("  edge energy: {:.8g} eV\n", result.edge_energy_eV);
    fmt::print("  principal masses: {:.8g}, {:.8g}, {:.8g} m0\n",
               result.principal_masses_m0[0],
               result.principal_masses_m0[1],
               result.principal_masses_m0[2]);
    fmt::print("  non-parabolicity: {:.8g} 1/eV\n", result.non_parabolicity_eV_inv);
    fmt::print("  RMS mass residual: {:.8g} meV over {} mass samples\n",
               result.mass_rms_error_meV,
               result.mass_sample_count);
    fmt::print("  RMS fixed-mass Kane residual: {:.8g} meV over {} alpha samples\n",
               result.alpha_rms_error_meV,
               result.alpha_sample_count);
    if (result.alpha_rms_error_meV > 5.0) {
        fmt::print("  Warning: alpha fit residual is large; shrink --radius or --alpha-max-energy.\n");
    }
    fmt::print("  principal axes are columns:\n");
    for (std::size_t row = 0; row < 3; ++row) {
        fmt::print("    [{: .8g}, {: .8g}, {: .8g}]\n",
                   result.principal_axes[row][0],
                   result.principal_axes[row][1],
                   result.principal_axes[row][2]);
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    TCLAP::CmdLine cmd("Fit effective masses and Kane non-parabolicity from local EPM samples.", ' ', "1.0");

    TCLAP::ValueArg<std::string> arg_material("m", "material", "Material symbol", false, "Si", "string", cmd);
    TCLAP::ValueArg<std::string>
                         arg_epm_set("d", "epm-set", "Named EPM parameter set", false, "local-cohen", "string", cmd);
    TCLAP::ValueArg<int> arg_band("b", "band", "Zero-based band index to fit", false, 4, "int", cmd);
    TCLAP::ValueArg<std::string> arg_k0("", "k0", "Reduced valley center kx,ky,kz", false, "0.85,0,0", "string", cmd);
    TCLAP::ValueArg<std::string> arg_edge("e", "edge", "Band edge kind: min or max", false, "min", "string", cmd);
    TCLAP::ValueArg<double>
        arg_radius("", "radius", "Reduced-coordinate alpha sampling radius", false, 0.04, "float", cmd);
    TCLAP::ValueArg<int> arg_shells("", "shells", "Number of alpha stencil shells", false, 4, "int", cmd);
    TCLAP::ValueArg<double>
        arg_mass_radius("", "mass-radius", "Near-edge mass sampling radius", false, 0.01, "float", cmd);
    TCLAP::ValueArg<int>    arg_mass_shells("", "mass-shells", "Number of mass stencil shells", false, 2, "int", cmd);
    TCLAP::ValueArg<double> arg_alpha_max_energy("",
                                                 "alpha-max-energy",
                                                 "Maximum kinetic energy used for alpha fit",
                                                 false,
                                                 0.3,
                                                 "float",
                                                 cmd);
    TCLAP::ValueArg<int>    arg_neighbors("n", "nearestNeighbors", "Number of EPM basis shells", false, 10, "int", cmd);
    TCLAP::ValueArg<int>    arg_threads("j", "nthreads", "Number of threads", false, 1, "int", cmd);
    TCLAP::ValueArg<std::string> arg_out("o", "out", "Optional CSV output file", false, "", "path", cmd);
    TCLAP::SwitchArg arg_nonlocal("C", "nonlocal-correction", "Enable non-local EPM correction", cmd, false);
    TCLAP::SwitchArg arg_soc("S", "soc", "Enable spin-orbit coupling", cmd, false);
    TCLAP::SwitchArg arg_allow_negative_alpha("",
                                              "allow-negative-alpha",
                                              "Do not clamp negative alpha to zero",
                                              cmd,
                                              false);

    cmd.parse(argc, argv);

    if (arg_band.getValue() < 0) {
        throw TCLAP::ArgException("band index must be non-negative", arg_band.getName());
    }
    if (arg_threads.getValue() <= 0) {
        throw TCLAP::ArgException("number of threads must be positive", arg_threads.getName());
    }

    const Vector3D<double> k0            = parse_kpoint(arg_k0.getValue());
    const auto             edge_kind     = uepm::pseudopotential::band_edge_kind_from_string(arg_edge.getValue());
    const auto             mass_kpoints  = make_stencil(k0, arg_mass_radius.getValue(), arg_mass_shells.getValue());
    const auto             alpha_kpoints = make_stencil(k0, arg_radius.getValue(), arg_shells.getValue());

    const uepm::physics::material_repository repository;
    uepm::pseudopotential::Materials         materials;
    materials.load_material(repository, arg_material.getValue(), arg_epm_set.getValue());
    const auto& material = materials.materials.at(arg_material.getValue());

    const auto mass_samples  = compute_samples(material,
                                              mass_kpoints,
                                              arg_band.getValue(),
                                              arg_neighbors.getValue(),
                                              arg_nonlocal.isSet(),
                                              arg_soc.isSet(),
                                              arg_threads.getValue());
    const auto alpha_samples = compute_samples(material,
                                               alpha_kpoints,
                                               arg_band.getValue(),
                                               arg_neighbors.getValue(),
                                               arg_nonlocal.isSet(),
                                               arg_soc.isSet(),
                                               arg_threads.getValue());

    const double edge_energy_eV = mass_samples.front().energy_eV;
    const auto   result =
        uepm::pseudopotential::fit_effective_mass_then_nonparabolicity(mass_samples,
                                                                       alpha_samples,
                                                                       k0,
                                                                       edge_energy_eV,
                                                                       material.get_lattice_constant_meter(),
                                                                       edge_kind,
                                                                       arg_alpha_max_energy.getValue(),
                                                                       !arg_allow_negative_alpha.isSet());

    print_result(result);
    if (!arg_out.getValue().empty()) {
        write_csv(arg_out.getValue(), result);
        fmt::print("Wrote {}\n", arg_out.getValue());
    }

    return 0;
}
