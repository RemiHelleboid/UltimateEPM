#include <fmt/format.h>
#include <tclap/CmdLine.h>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "BandStructure.h"
#include "Vector3D.h"
#include "epm_material.hpp"
#include "materials.hpp"

namespace {

struct k_sample {
    std::string      label;
    Vector3D<double> k_reduced;
};

struct band_edge {
    std::string      label;
    Vector3D<double> k_reduced;
    std::size_t      band_index = 0;
    double           energy_eV  = 0.0;
};

struct edge_report {
    band_edge vbm;
    band_edge cbm;
    band_edge gamma_vbm;
    band_edge gamma_cbm;
    band_edge x_cbm;
    band_edge l_cbm;
    band_edge delta_cbm;
    double    indirect_gap_eV = 0.0;
    double    gamma_gap_eV    = 0.0;
};

std::vector<k_sample> make_samples(std::size_t delta_samples) {
    if (delta_samples < 2) {
        throw std::invalid_argument("delta sample count must be at least 2");
    }

    std::vector<k_sample> samples{
        {"Gamma", Vector3D<double>(0.0, 0.0, 0.0)},
        {"X", Vector3D<double>(1.0, 0.0, 0.0)},
        {"L", Vector3D<double>(0.5, 0.5, 0.5)},
    };

    for (std::size_t i = 0; i < delta_samples; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(delta_samples - 1);
        samples.push_back({fmt::format("Delta_{:.8g}", t), Vector3D<double>(t, 0.0, 0.0)});
    }
    return samples;
}

std::size_t find_label(const std::vector<k_sample>& samples, const std::string& label) {
    for (std::size_t i = 0; i < samples.size(); ++i) {
        if (samples[i].label == label) {
            return i;
        }
    }
    throw std::runtime_error("internal error: missing sample label '" + label + "'");
}

band_edge make_edge(const std::vector<k_sample>&              samples,
                    const std::vector<std::vector<double>>& energies,
                    std::size_t                            sample_index,
                    std::size_t                            band_index) {
    return band_edge{samples.at(sample_index).label,
                     samples.at(sample_index).k_reduced,
                     band_index,
                     energies.at(sample_index).at(band_index)};
}

edge_report analyze_edges(const std::vector<k_sample>&              samples,
                          const std::vector<std::vector<double>>& energies,
                          std::size_t                            valence_band,
                          std::size_t                            conduction_band) {
    if (samples.empty() || energies.empty()) {
        throw std::runtime_error("no band samples available");
    }
    if (energies.front().size() <= conduction_band || valence_band >= conduction_band) {
        throw std::invalid_argument("invalid valence/conduction band indices");
    }

    edge_report report;
    report.vbm.energy_eV = std::numeric_limits<double>::lowest();
    report.cbm.energy_eV = std::numeric_limits<double>::infinity();
    report.delta_cbm.energy_eV = std::numeric_limits<double>::infinity();

    for (std::size_t i = 0; i < samples.size(); ++i) {
        const auto& bands = energies.at(i);
        for (std::size_t b = 0; b <= valence_band; ++b) {
            if (bands.at(b) > report.vbm.energy_eV) {
                report.vbm = make_edge(samples, energies, i, b);
            }
        }
        if (bands.at(conduction_band) < report.cbm.energy_eV) {
            report.cbm = make_edge(samples, energies, i, conduction_band);
        }
        if (samples[i].label.rfind("Delta_", 0) == 0 && bands.at(conduction_band) < report.delta_cbm.energy_eV) {
            report.delta_cbm = make_edge(samples, energies, i, conduction_band);
        }
    }

    const std::size_t gamma_index = find_label(samples, "Gamma");
    const std::size_t x_index     = find_label(samples, "X");
    const std::size_t l_index     = find_label(samples, "L");

    report.gamma_vbm = make_edge(samples, energies, gamma_index, valence_band);
    report.gamma_cbm = make_edge(samples, energies, gamma_index, conduction_band);
    report.x_cbm     = make_edge(samples, energies, x_index, conduction_band);
    report.l_cbm     = make_edge(samples, energies, l_index, conduction_band);

    report.indirect_gap_eV = report.cbm.energy_eV - report.vbm.energy_eV;
    report.gamma_gap_eV    = report.gamma_cbm.energy_eV - report.gamma_vbm.energy_eV;
    return report;
}

void print_edge(const std::string& name, const band_edge& edge, double vbm_eV) {
    fmt::print("  {}: band {} at {} k=({:.8g}, {:.8g}, {:.8g}) raw={:.8g} eV rel_vbm={:.8g} eV\n",
               name,
               edge.band_index,
               edge.label,
               edge.k_reduced.X,
               edge.k_reduced.Y,
               edge.k_reduced.Z,
               edge.energy_eV,
               edge.energy_eV - vbm_eV);
}

void print_report(const edge_report& report) {
    fmt::print("Band-edge report\n");
    fmt::print("  indirect gap: {:.8g} eV\n", report.indirect_gap_eV);
    fmt::print("  direct Gamma gap: {:.8g} eV\n", report.gamma_gap_eV);
    print_edge("VBM", report.vbm, report.vbm.energy_eV);
    print_edge("CBM", report.cbm, report.vbm.energy_eV);
    print_edge("Gamma conduction", report.gamma_cbm, report.vbm.energy_eV);
    print_edge("X conduction", report.x_cbm, report.vbm.energy_eV);
    print_edge("L conduction", report.l_cbm, report.vbm.energy_eV);
    print_edge("Delta-line conduction", report.delta_cbm, report.vbm.energy_eV);
}

void write_csv(const std::filesystem::path& path, const edge_report& report) {
    std::ofstream file(path);
    if (!file) {
        throw std::runtime_error("cannot open output file '" + path.string() + "'");
    }

    const auto write_edge = [&](const std::string& name, const band_edge& edge) {
        file << name << "_label," << edge.label << "\n";
        file << name << "_band," << edge.band_index << "\n";
        file << name << "_kx_reduced," << edge.k_reduced.X << "\n";
        file << name << "_ky_reduced," << edge.k_reduced.Y << "\n";
        file << name << "_kz_reduced," << edge.k_reduced.Z << "\n";
        file << name << "_raw_eV," << edge.energy_eV << "\n";
        file << name << "_rel_vbm_eV," << edge.energy_eV - report.vbm.energy_eV << "\n";
    };

    file << "quantity,value\n";
    file << "indirect_gap_eV," << report.indirect_gap_eV << "\n";
    file << "gamma_gap_eV," << report.gamma_gap_eV << "\n";
    write_edge("vbm", report.vbm);
    write_edge("cbm", report.cbm);
    write_edge("gamma_cbm", report.gamma_cbm);
    write_edge("x_cbm", report.x_cbm);
    write_edge("l_cbm", report.l_cbm);
    write_edge("delta_cbm", report.delta_cbm);
}

}  // namespace

int main(int argc, char* argv[]) {
    TCLAP::CmdLine cmd("Report EPM band-edge and gap targets.", ' ', "1.0");

    TCLAP::ValueArg<std::string> arg_material("m", "material", "Material symbol", false, "Si", "string", cmd);
    TCLAP::ValueArg<std::string> arg_epm_set("d", "epm-set", "Named EPM parameter set", false, "local-cohen", "string", cmd);
    TCLAP::ValueArg<int>         arg_nbands("b", "nbands", "Number of bands to compute", false, 8, "int", cmd);
    TCLAP::ValueArg<int> arg_vband("v", "valence-band", "Highest valence band index", false, 3, "int", cmd);
    TCLAP::ValueArg<int> arg_cband("c", "conduction-band", "Lowest conduction band index", false, 4, "int", cmd);
    TCLAP::ValueArg<int> arg_delta_samples("", "delta-samples", "Samples on Gamma-X line", false, 401, "int", cmd);
    TCLAP::ValueArg<int> arg_neighbors("n", "nearestNeighbors", "Number of EPM basis shells", false, 10, "int", cmd);
    TCLAP::ValueArg<int> arg_threads("j", "nthreads", "Number of threads", false, 1, "int", cmd);
    TCLAP::ValueArg<std::string> arg_out("o", "out", "Optional CSV output file", false, "", "path", cmd);
    TCLAP::SwitchArg             arg_nonlocal("C", "nonlocal-correction", "Enable non-local EPM correction", cmd, false);
    TCLAP::SwitchArg             arg_soc("S", "soc", "Enable spin-orbit coupling", cmd, false);

    cmd.parse(argc, argv);

    if (arg_nbands.getValue() <= arg_cband.getValue()) {
        throw TCLAP::ArgException("--nbands must be greater than --conduction-band", arg_nbands.getName());
    }
    if (arg_vband.getValue() < 0 || arg_cband.getValue() <= arg_vband.getValue()) {
        throw TCLAP::ArgException("--valence-band and --conduction-band must define a positive gap", arg_cband.getName());
    }
    if (arg_threads.getValue() <= 0) {
        throw TCLAP::ArgException("number of threads must be positive", arg_threads.getName());
    }

    const uepm::physics::material_repository repository;
    uepm::pseudopotential::Materials         materials;
    materials.load_material(repository, arg_material.getValue(), arg_epm_set.getValue());
    const auto& material = materials.materials.at(arg_material.getValue());

    const auto samples = make_samples(static_cast<std::size_t>(arg_delta_samples.getValue()));
    std::vector<Vector3D<double>> kpoints;
    kpoints.reserve(samples.size());
    for (const auto& sample : samples) {
        kpoints.push_back(sample.k_reduced);
    }

    uepm::pseudopotential::BandStructure band_structure;
    band_structure.Initialize(material,
                              static_cast<std::size_t>(arg_nbands.getValue()),
                              kpoints,
                              static_cast<unsigned int>(arg_neighbors.getValue()),
                              arg_nonlocal.isSet(),
                              arg_soc.isSet());
    band_structure.Compute_parallel(arg_threads.getValue());

    const auto report = analyze_edges(samples,
                                      band_structure.get_band_energies(),
                                      static_cast<std::size_t>(arg_vband.getValue()),
                                      static_cast<std::size_t>(arg_cband.getValue()));
    print_report(report);
    if (!arg_out.getValue().empty()) {
        write_csv(arg_out.getValue(), report);
        fmt::print("Wrote {}\n", arg_out.getValue());
    }

    return 0;
}
