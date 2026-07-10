/**
 * @file epsilon.cpp
 * @brief MPI dielectric-function validation and generation app.
 */

#include <mpi.h>
#include <tclap/CmdLine.h>
#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "BandStructure.h"
#include "DielectricFunction.hpp"
#include "Options.h"

namespace {

struct EpsilonAppConfig {
    std::string      material_name{"Si"};
    std::string      epm_parameter_set{"local-cohen"};
    std::string      output_prefix{"epsilon"};
    std::string      mode{"q-list"};
    std::string      q_file{};
    std::string      q_values{"1e-2,1e-3,1e-4"};
    Vector3D<double> direction{1.0, 0.0, 0.0};

    int         nb_nearest_neighbors{10};
    int         nb_bands{16};
    int         nkx{40};
    int         nky{40};
    int         nkz{40};
    std::string bz_sampling{"full"};
    int         q_count{40};

    double min_energy_eV{0.0};
    double max_energy_eV{20.0};
    double energy_step_eV{0.01};
    double eta_smearing_eV{0.05};
    double q_min{1e-2};
    double q_max{3.0};
    double small_q_diagnostic_threshold{5e-2};

    bool nonlocal_epm{false};
    bool enable_soc{false};
};

template <typename T>
void load_yaml_value(const YAML::Node& node, const char* key, T& value) {
    if (node[key]) {
        value = node[key].as<T>();
    }
}

void load_yaml_value_any_key(const YAML::Node& node, const std::vector<const char*>& keys, int& value) {
    for (const char* key : keys) {
        if (node[key]) {
            value = node[key].as<int>();
            return;
        }
    }
}

std::string lowercase_ascii(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return text;
}

std::string normalize_bz_sampling(std::string value) {
    value = lowercase_ascii(value);
    std::replace(value.begin(), value.end(), '_', '-');
    if (value == "1" || value == "full" || value == "full-bz") {
        return "full";
    }
    if (value == "8" || value == "q100-octant" || value == "100-octant" || value == "octant") {
        return "q100-octant";
    }
    if (value == "48" || value == "fcc-ibz" || value == "ibz" || value == "irreducible-wedge") {
        return "fcc-ibz";
    }
    throw std::invalid_argument("--bz-sampling must be one of full, q100-octant, fcc-ibz (legacy aliases: 1, 8, 48).");
}

uepm::pseudopotential::DielectricKPointSampling dielectric_sampling_mode(const std::string& value) {
    const std::string mode = normalize_bz_sampling(value);
    if (mode == "full") {
        return uepm::pseudopotential::DielectricKPointSampling::full_bz;
    }
    if (mode == "q100-octant") {
        return uepm::pseudopotential::DielectricKPointSampling::q100_octant;
    }
    return uepm::pseudopotential::DielectricKPointSampling::fcc_irreducible_wedge;
}

bool vector_along_100(const Vector3D<double>& vector) {
    constexpr double tolerance = 1.0e-12;
    return std::abs(vector.X) > tolerance && std::abs(vector.Y) <= tolerance && std::abs(vector.Z) <= tolerance;
}

bool all_qpoints_along_100(const std::vector<Vector3D<double>>& qpoints) {
    return std::all_of(qpoints.begin(), qpoints.end(), vector_along_100);
}

std::vector<double> parse_double_list(std::string text) {
    std::replace(text.begin(), text.end(), ':', ',');
    std::replace(text.begin(), text.end(), ';', ',');

    std::vector<double> values;
    std::stringstream   stream(text);
    std::string         token;
    while (std::getline(stream, token, ',')) {
        if (token.empty()) {
            continue;
        }
        values.push_back(std::stod(token));
    }
    return values;
}

Vector3D<double> parse_vector3(std::string text) {
    std::replace(text.begin(), text.end(), ';', ',');
    std::replace(text.begin(), text.end(), ' ', ',');
    const auto values = parse_double_list(text);
    if (values.size() != 3) {
        throw std::invalid_argument("Expected a vector with three components, e.g. 1,0,0.");
    }
    return Vector3D<double>(values[0], values[1], values[2]);
}

Vector3D<double> normalized_direction(Vector3D<double> direction) {
    const double norm = direction.Length();
    if (!(norm > 0.0)) {
        throw std::invalid_argument("Direction vector must be non-zero.");
    }
    return direction / norm;
}

std::vector<double> make_energy_grid(double emin, double emax, double estep) {
    if (!(estep > 0.0)) {
        throw std::invalid_argument("Energy step must be positive.");
    }
    if (emax < emin) {
        throw std::invalid_argument("Maximum energy must be larger than minimum energy.");
    }

    std::vector<double> energies;
    for (double energy = emin; energy <= emax + 0.5 * estep; energy += estep) {
        energies.push_back(energy);
    }
    return energies;
}

std::vector<Vector3D<double>> read_qpoint_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open q-point file: " + filename);
    }

    std::vector<Vector3D<double>> qpoints;
    std::string                   line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::replace(line.begin(), line.end(), ',', ' ');
        std::stringstream stream(line);
        double            qx = 0.0;
        double            qy = 0.0;
        double            qz = 0.0;
        if (stream >> qx >> qy >> qz) {
            qpoints.emplace_back(qx, qy, qz);
        }
    }
    return qpoints;
}

std::vector<Vector3D<double>> make_qpoints(const EpsilonAppConfig& config) {
    if (config.mode == "optical") {
        return {Vector3D<double>(0.0, 0.0, 0.0)};
    }

    if (!config.q_file.empty()) {
        auto qpoints = read_qpoint_file(config.q_file);
        if (qpoints.empty()) {
            throw std::runtime_error("The q-point file did not contain any valid q vectors.");
        }
        return qpoints;
    }

    const Vector3D<double> direction = normalized_direction(config.direction);
    std::vector<double>    q_norms;
    if (config.mode == "q-list") {
        q_norms = parse_double_list(config.q_values);
    } else if (config.mode == "q-line") {
        if (config.q_count < 2) {
            throw std::invalid_argument("q-line mode requires --q-count >= 2.");
        }
        if (config.q_max < config.q_min) {
            throw std::invalid_argument("q-line mode requires --q-max >= --q-min.");
        }
        const double step = (config.q_max - config.q_min) / static_cast<double>(config.q_count - 1);
        for (int i = 0; i < config.q_count; ++i) {
            q_norms.push_back(config.q_min + static_cast<double>(i) * step);
        }
    } else if (config.mode == "q-file") {
        throw std::invalid_argument("q-file mode requires --q-file or YAML file-list-q.");
    } else {
        throw std::invalid_argument("Unknown epsilon mode '" + config.mode +
                                    "'. Use optical, q-list, q-line, or q-file.");
    }

    if (q_norms.empty()) {
        throw std::invalid_argument("No q values were provided.");
    }

    std::vector<Vector3D<double>> qpoints;
    qpoints.reserve(q_norms.size());
    for (double q_norm : q_norms) {
        qpoints.push_back(direction * q_norm);
    }
    return qpoints;
}

void load_config_file(const std::string& filename, EpsilonAppConfig& config) {
    if (filename.empty()) {
        return;
    }

    const YAML::Node yaml = YAML::LoadFile(filename);
    load_yaml_value(yaml, "material", config.material_name);
    load_yaml_value(yaml, "epm-set", config.epm_parameter_set);
    load_yaml_value(yaml, "parameter-set", config.epm_parameter_set);
    load_yaml_value(yaml, "nb-bands", config.nb_bands);
    load_yaml_value_any_key(yaml, {"nearest-neighbors", "nearest-neigbors"}, config.nb_nearest_neighbors);
    load_yaml_value(yaml, "nonlocal", config.nonlocal_epm);
    load_yaml_value(yaml, "min-energy", config.min_energy_eV);
    load_yaml_value(yaml, "max-energy", config.max_energy_eV);
    load_yaml_value(yaml, "step-energy", config.energy_step_eV);
    load_yaml_value(yaml, "eta-smearing", config.eta_smearing_eV);
    load_yaml_value(yaml, "Nkx", config.nkx);
    load_yaml_value(yaml, "Nky", config.nky);
    load_yaml_value(yaml, "Nkz", config.nkz);
    load_yaml_value(yaml, "bz-sampling", config.bz_sampling);
    load_yaml_value(yaml, "file-list-q", config.q_file);
    if (yaml["outdir"]) {
        config.output_prefix = yaml["outdir"].as<std::string>() + "/" + config.material_name;
    }
}

void create_output_parent(const std::string& output_prefix) {
    const std::filesystem::path prefix_path(output_prefix);
    const auto                  parent = prefix_path.parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
}

bool has_nonempty_env(const char* name) {
    const char* value = std::getenv(name);
    return value != nullptr && value[0] != '\0';
}

bool launched_under_mpi() {
    return has_nonempty_env("OMPI_COMM_WORLD_SIZE") || has_nonempty_env("PMI_SIZE") || has_nonempty_env("PMIX_RANK") ||
           has_nonempty_env("MPI_LOCALNRANKS");
}

void normalize_mode(EpsilonAppConfig& config) {
    config.mode = lowercase_ascii(config.mode);
    std::replace(config.mode.begin(), config.mode.end(), '_', '-');
}

void print_config(const EpsilonAppConfig&              config,
                  const std::vector<double>&           energies,
                  const std::vector<Vector3D<double>>& qpoints,
                  int                                  number_processes) {
    std::cout << "EPSILON PROGRAM\n";
    std::cout << "MPI processes: " << number_processes << '\n';
    std::cout << "Material: " << config.material_name << '\n';
    std::cout << "EPM set: " << config.epm_parameter_set << '\n';
    std::cout << "Bands: " << config.nb_bands << '\n';
    std::cout << "Nearest neighbors: " << config.nb_nearest_neighbors << '\n';
    std::cout << "Nonlocal EPM: " << config.nonlocal_epm << '\n';
    std::cout << "k-grid: " << config.nkx << " x " << config.nky << " x " << config.nkz << '\n';
    std::cout << "BZ sampling: " << normalize_bz_sampling(config.bz_sampling) << '\n';
    std::cout << "Energy grid: " << energies.front() << " -> " << energies.back() << " eV, N=" << energies.size()
              << ", eta=" << config.eta_smearing_eV << " eV\n";
    std::cout << "Mode: " << config.mode << '\n';
    if (config.mode == "optical") {
        std::cout << "Optical polarization direction: " << normalized_direction(config.direction) << '\n';
    }
    std::cout << "q-points: " << qpoints.size() << '\n';
    if (!qpoints.empty() && config.mode != "optical") {
        auto [min_q, max_q] = std::minmax_element(qpoints.begin(), qpoints.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.Length() < rhs.Length();
        });
        std::cout << "q range: " << min_q->Length() << " -> " << max_q->Length() << '\n';
        std::cout << "First q: " << qpoints.front() << '\n';
        if (min_q->Length() < config.small_q_diagnostic_threshold) {
            std::cout << "Warning: small-q run. The finite-q density matrix element should scale as q and cancel the "
                         "Coulomb 1/q^2 factor. If epsilon diverges here, that is an implementation/gauge problem, "
                         "not proof that q is physically too small.\n";
        }
    }
    std::cout << "Output prefix: " << config.output_prefix << '\n';
}

std::vector<std::vector<std::vector<double>>> gather_dielectric_contributions(const std::vector<double>& flat_local,
                                                                              std::size_t                nb_qpoints,
                                                                              std::size_t                nb_energies,
                                                                              int number_processes,
                                                                              int process_rank) {
    const std::size_t qe      = nb_qpoints * nb_energies;
    const int         local_n = static_cast<int>(flat_local.size());

    std::vector<int> recvcounts;
    if (process_rank == 0) {
        recvcounts.resize(number_processes, 0);
    }
    MPI_Gather(&local_n, 1, MPI_INT, process_rank == 0 ? recvcounts.data() : nullptr, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int>    displacements;
    std::vector<double> all_flat;
    if (process_rank == 0) {
        displacements.resize(number_processes, 0);
        for (int p = 1; p < number_processes; ++p) {
            displacements[p] = displacements[p - 1] + recvcounts[p - 1];
        }
        const int total = number_processes > 0 ? displacements.back() + recvcounts.back() : 0;
        all_flat.resize(static_cast<std::size_t>(total));
    }

    MPI_Gatherv(flat_local.data(),
                local_n,
                MPI_DOUBLE,
                process_rank == 0 ? all_flat.data() : nullptr,
                process_rank == 0 ? recvcounts.data() : nullptr,
                process_rank == 0 ? displacements.data() : nullptr,
                MPI_DOUBLE,
                0,
                MPI_COMM_WORLD);

    if (process_rank != 0) {
        return {};
    }

    std::vector<std::vector<std::vector<double>>> by_process(
        number_processes,
        std::vector<std::vector<double>>(nb_qpoints, std::vector<double>(nb_energies, 0.0)));

    for (int p = 0; p < number_processes; ++p) {
        const int         count = recvcounts[p];
        const std::size_t base  = static_cast<std::size_t>(displacements[p]);
        if (count == static_cast<int>(qe)) {
            for (std::size_t q = 0; q < nb_qpoints; ++q) {
                const double* src = &all_flat[base + q * nb_energies];
                std::copy(src, src + nb_energies, by_process[p][q].data());
            }
        } else if (count % static_cast<int>(qe) == 0) {
            const int local_k_blocks = count / static_cast<int>(qe);
            for (int k = 0; k < local_k_blocks; ++k) {
                const std::size_t block = base + static_cast<std::size_t>(k) * qe;
                for (std::size_t q = 0; q < nb_qpoints; ++q) {
                    const double* src = &all_flat[block + q * nb_energies];
                    double*       dst = by_process[p][q].data();
                    for (std::size_t e = 0; e < nb_energies; ++e) {
                        dst[e] += src[e];
                    }
                }
            }
        } else {
            throw std::runtime_error("Unexpected dielectric payload size from rank " + std::to_string(p));
        }
    }
    return by_process;
}

std::vector<std::vector<std::vector<double>>> local_dielectric_contribution_as_single_process(
    const std::vector<double>& flat_local,
    std::size_t                nb_qpoints,
    std::size_t                nb_energies) {
    const std::size_t expected_size = nb_qpoints * nb_energies;
    if (flat_local.size() != expected_size) {
        throw std::runtime_error("Unexpected serial dielectric payload size.");
    }

    std::vector<std::vector<std::vector<double>>> result(
        1,
        std::vector<std::vector<double>>(nb_qpoints, std::vector<double>(nb_energies, 0.0)));
    for (std::size_t q = 0; q < nb_qpoints; ++q) {
        const double* src = &flat_local[q * nb_energies];
        std::copy(src, src + nb_energies, result[0][q].data());
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    TCLAP::CmdLine               cmd("Compute dynamic dielectric functions epsilon(q,E).", ' ', "0.2");
    TCLAP::ValueArg<std::string> arg_config("c", "config", "Optional YAML config file.", false, "", "path");
    TCLAP::ValueArg<std::string> arg_material("m", "material", "Material symbol.", false, "", "symbol");
    TCLAP::ValueArg<std::string> arg_epm_set("", "epm-set", "Named EPM parameter set.", false, "", "name");
    TCLAP::ValueArg<std::string> arg_mode("", "mode", "Mode: optical, q-list, q-line, or q-file.", false, "", "mode");
    TCLAP::ValueArg<std::string> arg_q_values("",
                                              "q-values",
                                              "Comma-separated q magnitudes in reduced units.",
                                              false,
                                              "",
                                              "list");
    TCLAP::ValueArg<std::string> arg_direction("", "direction", "q direction, e.g. 1,0,0.", false, "", "vector");
    TCLAP::ValueArg<std::string> arg_q_file("",
                                            "q-file",
                                            "File containing qx qy qz rows in reduced units.",
                                            false,
                                            "",
                                            "path");
    TCLAP::ValueArg<std::string> arg_out("o", "out", "Output file prefix.", false, "", "prefix");
    TCLAP::ValueArg<int>         arg_bands("b", "bands", "Number of EPM bands.", false, -1, "int");
    TCLAP::ValueArg<int> arg_neighbors("", "nearest-neighbors", "Number of reciprocal shells.", false, -1, "int");
    TCLAP::ValueArg<int> arg_nkx("", "Nkx", "k-grid count in x.", false, -1, "int");
    TCLAP::ValueArg<int> arg_nky("", "Nky", "k-grid count in y.", false, -1, "int");
    TCLAP::ValueArg<int> arg_nkz("", "Nkz", "k-grid count in z.", false, -1, "int");
    TCLAP::ValueArg<std::string> arg_bz_sampling(
        "",
        "bz-sampling",
        "BZ sampling: full, q100-octant, or fcc-ibz. Legacy aliases: 1, 8, 48.",
        false,
        "",
        "mode");
    TCLAP::ValueArg<int>    arg_q_count("", "q-count", "Number of q samples for q-line mode.", false, -1, "int");
    TCLAP::ValueArg<double> arg_emin("", "emin", "Minimum energy in eV.", false, std::nan(""), "eV");
    TCLAP::ValueArg<double> arg_emax("", "emax", "Maximum energy in eV.", false, std::nan(""), "eV");
    TCLAP::ValueArg<double> arg_estep("", "estep", "Energy step in eV.", false, std::nan(""), "eV");
    TCLAP::ValueArg<double> arg_eta("", "eta", "Smearing in eV.", false, std::nan(""), "eV");
    TCLAP::ValueArg<double> arg_q_min("", "q-min", "Minimum q magnitude for q-line mode.", false, std::nan(""), "q");
    TCLAP::ValueArg<double> arg_q_max("", "q-max", "Maximum q magnitude for q-line mode.", false, std::nan(""), "q");
    TCLAP::SwitchArg        arg_nonlocal("", "nonlocal", "Use nonlocal EPM parameters.", false);

    cmd.add(arg_config);
    cmd.add(arg_material);
    cmd.add(arg_epm_set);
    cmd.add(arg_mode);
    cmd.add(arg_q_values);
    cmd.add(arg_direction);
    cmd.add(arg_q_file);
    cmd.add(arg_out);
    cmd.add(arg_bands);
    cmd.add(arg_neighbors);
    cmd.add(arg_nkx);
    cmd.add(arg_nky);
    cmd.add(arg_nkz);
    cmd.add(arg_bz_sampling);
    cmd.add(arg_q_count);
    cmd.add(arg_emin);
    cmd.add(arg_emax);
    cmd.add(arg_estep);
    cmd.add(arg_eta);
    cmd.add(arg_q_min);
    cmd.add(arg_q_max);
    cmd.add(arg_nonlocal);

    if (argc == 1) {
        std::cout << "Usage example:\n"
                  << "  epsilon.epm --material Si --bands 16 --mode optical "
                     "--direction 1,0,0 --emin 0 --emax 20 --estep 0.05 --eta 0.05 "
                     "--Nkx 20 --Nky 20 --Nkz 20 --out results/epsilon_si_optical\n\n"
                  << "Run epsilon.epm --help for all options.\n";
        return 0;
    }

    cmd.parse(argc, argv);

    EpsilonAppConfig config;
    load_config_file(arg_config.getValue(), config);

    if (arg_material.isSet()) {
        config.material_name = arg_material.getValue();
    }
    if (arg_epm_set.isSet()) {
        config.epm_parameter_set = arg_epm_set.getValue();
    }
    if (arg_mode.isSet()) {
        config.mode = arg_mode.getValue();
    }
    if (arg_q_values.isSet()) {
        config.q_values = arg_q_values.getValue();
    }
    if (arg_direction.isSet()) {
        config.direction = parse_vector3(arg_direction.getValue());
    }
    if (arg_q_file.isSet()) {
        config.q_file = arg_q_file.getValue();
        config.mode   = "q-file";
    }
    if (arg_out.isSet()) {
        config.output_prefix = arg_out.getValue();
    }
    if (arg_bands.isSet()) {
        config.nb_bands = arg_bands.getValue();
    }
    if (arg_neighbors.isSet()) {
        config.nb_nearest_neighbors = arg_neighbors.getValue();
    }
    if (arg_nkx.isSet()) {
        config.nkx = arg_nkx.getValue();
    }
    if (arg_nky.isSet()) {
        config.nky = arg_nky.getValue();
    }
    if (arg_nkz.isSet()) {
        config.nkz = arg_nkz.getValue();
    }
    if (arg_bz_sampling.isSet()) {
        config.bz_sampling = arg_bz_sampling.getValue();
    }
    if (arg_q_count.isSet()) {
        config.q_count = arg_q_count.getValue();
    }
    if (arg_emin.isSet()) {
        config.min_energy_eV = arg_emin.getValue();
    }
    if (arg_emax.isSet()) {
        config.max_energy_eV = arg_emax.getValue();
    }
    if (arg_estep.isSet()) {
        config.energy_step_eV = arg_estep.getValue();
    }
    if (arg_eta.isSet()) {
        config.eta_smearing_eV = arg_eta.getValue();
    }
    if (arg_q_min.isSet()) {
        config.q_min = arg_q_min.getValue();
    }
    if (arg_q_max.isSet()) {
        config.q_max = arg_q_max.getValue();
    }
    if (arg_nonlocal.isSet()) {
        config.nonlocal_epm = true;
        if (!arg_epm_set.isSet()) {
            config.epm_parameter_set = "potz-vogl";
        }
    }
    if ((arg_mode.isSet() || arg_q_values.isSet()) && !arg_q_file.isSet() && config.mode != "q-file") {
        config.q_file.clear();
    }
    normalize_mode(config);
    if (config.mode != "optical" && config.mode != "q-list" && config.mode != "q-line" && config.mode != "q-file") {
        throw std::invalid_argument("Unknown epsilon mode '" + config.mode +
                                    "'. Use optical, q-list, q-line, or q-file.");
    }
    config.bz_sampling = normalize_bz_sampling(config.bz_sampling);

    const bool use_mpi          = launched_under_mpi();
    int        number_processes = 1;
    int        process_rank     = 0;
    if (use_mpi) {
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &number_processes);
        MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

        char processor_name[MPI_MAX_PROCESSOR_NAME];
        int  processor_name_len = 0;
        MPI_Get_processor_name(processor_name, &processor_name_len);
        std::cout << "Rank " << process_rank << "/" << number_processes << " on " << processor_name << '\n';
    } else {
        std::cout << "Running in serial mode. Launch with mpirun to use MPI.\n";
    }

    const auto energies = make_energy_grid(config.min_energy_eV, config.max_energy_eV, config.energy_step_eV);
    const auto qpoints  = make_qpoints(config);

    if (process_rank == 0) {
        create_output_parent(config.output_prefix);
        print_config(config, energies, qpoints, number_processes);
        if (config.bz_sampling == "fcc-ibz") {
            std::cout
                << "Warning: fcc-ibz sampling is not generally valid for finite-q or direction-resolved optical "
                   "dielectric functions unless symmetry weights and q/polarization-star handling are implemented. "
                   "Full-BZ sampling (--bz-sampling full) is recommended.\n";
        }
        const bool q100_direction_ok =
            config.mode == "optical" ? vector_along_100(config.direction) : all_qpoints_along_100(qpoints);
        if (config.bz_sampling == "q100-octant" && !q100_direction_ok) {
            std::cout << "Warning: --bz-sampling q100-octant assumes q is parallel to the [100] direction. "
                         "Use --bz-sampling full for arbitrary q directions.\n";
        }
    }

    uepm::pseudopotential::Materials         materials;
    const uepm::physics::material_repository material_repository;
    materials.load_material(material_repository, config.material_name, config.epm_parameter_set);
    uepm::pseudopotential::epm_material current_material = materials.materials.at(config.material_name);

    uepm::pseudopotential::BandStructure band_structure{};
    band_structure.Initialize(current_material,
                              config.nb_bands,
                              {},
                              config.nb_nearest_neighbors,
                              config.nonlocal_epm,
                              config.enable_soc);

    uepm::pseudopotential::DielectricFunction dielectric(current_material,
                                                         band_structure.get_basis_vectors(),
                                                         config.nb_bands);

    dielectric.generate_k_points_grid(config.nkx,
                                      config.nky,
                                      config.nkz,
                                      0.0,
                                      dielectric_sampling_mode(config.bz_sampling));
    const std::size_t nb_kpoints = dielectric.get_kpoints().size();

    if (process_rank == 0) {
        dielectric.export_kpoints(config.output_prefix + "_kpoints.csv");
        std::cout << "Generated k-points: " << nb_kpoints << '\n';
    }

    const int total_kpoints = static_cast<int>(nb_kpoints);
    const int base          = total_kpoints / number_processes;
    const int remainder     = total_kpoints % number_processes;

    std::vector<int> counts_kpoints_per_process(number_processes, 0);
    std::vector<int> displacements_kpoints_per_process(number_processes, 0);
    for (int p = 0; p < number_processes; ++p) {
        counts_kpoints_per_process[p] = base + (p < remainder ? 1 : 0);
        if (p > 0) {
            displacements_kpoints_per_process[p] =
                displacements_kpoints_per_process[p - 1] + counts_kpoints_per_process[p - 1];
        }
    }

    std::cout << "Rank " << process_rank << " handles " << counts_kpoints_per_process[process_rank] << " k-points\n";

    dielectric.set_export_prefix(config.output_prefix);
    dielectric.set_qpoints(qpoints);
    if (config.mode == "optical") {
        dielectric.set_response_mode(uepm::pseudopotential::DielectricResponseMode::optical_limit);
        dielectric.set_optical_direction(config.direction);
    }
    dielectric.set_energies(energies);
    dielectric.set_offset_k_index(static_cast<std::size_t>(displacements_kpoints_per_process[process_rank]));
    dielectric.set_nb_kpoints(static_cast<std::size_t>(counts_kpoints_per_process[process_rank]));
    dielectric.set_non_local_epm(config.nonlocal_epm);

    if (use_mpi) {
        MPI_Barrier(MPI_COMM_WORLD);
    }
    const auto start = std::chrono::high_resolution_clock::now();
    dielectric.compute_dielectric_function(config.eta_smearing_eV, process_rank);
    dielectric.clear_eigen_states();
    if (use_mpi) {
        MPI_Barrier(MPI_COMM_WORLD);
    }
    const auto end = std::chrono::high_resolution_clock::now();
    if (process_rank == 0) {
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Total computational time: " << elapsed.count() << " s\n";
    }

    const auto local_flat_real = dielectric.get_flat_dielectric_function();
    const auto local_flat_imag = dielectric.get_flat_dielectric_function_imaginary();
    auto       dielectric_real_by_process =
        use_mpi ? gather_dielectric_contributions(local_flat_real,
                                                  qpoints.size(),
                                                  energies.size(),
                                                  number_processes,
                                                  process_rank)
                      : local_dielectric_contribution_as_single_process(local_flat_real, qpoints.size(), energies.size());
    auto dielectric_imag_by_process =
        use_mpi ? gather_dielectric_contributions(local_flat_imag,
                                                  qpoints.size(),
                                                  energies.size(),
                                                  number_processes,
                                                  process_rank)
                : local_dielectric_contribution_as_single_process(local_flat_imag, qpoints.size(), energies.size());

    if (process_rank == 0) {
        auto merged = uepm::pseudopotential::DielectricFunction::merge_results(dielectric,
                                                                               dielectric_real_by_process,
                                                                               dielectric_imag_by_process,
                                                                               counts_kpoints_per_process);
        merged.export_dielectric_function("", true);
        std::cout << "Wrote epsilon spectra with prefix '" << config.output_prefix << "'.\n";
    }

    if (use_mpi) {
        MPI_Finalize();
    }
    return 0;
}
