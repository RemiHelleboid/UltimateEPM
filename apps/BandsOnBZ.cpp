/**
 * @file compute_bands_on_bz_mesh.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-07-07
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <tclap/CmdLine.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <vector>

#include "BandStructure.h"
#include "Material.h"
#include "Options.h"
#include "bz_mesh.hpp"
#include "bz_meshfile.hpp"

namespace {

using uepm::mesh_bz::MeshBZ;
using uepm::mesh_bz::MeshParticleType;
using uepm::mesh_bz::Tetra;
using uepm::mesh_bz::vector3;

struct SpreadSample {
    double      spread{0.0};
    double      volume{0.0};
    std::size_t tetraIndex{0};
};

struct RefinementCandidate {
    double      error{0.0};
    double      targetH{0.0};
    double      volume{0.0};
    double      quality{0.0};
    std::size_t bandIndex{0};
    std::size_t tetraIndex{0};
    vector3     barycenterReduced{};
};

struct AdaptiveSummary {
    std::size_t ibzTetrahedra{0};
    std::size_t violatingTetrahedra{0};
    double      violatingVolumeFraction{0.0};
    double      p50Error{0.0};
    double      p90Error{0.0};
    double      p99Error{0.0};
    double      maxError{0.0};
};

double percentile(const std::vector<double>& sortedValues, double p) {
    if (sortedValues.empty()) {
        return 0.0;
    }

    const double x = p * static_cast<double>(sortedValues.size() - 1);
    const auto   i = static_cast<std::size_t>(std::floor(x));
    const auto   j = std::min(i + 1, sortedValues.size() - 1);
    const double t = x - static_cast<double>(i);
    return (1.0 - t) * sortedValues[i] + t * sortedValues[j];
}

double volume_weighted_percentile(std::vector<SpreadSample> samples, double p) {
    if (samples.empty()) {
        return 0.0;
    }

    std::sort(samples.begin(), samples.end(), [](const auto& a, const auto& b) { return a.spread < b.spread; });
    const double totalVolume =
        std::accumulate(samples.begin(), samples.end(), 0.0, [](double sum, const auto& s) { return sum + s.volume; });
    const double target = p * totalVolume;
    double       cumulative = 0.0;

    for (const auto& sample : samples) {
        cumulative += sample.volume;
        if (cumulative >= target) {
            return sample.spread;
        }
    }

    return samples.back().spread;
}

double tetra_quality(const Tetra& tetra) {
    const auto& vertices = tetra.get_list_vertices();
    double      sumEdgeSquared = 0.0;

    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = i + 1; j < 4; ++j) {
            sumEdgeSquared += (vertices[i]->get_position() - vertices[j]->get_position()).norm_squared();
        }
    }

    const double volume = std::abs(tetra.get_signed_volume());
    if (volume <= 0.0 || sumEdgeSquared <= 0.0) {
        return 0.0;
    }

    return 12.0 * std::pow(3.0 * volume, 2.0 / 3.0) / sumEdgeSquared;
}

double max_edge_length_reduced(const MeshBZ& mesh, const Tetra& tetra) {
    const auto& vertices = tetra.get_list_vertices();
    double      maxLength = 0.0;

    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = i + 1; j < 4; ++j) {
            const vector3 edgeReduced = mesh.si_to_reduced_k(vertices[i]->get_position() - vertices[j]->get_position());
            maxLength                 = std::max(maxLength, edgeReduced.norm());
        }
    }

    return maxLength;
}

std::vector<RefinementCandidate> spatially_diverse_candidates(std::vector<RefinementCandidate> candidates,
                                                               std::size_t                      maxPoints) {
    std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) { return a.error > b.error; });
    if (candidates.size() <= maxPoints) {
        return candidates;
    }

    std::vector<RefinementCandidate> selected;
    selected.reserve(maxPoints);

    std::vector<double> minDistanceSquared(candidates.size(), std::numeric_limits<double>::max());
    std::vector<bool>   wasSelected(candidates.size(), false);
    std::size_t         selectedIndex = 0;
    const double        maxError      = candidates.front().error;

    for (std::size_t selectedCount = 0; selectedCount < maxPoints; ++selectedCount) {
        wasSelected[selectedIndex] = true;
        selected.push_back(candidates[selectedIndex]);

        const vector3 selectedPosition = candidates[selectedIndex].barycenterReduced;
        double        bestScore        = -1.0;
        std::size_t   bestIndex        = candidates.size();

        for (std::size_t i = 0; i < candidates.size(); ++i) {
            if (wasSelected[i]) {
                continue;
            }

            const vector3 delta = candidates[i].barycenterReduced - selectedPosition;
            minDistanceSquared[i] = std::min(minDistanceSquared[i], delta.norm_squared());

            // Error remains important, but distance prevents the capped map
            // from repeatedly selecting a dense cluster of adjacent tetrahedra.
            const double errorWeight = std::sqrt(candidates[i].error / maxError);
            const double score       = minDistanceSquared[i] * errorWeight;
            if (score > bestScore) {
                bestScore = score;
                bestIndex = i;
            }
        }

        if (bestIndex == candidates.size()) {
            break;
        }
        selectedIndex = bestIndex;
    }

    return selected;
}

std::filesystem::path adaptive_summary_path(const std::filesystem::path& refinementMap) {
    return refinementMap.parent_path() / (refinementMap.stem().string() + "_summary.csv");
}

void write_adaptive_summary(const std::filesystem::path& refinementMap,
                            double                       targetEnergy,
                            const AdaptiveSummary&       summary,
                            std::size_t                  exportedPoints) {
    const std::filesystem::path summaryPath = adaptive_summary_path(refinementMap);
    std::ofstream               out(summaryPath);
    if (!out) {
        throw std::runtime_error("Could not open adaptive summary " + summaryPath.string());
    }

    out << "target_eV,ibz_tetrahedra,violating_tetrahedra,violating_volume_fraction,"
           "p50_error_eV,p90_error_eV,p99_error_eV,max_error_eV,exported_points\n";
    out << std::setprecision(17) << targetEnergy << ',' << summary.ibzTetrahedra << ','
        << summary.violatingTetrahedra << ',' << summary.violatingVolumeFraction << ',' << summary.p50Error << ','
        << summary.p90Error << ',' << summary.p99Error << ',' << summary.maxError << ',' << exportedPoints << '\n';

    fmt::print("Wrote adaptive summary to {}\n", summaryPath.string());
}

std::vector<std::size_t> selected_export_bands(const MeshBZ& mesh, std::size_t nbValence, std::size_t nbConduction) {
    std::vector<std::size_t> bands;

    const auto valence    = mesh.get_band_indices(MeshParticleType::valence);
    const auto conduction = mesh.get_band_indices(MeshParticleType::conduction);

    nbValence    = std::min(nbValence, valence.size());
    nbConduction = std::min(nbConduction, conduction.size());

    bands.insert(bands.end(), valence.end() - static_cast<std::ptrdiff_t>(nbValence), valence.end());
    bands.insert(bands.end(), conduction.begin(), conduction.begin() + static_cast<std::ptrdiff_t>(nbConduction));
    return bands;
}

void write_energy_adaptive_diagnostics(const MeshBZ&                 mesh,
                                       const std::vector<std::size_t>& bands,
                                       double                        targetEnergy,
                                       std::size_t                   maxRefinementPoints,
                                       double                        backgroundH,
                                       const std::filesystem::path&  refinementMap) {
    if (bands.empty()) {
        throw std::runtime_error("No bands were selected for adaptive diagnostics.");
    }
    if (!(targetEnergy > 0.0)) {
        throw std::runtime_error("Adaptive energy target must be positive.");
    }

    const auto& tetrahedra = mesh.get_list_tetrahedra();
    std::vector<RefinementCandidate> candidates;
    candidates.reserve(std::min<std::size_t>(tetrahedra.size(), maxRefinementPoints * 4));
    std::vector<SpreadSample> worstBandSamples;
    worstBandSamples.reserve(tetrahedra.size() / 48 + 1);

    constexpr std::array<double, 4> thresholds{0.010, 0.020, 0.050, 0.100};

    for (const std::size_t band : bands) {
        std::vector<SpreadSample> samples;
        samples.reserve(tetrahedra.size() / 48 + 1);

        for (const auto& tetra : tetrahedra) {
            if (!tetra.lies_in_irreducible_wedge()) {
                continue;
            }

            const auto energies = tetra.get_band_energies_at_vertices(band);
            const auto minmax   = std::minmax_element(energies.begin(), energies.end());
            const double spread = *minmax.second - *minmax.first;
            const double volume = std::abs(tetra.get_signed_volume());
            samples.push_back({spread, volume, tetra.get_index()});
        }

        if (samples.empty()) {
            throw std::runtime_error("No irreducible-wedge tetrahedra found for adaptive diagnostics.");
        }

        std::vector<double> sortedSpreads;
        sortedSpreads.reserve(samples.size());
        for (const auto& sample : samples) {
            sortedSpreads.push_back(sample.spread);
        }
        std::sort(sortedSpreads.begin(), sortedSpreads.end());

        const double totalVolume =
            std::accumulate(samples.begin(), samples.end(), 0.0, [](double sum, const auto& s) { return sum + s.volume; });
        const auto worst = std::max_element(samples.begin(), samples.end(), [](const auto& a, const auto& b) {
            return a.spread < b.spread;
        });
        const Tetra&  worstTetra      = tetrahedra[worst->tetraIndex];
        const vector3 worstBarycenter = mesh.si_to_reduced_k(worstTetra.get_barycenter());

        fmt::print("Band {} IBZ deltaE: p50={:.6f}, p90={:.6f}, p99={:.6f}, p99.9={:.6f}, max={:.6f} eV\n",
                   band,
                   percentile(sortedSpreads, 0.50),
                   percentile(sortedSpreads, 0.90),
                   percentile(sortedSpreads, 0.99),
                   percentile(sortedSpreads, 0.999),
                   sortedSpreads.back());
        fmt::print("Band {} volume-weighted deltaE: p50={:.6f}, p90={:.6f}, p99={:.6f}, p99.9={:.6f} eV\n",
                   band,
                   volume_weighted_percentile(samples, 0.50),
                   volume_weighted_percentile(samples, 0.90),
                   volume_weighted_percentile(samples, 0.99),
                   volume_weighted_percentile(samples, 0.999));

        for (const double threshold : thresholds) {
            const std::size_t count = static_cast<std::size_t>(
                std::count_if(samples.begin(), samples.end(), [threshold](const auto& s) { return s.spread > threshold; }));
            const double volumeAbove =
                std::accumulate(samples.begin(), samples.end(), 0.0, [threshold](double sum, const auto& s) {
                    return sum + ((s.spread > threshold) ? s.volume : 0.0);
                });
            fmt::print("  >{:.0f} meV: {} tets ({:.4f}% of IBZ volume)\n",
                       1000.0 * threshold,
                       count,
                       100.0 * volumeAbove / totalVolume);
        }

        fmt::print("  worst tetra={} barycenter=({:.9f},{:.9f},{:.9f}) volume={:.6e} quality={:.6f}\n",
                   worst->tetraIndex,
                   worstBarycenter.x(),
                   worstBarycenter.y(),
                   worstBarycenter.z(),
                   worst->volume,
                   tetra_quality(worstTetra));
    }

    for (const auto& tetra : tetrahedra) {
        if (!tetra.lies_in_irreducible_wedge()) {
            continue;
        }

        double      worstSpread = 0.0;
        std::size_t worstBand   = bands.front();
        for (const std::size_t band : bands) {
            const auto energies = tetra.get_band_energies_at_vertices(band);
            const auto minmax   = std::minmax_element(energies.begin(), energies.end());
            const double spread = *minmax.second - *minmax.first;
            if (spread > worstSpread) {
                worstSpread = spread;
                worstBand   = band;
            }
        }

        const double volume = std::abs(tetra.get_signed_volume());
        worstBandSamples.push_back({worstSpread, volume, tetra.get_index()});
        if (worstSpread <= targetEnergy) {
            continue;
        }

        const double currentH = max_edge_length_reduced(mesh, tetra);
        const double scale    = std::clamp(std::sqrt(targetEnergy / worstSpread), 0.25, 0.80);
        const double targetH =
            (backgroundH > 0.0) ? std::min(currentH * scale, 0.80 * backgroundH) : currentH * scale;
        candidates.push_back({worstSpread,
                              targetH,
                              volume,
                              tetra_quality(tetra),
                              worstBand,
                              tetra.get_index(),
                              mesh.si_to_reduced_k(tetra.get_barycenter())});
    }

    AdaptiveSummary summary;
    summary.ibzTetrahedra = static_cast<std::size_t>(
        std::count_if(tetrahedra.begin(), tetrahedra.end(), [](const auto& tetra) {
            return tetra.lies_in_irreducible_wedge();
        }));
    summary.violatingTetrahedra = candidates.size();

    if (!worstBandSamples.empty()) {
        std::vector<double> sortedErrors;
        sortedErrors.reserve(worstBandSamples.size());
        for (const auto& sample : worstBandSamples) {
            sortedErrors.push_back(sample.spread);
        }
        std::sort(sortedErrors.begin(), sortedErrors.end());

        const double violatingVolume =
            std::accumulate(candidates.begin(), candidates.end(), 0.0, [](double sum, const auto& candidate) {
                return sum + candidate.volume;
            });
        const double ibzVolume = std::accumulate(
            tetrahedra.begin(), tetrahedra.end(), 0.0, [](double sum, const auto& tetra) {
                return sum + (tetra.lies_in_irreducible_wedge() ? std::abs(tetra.get_signed_volume()) : 0.0);
            });

        summary.violatingVolumeFraction = violatingVolume / ibzVolume;
        summary.p50Error                = percentile(sortedErrors, 0.50);
        summary.p90Error                = percentile(sortedErrors, 0.90);
        summary.p99Error                = percentile(sortedErrors, 0.99);
        summary.maxError                = sortedErrors.back();
    }

    const std::size_t totalCandidates = candidates.size();
    candidates = spatially_diverse_candidates(std::move(candidates), maxRefinementPoints);

    if (refinementMap.has_parent_path()) {
        std::filesystem::create_directories(refinementMap.parent_path());
    }
    std::ofstream out(refinementMap);
    if (!out) {
        throw std::runtime_error("Could not open refinement map " + refinementMap.string());
    }

    out << "x,y,z,error_eV,target_h,band,tetra,volume,quality\n";
    out << std::setprecision(17);
    for (const auto& candidate : candidates) {
        out << candidate.barycenterReduced.x() << ',' << candidate.barycenterReduced.y() << ','
            << candidate.barycenterReduced.z() << ',' << candidate.error << ',' << candidate.targetH << ','
            << candidate.bandIndex << ',' << candidate.tetraIndex << ',' << candidate.volume << ',' << candidate.quality
            << '\n';
    }

    write_adaptive_summary(refinementMap, targetEnergy, summary, candidates.size());
    fmt::print("Adaptive target: {} / {} IBZ tetrahedra violate {:.6f} eV ({:.4f}% of IBZ volume)\n",
               totalCandidates,
               summary.ibzTetrahedra,
               targetEnergy,
               100.0 * summary.violatingVolumeFraction);
    fmt::print("Wrote {} spatially diverse adaptive refinement points to {}\n",
               candidates.size(),
               refinementMap.string());
}

}  // namespace

int main(int argc, char* argv[]) {
    TCLAP::CmdLine               cmd("EPP PROGRAM. COMPUTE BAND STRUCTURE ON A BZ MESH.", ' ', "1.0");
    TCLAP::ValueArg<std::string> arg_mesh_file("f", "meshfile", "Name to print", true, "bz.msh", "string");
    TCLAP::ValueArg<std::string> arg_material("m", "material", "Symbol of the material to use (Si, Ge, GaAs, ...)", true, "Si", "string");
    TCLAP::ValueArg<std::string> arg_data_mat("d", "file-data", "Material data file", false, "materials-local-cohen.yaml", "string");
    TCLAP::ValueArg<std::string> arg_outfile("o", "outfile", "Name of the output file", false, "", "string");
    TCLAP::ValueArg<int>         arg_nb_valence_bands("v", "nvbands", "Number of valence bands to export", false, 4, "int");
    TCLAP::ValueArg<int>         arg_nb_conduction_bands("c", "ncbands", "Number of conduction bands to export", false, 12, "int");
    TCLAP::ValueArg<int>         arg_nearest_neighbors("n",
                                               "nearestNeighbors",
                                               "number of nearest neiborgs to consider for the EPP calculation.",
                                               false,
                                               10,
                                               "int");
    TCLAP::SwitchArg arg_enable_nonlocal_correction("C", "nonlocal-correction", "Enable the non-local-correction for the EPM model", false);
    TCLAP::SwitchArg arg_enable_soc("s", "soc", "Enable the spin-orbit coupling for the EPM model", false);
    TCLAP::SwitchArg arg_cond_band_zero("z", "MinCondZero", "Shift the conduction band minimum to 0 eV", false);
    TCLAP::SwitchArg arg_irr_wedge("w", "IrrWedge", "Compute bands only in the irreducible wedge of the BZ", false);
    TCLAP::ValueArg<int> arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    TCLAP::ValueArg<std::string> arg_refinement_map("",
                                                    "refinement-map",
                                                    "Adaptive refinement CSV output (default: <mesh>_refinement.csv)",
                                                    false,
                                                    "",
                                                    "string");
    TCLAP::ValueArg<double> arg_adaptive_energy_target("",
                                                       "adaptive-energy-target",
                                                       "Tetra energy-spread target in eV",
                                                       false,
                                                       0.020,
                                                       "float");
    TCLAP::ValueArg<int> arg_adaptive_max_points("",
                                                 "adaptive-max-points",
                                                 "Maximum spatially diverse refinement points to export",
                                                 false,
                                                 5000,
                                                 "int");
    TCLAP::ValueArg<double> arg_adaptive_background_h(
        "",
        "adaptive-background-h",
        "Background mesh size; exported targets are kept below this value (0 disables)",
        false,
        0.0,
        "float");
    cmd.add(arg_mesh_file);
    cmd.add(arg_material);
    cmd.add(arg_nb_valence_bands);
    cmd.add(arg_nb_conduction_bands);
    cmd.add(arg_outfile);
    cmd.add(arg_nearest_neighbors);
    cmd.add(arg_nb_threads);
    cmd.add(arg_enable_nonlocal_correction);
    cmd.add(arg_enable_soc);
    cmd.add(arg_cond_band_zero);
    cmd.add(arg_data_mat);
    cmd.add(arg_irr_wedge);
    cmd.add(arg_refinement_map);
    cmd.add(arg_adaptive_energy_target);
    cmd.add(arg_adaptive_max_points);
    cmd.add(arg_adaptive_background_h);

    cmd.parse(argc, argv);

    uepm::pseudopotential::Materials materials;

    std::string file_material_parameters = arg_data_mat.getValue();
    if (!std::filesystem::exists(file_material_parameters)) {
        std::filesystem::path p_try = std::filesystem::path(PROJECT_SRC_DIR) / "parameter_files" / file_material_parameters;
        if (std::filesystem::exists(p_try)) {
            file_material_parameters = p_try.string();
        } else {
            std::cerr << "Error: material data file " << file_material_parameters << " does not exist!" << std::endl;
            return -1;
        }
    }
    fmt::print("Loading material parameters from file: {}\n", file_material_parameters);

    materials.load_material_parameters(file_material_parameters);
    bool enable_nonlocal_correction = arg_enable_nonlocal_correction.isSet();
    bool enable_soc                 = arg_enable_soc.isSet();

    Options my_options;
    my_options.materialName = arg_material.getValue();

    const int max_valence_bands = enable_soc ? 8 : 4;

    my_options.nrLevels         = arg_nb_valence_bands.getValue() + arg_nb_conduction_bands.getValue() +
                          max_valence_bands;  // add extra bands for the calculations (won't be exported)
    my_options.nearestNeighbors = arg_nearest_neighbors.getValue();
    my_options.nrThreads        = arg_nb_threads.getValue();
    my_options.print_options();

    uepm::pseudopotential::Material mat = materials.materials.at(my_options.materialName);

    const std::string     mesh_filename = arg_mesh_file.getValue();
    uepm::mesh_bz::MeshBZ my_bz_mesh{mat};
    my_bz_mesh.set_number_threads_mesh_ops(my_options.nrThreads);

    my_bz_mesh.read_mesh_geometry_from_msh_file(mesh_filename);

    auto                                 start = std::chrono::high_resolution_clock::now();
    std::vector<Vector3D<double>>        mesh_kpoints{};
    uepm::pseudopotential::BandStructure my_bandstructure;
    my_bandstructure
        .Initialize(mat, my_options.nrLevels, mesh_kpoints, my_options.nearestNeighbors, enable_nonlocal_correction, arg_enable_soc);
    bool use_iwedge = arg_irr_wedge.isSet();
    my_bz_mesh.compute_band_structure_over_mesh(my_bandstructure, use_iwedge);

    const std::filesystem::path mesh_path(mesh_filename);
    const std::filesystem::path refinement_map =
        arg_refinement_map.isSet()
            ? std::filesystem::path(arg_refinement_map.getValue())
            : mesh_path.parent_path() / (mesh_path.stem().string() + "_refinement.csv");
    const auto diagnostic_bands = selected_export_bands(my_bz_mesh,
                                                        static_cast<std::size_t>(std::max(0, arg_nb_valence_bands.getValue())),
                                                        static_cast<std::size_t>(std::max(0, arg_nb_conduction_bands.getValue())));
    write_energy_adaptive_diagnostics(my_bz_mesh,
                                      diagnostic_bands,
                                      arg_adaptive_energy_target.getValue(),
                                      static_cast<std::size_t>(std::max(0, arg_adaptive_max_points.getValue())),
                                      arg_adaptive_background_h.getValue(),
                                      refinement_map);

    auto end              = std::chrono::high_resolution_clock::now();
    auto total_time_count = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Total computation time: " << total_time_count / double(1000) << std::endl;

    std::filesystem::path in_path(mesh_filename);
    std::string           out_file_bands = in_path.stem().replace_extension("").string() + "_" + my_bandstructure.path_band_filename();

    if (arg_outfile.isSet()) {
        out_file_bands = arg_outfile.getValue();
    }
    std::cout << "Exporting " << arg_nb_valence_bands.getValue() << " valence bands and " << arg_nb_conduction_bands.getValue()
              << " conduction bands to file: " << out_file_bands << std::endl;
    bool highest_valence_as_band0 = true;
    bool write_gradients         = true;
    my_bz_mesh.export_selected_bands_to_gmsh(out_file_bands,
                                             arg_nb_valence_bands.getValue(),
                                             arg_nb_conduction_bands.getValue(),
                                             highest_valence_as_band0,
                                             mesh_filename,
                                             write_gradients);

    const std::string vtk_file = in_path.stem().replace_extension("").string() + "_bands.vtk";
    my_bz_mesh.export_energies_and_gradients_to_vtk(vtk_file);

    return 0;
}
