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
#include "Options.h"
#include "bz_mesh.hpp"
#include "bz_meshfile.hpp"
#include "epm_material.hpp"
#include "pbmc_material_model.hpp"
#include "physical_constants.hpp"

namespace {

using uepm::mesh_bz::MeshBZ;
using uepm::mesh_bz::MeshParticleType;
using uepm::mesh_bz::Tetra;
using uepm::mesh_bz::vector3;

struct PbmcBandData {
    std::vector<double>  heavyHoleEnergy;
    std::vector<vector3> heavyHoleGradient;
    std::vector<double>  lightHoleEnergy;
    std::vector<vector3> lightHoleGradient;
    std::vector<double>  electronEnergy;
    std::vector<vector3> electronGradient;
    std::vector<double>  electronValley;
};

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
    const double target     = p * totalVolume;
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
    const auto& vertices       = tetra.get_list_vertices();
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
    const auto& vertices  = tetra.get_list_vertices();
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

            const vector3 delta   = candidates[i].barycenterReduced - selectedPosition;
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
    out << std::setprecision(17) << targetEnergy << ',' << summary.ibzTetrahedra << ',' << summary.violatingTetrahedra
        << ',' << summary.violatingVolumeFraction << ',' << summary.p50Error << ',' << summary.p90Error << ','
        << summary.p99Error << ',' << summary.maxError << ',' << exportedPoints << '\n';

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

uepm::PBMC::valley_model::vector3 to_pbmc_vector(const vector3& value) { return {value.x(), value.y(), value.z()}; }

vector3 to_bz_vector(const uepm::PBMC::valley_model::vector3& value) { return {value.x(), value.y(), value.z()}; }

vector3 delta_valley_center_reduced(const std::string& name, double deltaPosition) {
    if (name == "Delta_x_plus") {
        return {deltaPosition, 0.0, 0.0};
    }
    if (name == "Delta_x_minus") {
        return {-deltaPosition, 0.0, 0.0};
    }
    if (name == "Delta_y_plus") {
        return {0.0, deltaPosition, 0.0};
    }
    if (name == "Delta_y_minus") {
        return {0.0, -deltaPosition, 0.0};
    }
    if (name == "Delta_z_plus") {
        return {0.0, 0.0, deltaPosition};
    }
    if (name == "Delta_z_minus") {
        return {0.0, 0.0, -deltaPosition};
    }
    throw std::runtime_error("Unsupported PBMC electron valley name '" + name + "'.");
}

PbmcBandData evaluate_pbmc_bands(const MeshBZ&                          mesh,
                                 const uepm::PBMC::pbmc_material_model& material,
                                 double                                 deltaPosition,
                                 double                                 valenceEdge,
                                 double                                 conductionEdge,
                                 int                                    numberThreads) {
    if (material.m_electron_valleys.size() != 6 || material.m_hole_bands.size() != 2) {
        throw std::runtime_error("BandsOnBZ PBMC overlay expects six electron valleys and two hole bands.");
    }
    if (material.m_hole_bands[0].name() != "heavy_hole" || material.m_hole_bands[1].name() != "light_hole") {
        throw std::runtime_error("BandsOnBZ PBMC overlay expects heavy_hole followed by light_hole.");
    }

    PbmcBandData result;
    const auto   vertexCount = mesh.get_number_vertices();
    result.heavyHoleEnergy.resize(vertexCount);
    result.heavyHoleGradient.resize(vertexCount);
    result.lightHoleEnergy.resize(vertexCount);
    result.lightHoleGradient.resize(vertexCount);
    result.electronEnergy.resize(vertexCount);
    result.electronGradient.resize(vertexCount);
    result.electronValley.resize(vertexCount);

    std::array<vector3, 6> valleyCentersSI;
    for (std::size_t valleyIndex = 0; valleyIndex < material.m_electron_valleys.size(); ++valleyIndex) {
        valleyCentersSI[valleyIndex] = mesh.reduced_to_si_k(
            delta_valley_center_reduced(material.m_electron_valleys[valleyIndex].name(), deltaPosition));
    }

#pragma omp parallel for schedule(static) num_threads(numberThreads)
    for (std::size_t vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex) {
        const vector3 k = mesh.get_vertex_position(vertexIndex);

        for (std::size_t holeIndex = 0; holeIndex < material.m_hole_bands.size(); ++holeIndex) {
            const auto&   band       = material.m_hole_bands[holeIndex];
            const vector3 localK     = mesh.retrieve_k_inside_mesh_geometry(k);
            const auto    localPbmcK = to_pbmc_vector(localK);
            const double  energy     = valenceEdge - band.total_energy_from_k(localPbmcK);
            const vector3 gradient   = -uepm::constants::h_bar_eV * to_bz_vector(band.velocity_from_k(localPbmcK));

            if (holeIndex == 0) {
                result.heavyHoleEnergy[vertexIndex]   = energy;
                result.heavyHoleGradient[vertexIndex] = gradient;
            } else {
                result.lightHoleEnergy[vertexIndex]   = energy;
                result.lightHoleGradient[vertexIndex] = gradient;
            }
        }

        double      minimumEnergy = std::numeric_limits<double>::infinity();
        vector3     selectedGradient{};
        std::size_t selectedValley = 0;
        for (std::size_t valleyIndex = 0; valleyIndex < material.m_electron_valleys.size(); ++valleyIndex) {
            const auto&   valley     = material.m_electron_valleys[valleyIndex];
            const vector3 localK     = mesh.retrieve_k_inside_mesh_geometry(k - valleyCentersSI[valleyIndex]);
            const auto    localPbmcK = to_pbmc_vector(localK);
            const double  energy     = conductionEdge + valley.total_energy_from_k(localPbmcK);

            if (energy < minimumEnergy) {
                minimumEnergy    = energy;
                selectedGradient = uepm::constants::h_bar_eV * to_bz_vector(valley.velocity_from_k(localPbmcK));
                selectedValley   = valleyIndex;
            }
        }

        result.electronEnergy[vertexIndex]   = minimumEnergy;
        result.electronGradient[vertexIndex] = selectedGradient;
        result.electronValley[vertexIndex]   = static_cast<double>(selectedValley);
    }

    const double sampledValenceEdge =
        std::max(*std::max_element(result.heavyHoleEnergy.begin(), result.heavyHoleEnergy.end()),
                 *std::max_element(result.lightHoleEnergy.begin(), result.lightHoleEnergy.end()));
    const double sampledConductionEdge = *std::min_element(result.electronEnergy.begin(), result.electronEnergy.end());
    const double valenceShift          = valenceEdge - sampledValenceEdge;
    const double conductionShift       = conductionEdge - sampledConductionEdge;

    for (double& energy : result.heavyHoleEnergy) {
        energy += valenceShift;
    }
    for (double& energy : result.lightHoleEnergy) {
        energy += valenceShift;
    }
    for (double& energy : result.electronEnergy) {
        energy += conductionShift;
    }

    return result;
}

std::pair<double, double> epm_band_edges(const MeshBZ& mesh) {
    const auto valenceBands    = mesh.get_band_indices(MeshParticleType::valence);
    const auto conductionBands = mesh.get_band_indices(MeshParticleType::conduction);
    if (valenceBands.empty() || conductionBands.empty()) {
        throw std::runtime_error("EPM mesh must contain both valence and conduction bands for PBMC alignment.");
    }

    double valenceEdge    = std::numeric_limits<double>::lowest();
    double conductionEdge = std::numeric_limits<double>::infinity();
    for (const std::size_t band : valenceBands) {
        valenceEdge = std::max(valenceEdge, mesh.get_min_max_energy_at_band(static_cast<int>(band)).second);
    }
    for (const std::size_t band : conductionBands) {
        conductionEdge = std::min(conductionEdge, mesh.get_min_max_energy_at_band(static_cast<int>(band)).first);
    }
    return {valenceEdge, conductionEdge};
}

void export_pbmc_mesh(const std::string&                         meshFilename,
                      const uepm::pseudopotential::epm_material& epmMaterial,
                      const uepm::PBMC::pbmc_material_model&     pbmcMaterial,
                      const MeshBZ&                              epmMesh,
                      uepm::mesh_bz::BZDomainMode                domainMode,
                      double                                     deltaPosition,
                      int                                        numberThreads,
                      const std::filesystem::path&               gmshOutput,
                      const std::filesystem::path&               vtkOutput) {
    MeshBZ pbmcMesh{epmMaterial};
    pbmcMesh.set_domain_mode(domainMode);
    pbmcMesh.set_number_threads_mesh_ops(numberThreads);
    pbmcMesh.read_mesh_geometry_from_msh_file(meshFilename);

    const auto [valenceEdge, conductionEdge] = epm_band_edges(epmMesh);
    const PbmcBandData bands                 = evaluate_pbmc_bands(pbmcMesh,
                                                   pbmcMaterial,
                                                   deltaPosition,
                                                   valenceEdge,
                                                   conductionEdge,
                                                   std::max(1, numberThreads));

    pbmcMesh.append_band(MeshParticleType::valence, bands.lightHoleEnergy, bands.lightHoleGradient);
    pbmcMesh.append_band(MeshParticleType::valence, bands.heavyHoleEnergy, bands.heavyHoleGradient);
    pbmcMesh.append_band(MeshParticleType::conduction, bands.electronEnergy, bands.electronGradient);
    pbmcMesh.recompute_energies_data_and_sync(true, false, false, 0.0, 0.0);

    if (gmshOutput.has_parent_path()) {
        std::filesystem::create_directories(gmshOutput.parent_path());
    }
    if (vtkOutput.has_parent_path()) {
        std::filesystem::create_directories(vtkOutput.parent_path());
    }
    std::filesystem::remove(gmshOutput);
    pbmcMesh.export_selected_bands_to_gmsh(gmshOutput.string(), 2, 1, true, meshFilename, true);

    uepm::mesh_bz::MapStringToDoubles scalars{
        {"pbmc_hh_energy_eV", bands.heavyHoleEnergy},
        {"pbmc_lh_energy_eV", bands.lightHoleEnergy},
        {"pbmc_electron_energy_eV", bands.electronEnergy},
        {"pbmc_electron_valley", bands.electronValley},
    };
    uepm::mesh_bz::MapStringToVectors vectors{
        {"pbmc_hh_gradient_eV_m", bands.heavyHoleGradient},
        {"pbmc_lh_gradient_eV_m", bands.lightHoleGradient},
        {"pbmc_electron_gradient_eV_m", bands.electronGradient},
    };
    pbmcMesh.export_to_vtk(vtkOutput.string(), scalars, vectors, {}, {});

    fmt::print("PBMC bands aligned to E_v,max={:.9f} eV and E_c,min={:.9f} eV\n", valenceEdge, conductionEdge);
    fmt::print("Wrote PBMC Gmsh mesh to {}\n", gmshOutput.string());
    fmt::print("Wrote PBMC VTK mesh to {}\n", vtkOutput.string());
}

void write_energy_adaptive_diagnostics(const MeshBZ&                   mesh,
                                       const std::vector<std::size_t>& bands,
                                       double                          targetEnergy,
                                       std::size_t                     maxRefinementPoints,
                                       double                          backgroundH,
                                       const std::filesystem::path&    refinementMap) {
    if (bands.empty()) {
        throw std::runtime_error("No bands were selected for adaptive diagnostics.");
    }
    if (!(targetEnergy > 0.0)) {
        throw std::runtime_error("Adaptive energy target must be positive.");
    }

    const auto&                      tetrahedra = mesh.get_list_tetrahedra();
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

            const auto   energies = tetra.get_band_energies_at_vertices(band);
            const auto   minmax   = std::minmax_element(energies.begin(), energies.end());
            const double spread   = *minmax.second - *minmax.first;
            const double volume   = std::abs(tetra.get_signed_volume());
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

        const double  totalVolume = std::accumulate(samples.begin(), samples.end(), 0.0, [](double sum, const auto& s) {
            return sum + s.volume;
        });
        const auto    worst       = std::max_element(samples.begin(), samples.end(), [](const auto& a, const auto& b) {
            return a.spread < b.spread;
        });
        const Tetra&  worstTetra  = tetrahedra[worst->tetraIndex];
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
            const std::size_t count =
                static_cast<std::size_t>(std::count_if(samples.begin(), samples.end(), [threshold](const auto& s) {
                    return s.spread > threshold;
                }));
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
            const auto   energies = tetra.get_band_energies_at_vertices(band);
            const auto   minmax   = std::minmax_element(energies.begin(), energies.end());
            const double spread   = *minmax.second - *minmax.first;
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
        const double targetH  = (backgroundH > 0.0) ? std::min(currentH * scale, 0.80 * backgroundH) : currentH * scale;
        candidates.push_back({worstSpread,
                              targetH,
                              volume,
                              tetra_quality(tetra),
                              worstBand,
                              tetra.get_index(),
                              mesh.si_to_reduced_k(tetra.get_barycenter())});
    }

    AdaptiveSummary summary;
    summary.ibzTetrahedra =
        static_cast<std::size_t>(std::count_if(tetrahedra.begin(), tetrahedra.end(), [](const auto& tetra) {
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
        const double ibzVolume =
            std::accumulate(tetrahedra.begin(), tetrahedra.end(), 0.0, [](double sum, const auto& tetra) {
                return sum + (tetra.lies_in_irreducible_wedge() ? std::abs(tetra.get_signed_volume()) : 0.0);
            });

        summary.violatingVolumeFraction = violatingVolume / ibzVolume;
        summary.p50Error                = percentile(sortedErrors, 0.50);
        summary.p90Error                = percentile(sortedErrors, 0.90);
        summary.p99Error                = percentile(sortedErrors, 0.99);
        summary.maxError                = sortedErrors.back();
    }

    const std::size_t totalCandidates = candidates.size();
    candidates                        = spatially_diverse_candidates(std::move(candidates), maxRefinementPoints);

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
    TCLAP::ValueArg<std::string> arg_material("m",
                                              "material",
                                              "Symbol of the material to use (Si, Ge, GaAs, ...)",
                                              true,
                                              "Si",
                                              "string");
    TCLAP::ValueArg<std::string> arg_epm_set("d", "epm-set", "Named EPM parameter set", false, "local-cohen", "string");
    TCLAP::ValueArg<std::string> arg_epm_file("", "epm-file", "External EPM YAML parameter file", false, "", "path");
    TCLAP::ValueArg<std::string> arg_outfile("o", "outfile", "Name of the output file", false, "", "string");
    TCLAP::ValueArg<int> arg_nb_valence_bands("v", "nvbands", "Number of valence bands to export", false, 4, "int");
    TCLAP::ValueArg<int> arg_nb_conduction_bands("c",
                                                 "ncbands",
                                                 "Number of conduction bands to export",
                                                 false,
                                                 12,
                                                 "int");
    TCLAP::ValueArg<int> arg_nearest_neighbors("n",
                                               "nearestNeighbors",
                                               "number of nearest neiborgs to consider for the EPP calculation.",
                                               false,
                                               10,
                                               "int");
    TCLAP::SwitchArg     arg_enable_nonlocal_correction("C",
                                                    "nonlocal-correction",
                                                    "Enable the non-local-correction for the EPM model",
                                                    false);
    TCLAP::SwitchArg     arg_enable_soc("s", "soc", "Enable the spin-orbit coupling for the EPM model", false);
    TCLAP::SwitchArg     arg_cond_band_zero("z", "MinCondZero", "Shift the conduction band minimum to 0 eV", false);
    TCLAP::SwitchArg     arg_irr_wedge("w", "IrrWedge", "Compute bands only in the irreducible wedge of the BZ", false);
    TCLAP::ValueArg<int> arg_nb_threads("j", "nthreads", "number of threads to use.", false, 1, "int");
    TCLAP::ValueArg<std::string> arg_refinement_map("",
                                                    "refinement-map",
                                                    "Adaptive refinement CSV output (default: <mesh>_refinement.csv)",
                                                    false,
                                                    "",
                                                    "string");
    TCLAP::ValueArg<double>      arg_adaptive_energy_target("",
                                                       "adaptive-energy-target",
                                                       "Tetra energy-spread target in eV",
                                                       false,
                                                       0.020,
                                                       "float");
    TCLAP::ValueArg<int>         arg_adaptive_max_points("",
                                                 "adaptive-max-points",
                                                 "Maximum spatially diverse refinement points to export",
                                                 false,
                                                 5000,
                                                 "int");
    TCLAP::ValueArg<double>      arg_adaptive_background_h(
        "",
        "adaptive-background-h",
        "Background mesh size; exported targets are kept below this value (0 disables)",
        false,
        0.0,
        "float");
    TCLAP::SwitchArg             arg_pbmc_overlay("",
                                      "pbmc-overlay",
                                      "Export the PBMC six-valley electron and two-band hole model on the mesh",
                                      false);
    TCLAP::ValueArg<std::string> arg_pbmc_set("", "pbmc-set", "Named PBMC parameter set", false, "default", "string");
    TCLAP::ValueArg<double>      arg_delta_position("",
                                               "pbmc-delta-position",
                                               "Delta-valley position along Gamma-X in reduced coordinates",
                                               false,
                                               0.85,
                                               "float");
    TCLAP::ValueArg<std::string> arg_pbmc_outfile("",
                                                  "pbmc-outfile",
                                                  "PBMC Gmsh output (default: <mesh>_pbmc_bands.msh)",
                                                  false,
                                                  "",
                                                  "string");
    TCLAP::ValueArg<std::string> arg_pbmc_vtk("",
                                              "pbmc-vtk",
                                              "PBMC VTK output (default: <mesh>_pbmc_bands.vtk)",
                                              false,
                                              "",
                                              "string");
    TCLAP::ValueArg<std::string> arg_bz_domain("",
                                               "bz-domain",
                                               "Stored BZ domain: full or octant",
                                               false,
                                               "full",
                                               "string");
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
    cmd.add(arg_epm_set);
    cmd.add(arg_epm_file);
    cmd.add(arg_irr_wedge);
    cmd.add(arg_refinement_map);
    cmd.add(arg_adaptive_energy_target);
    cmd.add(arg_adaptive_max_points);
    cmd.add(arg_adaptive_background_h);
    cmd.add(arg_pbmc_overlay);
    cmd.add(arg_pbmc_set);
    cmd.add(arg_delta_position);
    cmd.add(arg_pbmc_outfile);
    cmd.add(arg_pbmc_vtk);
    cmd.add(arg_bz_domain);

    cmd.parse(argc, argv);

    uepm::pseudopotential::Materials         materials;
    const uepm::physics::material_repository material_repository;
    if (arg_epm_file.isSet()) {
        materials.load_material_file(material_repository, arg_material.getValue(), arg_epm_file.getValue());
        fmt::print("Loaded EPM parameter file '{}' for {}\n", arg_epm_file.getValue(), arg_material.getValue());
    } else {
        materials.load_material(material_repository, arg_material.getValue(), arg_epm_set.getValue());
        fmt::print("Loaded EPM parameter set '{}' for {}\n", arg_epm_set.getValue(), arg_material.getValue());
    }
    bool enable_nonlocal_correction = arg_enable_nonlocal_correction.isSet();
    bool enable_soc                 = arg_enable_soc.isSet();

    Options my_options;
    my_options.materialName = arg_material.getValue();

    const int max_valence_bands = enable_soc ? 8 : 4;

    my_options.nrLevels = arg_nb_valence_bands.getValue() + arg_nb_conduction_bands.getValue() +
                          max_valence_bands;  // add extra bands for the calculations (won't be exported)
    my_options.nearestNeighbors = arg_nearest_neighbors.getValue();
    my_options.nrThreads        = arg_nb_threads.getValue();
    my_options.print_options();

    uepm::pseudopotential::epm_material mat         = materials.materials.at(my_options.materialName);
    const uepm::mesh_bz::BZDomainMode   domain_mode = [&]() {
        if (arg_bz_domain.getValue() == "full") {
            return uepm::mesh_bz::BZDomainMode::full;
        }
        if (arg_bz_domain.getValue() == "octant" || arg_bz_domain.getValue() == "positive-octant") {
            return uepm::mesh_bz::BZDomainMode::positive_octant;
        }
        throw std::invalid_argument("--bz-domain must be 'full' or 'octant'");
    }();

    const std::string     mesh_filename = arg_mesh_file.getValue();
    uepm::mesh_bz::MeshBZ my_bz_mesh{mat};
    my_bz_mesh.set_domain_mode(domain_mode);
    my_bz_mesh.set_number_threads_mesh_ops(my_options.nrThreads);

    my_bz_mesh.read_mesh_geometry_from_msh_file(mesh_filename);

    auto                                 start = std::chrono::high_resolution_clock::now();
    std::vector<Vector3D<double>>        mesh_kpoints{};
    uepm::pseudopotential::BandStructure my_bandstructure;
    my_bandstructure.Initialize(mat,
                                my_options.nrLevels,
                                mesh_kpoints,
                                my_options.nearestNeighbors,
                                enable_nonlocal_correction,
                                arg_enable_soc);
    bool use_iwedge = arg_irr_wedge.isSet();
    my_bz_mesh.compute_band_structure_over_mesh(my_bandstructure, use_iwedge);

    const std::filesystem::path mesh_path(mesh_filename);
    const std::filesystem::path refinement_map =
        arg_refinement_map.isSet() ? std::filesystem::path(arg_refinement_map.getValue())
                                   : mesh_path.parent_path() / (mesh_path.stem().string() + "_refinement.csv");
    const auto diagnostic_bands =
        selected_export_bands(my_bz_mesh,
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
    std::string           out_file_bands =
        in_path.stem().replace_extension("").string() + "_" + my_bandstructure.path_band_filename();

    if (arg_outfile.isSet()) {
        out_file_bands = arg_outfile.getValue();
    }
    std::cout << "Exporting " << arg_nb_valence_bands.getValue() << " valence bands and "
              << arg_nb_conduction_bands.getValue() << " conduction bands to file: " << out_file_bands << std::endl;
    bool highest_valence_as_band0 = true;
    bool write_gradients          = true;
    my_bz_mesh.export_selected_bands_to_gmsh(out_file_bands,
                                             arg_nb_valence_bands.getValue(),
                                             arg_nb_conduction_bands.getValue(),
                                             highest_valence_as_band0,
                                             mesh_filename,
                                             write_gradients);

    const std::string vtk_file = in_path.stem().replace_extension("").string() + "_bands.vtk";
    my_bz_mesh.export_energies_and_gradients_to_vtk(vtk_file);

    if (arg_pbmc_overlay.isSet()) {
        if (!(arg_delta_position.getValue() > 0.0 && arg_delta_position.getValue() < 1.0)) {
            throw std::invalid_argument("PBMC Delta-valley position must lie strictly between Gamma and X.");
        }

        const auto pbmc_material =
            uepm::PBMC::load_pbmc_material_model(material_repository, arg_material.getValue(), arg_pbmc_set.getValue());
        const std::filesystem::path pbmc_gmsh =
            arg_pbmc_outfile.isSet() ? std::filesystem::path(arg_pbmc_outfile.getValue())
                                     : in_path.parent_path() / (in_path.stem().string() + "_pbmc_bands.msh");
        const std::filesystem::path pbmc_vtk =
            arg_pbmc_vtk.isSet() ? std::filesystem::path(arg_pbmc_vtk.getValue())
                                 : in_path.parent_path() / (in_path.stem().string() + "_pbmc_bands.vtk");

        export_pbmc_mesh(mesh_filename,
                         mat,
                         pbmc_material,
                         my_bz_mesh,
                         domain_mode,
                         arg_delta_position.getValue(),
                         my_options.nrThreads,
                         pbmc_gmsh,
                         pbmc_vtk);
    }

    return 0;
}
