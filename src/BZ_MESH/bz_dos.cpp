/**
 * @file bz_dos.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-18
 *
 * @copyright Copyright (c) 2026
 *
 */

#include "bz_dos.hpp"

#include <fmt/core.h>

#include <algorithm>
#include <numeric>
#include <stdexcept>

#include "bz_mesh.hpp"
#include "physical_constants.hpp"

namespace uepm::mesh_bz {

double BZDos::iso_surface(const MeshBZ& mesh, double energy, int band_index) {
    double result = 0.0;
    for (const Tetra& tetra : mesh.get_list_tetrahedra()) {
        result += tetra.compute_tetra_dos_energy_band(energy, band_index);
    }
    return result;
}

double BZDos::density(const MeshBZ& mesh,
                      double        energy,
                      int           band_index,
                      bool          use_interpolation,
                      bool          irreducible_wedge_only) {
    if (mesh.stores_positive_octant() && irreducible_wedge_only) {
        throw std::invalid_argument("Cannot combine positive-octant storage with irreducible-wedge DOS");
    }
    double result = 0.0;
    for (const Tetra& tetra : mesh.get_list_tetrahedra()) {
        if (irreducible_wedge_only && !tetra.lies_in_irreducible_wedge()) {
            continue;
        }
        result += use_interpolation ? tetra.interpolate_dos_at_energy_per_band(energy, band_index)
                                    : tetra.compute_tetra_dos_energy_band(energy, band_index);
    }

    result *= mesh.get_bz_volume_correction();
    result *= mesh.get_spin_degeneracy();
    if (irreducible_wedge_only) {
        result *= uepm::constants::irreducible_wedge_factor_fcc;
    }
    return result;
}

std::vector<std::vector<double>> BZDos::curve(const MeshBZ& mesh,
                                              int           band_index,
                                              double        minimum_energy,
                                              double        maximum_energy,
                                              std::size_t   number_of_points,
                                              bool          use_interpolation,
                                              bool          irreducible_wedge_only) {
    if (number_of_points < 2) {
        throw std::invalid_argument("DOS curve requires at least two energy points");
    }

    const double        energy_step = (maximum_energy - minimum_energy) / static_cast<double>(number_of_points - 1);
    std::vector<double> energies(number_of_points);
    std::vector<double> values(number_of_points);
#pragma omp parallel for schedule(dynamic) num_threads(mesh.get_number_threads_mesh_ops())
    for (std::size_t index = 0; index < number_of_points; ++index) {
        energies[index] = minimum_energy + static_cast<double>(index) * energy_step;
        values[index]   = density(mesh, energies[index], band_index, use_interpolation, irreducible_wedge_only);
    }
    return {std::move(energies), std::move(values)};
}

std::vector<std::vector<double>> BZDos::automatic_curve(const MeshBZ& mesh,
                                                        int           band_index,
                                                        std::size_t   number_of_points,
                                                        bool          use_interpolation,
                                                        bool          irreducible_wedge_only) {
    const auto [minimum, maximum] = mesh.get_min_max_energy_at_band(band_index);
    constexpr double margin       = 0.1;
    return curve(mesh,
                 band_index,
                 minimum - margin,
                 maximum + margin,
                 number_of_points,
                 use_interpolation,
                 irreducible_wedge_only);
}

std::size_t BZStateSampler::draw_tetrahedron_at_energy(const MeshBZ& mesh,
                                                       double        energy,
                                                       std::size_t   band_index,
                                                       std::mt19937& rng) {
    std::vector<double> weights;
    weights.reserve(mesh.get_list_tetrahedra().size());
    for (const Tetra& tetra : mesh.get_list_tetrahedra()) {
        const double value = tetra.compute_tetra_dos_energy_band(energy, band_index);
        weights.push_back(std::isfinite(value) ? std::max(0.0, value) : 0.0);
    }
    if (std::accumulate(weights.begin(), weights.end(), 0.0) <= 0.0) {
        throw std::runtime_error("Total DOS is zero at the requested energy");
    }
    std::discrete_distribution<std::size_t> distribution(weights.begin(), weights.end());
    return distribution(rng);
}

vector3 BZStateSampler::draw_k_at_energy(const MeshBZ& mesh, double energy, std::size_t band_index, std::mt19937& rng) {
    const auto [minimum, maximum] = mesh.get_min_max_energy_at_band(static_cast<int>(band_index));
    if (energy < minimum || energy > maximum) {
        throw std::runtime_error(fmt::format("Energy {:.6f} eV is outside band {} [{:.6f}, {:.6f}] eV",
                                             energy,
                                             band_index,
                                             minimum,
                                             maximum));
    }
    const std::size_t tetra_index = draw_tetrahedron_at_energy(mesh, energy, band_index, rng);
    const vector3     representative =
        mesh.get_list_tetrahedra().at(tetra_index).draw_random_uniform_point_at_energy(energy, band_index, rng);
    if (!mesh.stores_positive_octant()) {
        return representative;
    }
    std::uniform_int_distribution<std::size_t> image_distribution(0, positive_octant_images.size() - 1);
    return apply_sign_image(representative, positive_octant_images[image_distribution(rng)]);
}

std::pair<vector3, std::size_t> BZStateSampler::draw_k_at_energy(const MeshBZ& mesh, double energy, std::mt19937& rng) {
    std::vector<std::size_t> candidate_bands;
    std::vector<double>      weights;
    for (std::size_t band_index = 0; band_index < mesh.get_number_bands_total(); ++band_index) {
        const auto [minimum, maximum] = mesh.get_min_max_energy_at_band(static_cast<int>(band_index));
        if (energy >= minimum && energy <= maximum) {
            candidate_bands.push_back(band_index);
            weights.push_back(std::max(0.0, BZDos::density(mesh, energy, static_cast<int>(band_index))));
        }
    }
    if (candidate_bands.empty() || std::accumulate(weights.begin(), weights.end(), 0.0) <= 0.0) {
        throw std::runtime_error("No band has non-zero DOS at the requested energy");
    }

    std::discrete_distribution<std::size_t> distribution(weights.begin(), weights.end());
    const std::size_t                       band_index = candidate_bands[distribution(rng)];
    return {draw_k_at_energy(mesh, energy, band_index, rng), band_index};
}

double MeshBZ::compute_iso_surface(double iso_energy, int band_index) const {
    return BZDos::iso_surface(*this, iso_energy, band_index);
}

double MeshBZ::compute_dos_at_energy_and_band(double iso_energy, int band_index, bool use_interp, bool use_iw) const {
    return BZDos::density(*this, iso_energy, band_index, use_interp, use_iw);
}

std::vector<std::vector<double>> MeshBZ::compute_dos_band_at_band(int         band_index,
                                                                  double      min_energy,
                                                                  double      max_energy,
                                                                  std::size_t nb_points,
                                                                  bool        use_interp,
                                                                  bool        use_iw) const {
    return BZDos::curve(*this, band_index, min_energy, max_energy, nb_points, use_interp, use_iw);
}

std::vector<std::vector<double>> MeshBZ::compute_dos_band_at_band_auto(int         band_index,
                                                                       std::size_t nb_points,
                                                                       bool        use_interp,
                                                                       bool        use_iw) const {
    return BZDos::automatic_curve(*this, band_index, nb_points, use_interp, use_iw);
}

std::size_t MeshBZ::draw_random_tetrahedron_index_with_dos_probability(double        energy,
                                                                       std::size_t   idx_band,
                                                                       std::mt19937& random_generator) const {
    return BZStateSampler::draw_tetrahedron_at_energy(*this, energy, idx_band, random_generator);
}

vector3 MeshBZ::draw_random_k_point_at_energy(double        energy,
                                              std::size_t   idx_band,
                                              std::mt19937& random_generator) const {
    return BZStateSampler::draw_k_at_energy(*this, energy, idx_band, random_generator);
}

std::pair<vector3, std::size_t> MeshBZ::draw_random_k_point_at_energy(double energy, std::mt19937& rng) const {
    return BZStateSampler::draw_k_at_energy(*this, energy, rng);
}

}  // namespace uepm::mesh_bz
