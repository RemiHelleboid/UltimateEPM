/**
 * @file bz_dos.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2026-06-18
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include "vector_bz.hpp"

namespace uepm::mesh_bz {

class MeshBZ;

class BZDos {
 public:
    static double iso_surface(const MeshBZ& mesh, double energy, int band_index);
    static double density(const MeshBZ& mesh,
                          double        energy,
                          int           band_index,
                          bool          use_interpolation      = false,
                          bool          irreducible_wedge_only = false);

    static std::vector<std::vector<double>> curve(const MeshBZ& mesh,
                                                  int           band_index,
                                                  double        minimum_energy,
                                                  double        maximum_energy,
                                                  std::size_t   number_of_points,
                                                  bool          use_interpolation      = false,
                                                  bool          irreducible_wedge_only = false);

    static std::vector<std::vector<double>> automatic_curve(const MeshBZ& mesh,
                                                            int           band_index,
                                                            std::size_t   number_of_points,
                                                            bool          use_interpolation      = false,
                                                            bool          irreducible_wedge_only = false);
};

class BZStateSampler {
 public:
    static std::size_t draw_tetrahedron_at_energy(const MeshBZ& mesh,
                                                  double        energy,
                                                  std::size_t   band_index,
                                                  std::mt19937& rng);

    static vector3 draw_k_at_energy(const MeshBZ& mesh, double energy, std::size_t band_index, std::mt19937& rng);

    static std::pair<vector3, std::size_t> draw_k_at_energy(const MeshBZ& mesh, double energy, std::mt19937& rng);
};

}  // namespace uepm::mesh_bz
