/**
 * @file impact_ionization_state.hpp
 * @brief EPM state access layer for ab initio impact ionization.
 */

#pragma once

#include <Eigen/Dense>
#include <cstddef>

#include "bz_states.hpp"
#include "vector_bz.hpp"

namespace uepm::mesh_bz {

struct ImpactIonizationBandState {
    std::size_t k_index{};
    std::size_t band_index{};
    vector3     k_SI{};
    double      energy_eV{};
};

/**
 * @brief Read-only facade exposing the EPM quantities needed by impact ionization.
 */
class ImpactIonizationStateTable {
 private:
    const BZ_States* m_states = nullptr;

 public:
    explicit ImpactIonizationStateTable(const BZ_States& states);

    std::size_t number_kpoints() const noexcept;
    std::size_t number_stored_bands() const;

    ImpactIonizationBandState state(std::size_t idx_k, std::size_t idx_band) const;
    const Eigen::VectorXd&    eigenvalues_at_k(std::size_t idx_k) const;
    const Eigen::MatrixXcd&   eigenvectors_at_k(std::size_t idx_k) const;
    Eigen::VectorXcd          coefficients(std::size_t idx_k, std::size_t idx_band) const;
};

}  // namespace uepm::mesh_bz
