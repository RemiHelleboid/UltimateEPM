/**
 * @file impact_ionization_state.cpp
 * @brief EPM state access layer for ab initio impact ionization.
 */

#include "impact_ionization_state.hpp"

#include <stdexcept>

namespace uepm::mesh_bz {

ImpactIonizationStateTable::ImpactIonizationStateTable(const BZ_States& states) : m_states(&states) {
    if (!states.has_eigenstates()) {
        throw std::logic_error("ImpactIonizationStateTable requires computed BZ eigenstates.");
    }
}

std::size_t ImpactIonizationStateTable::number_kpoints() const noexcept { return m_states->get_number_vertices(); }

std::size_t ImpactIonizationStateTable::number_stored_bands() const {
    if (number_kpoints() == 0) {
        return 0;
    }
    return static_cast<std::size_t>(m_states->get_eigenvectors_at_k(0).cols());
}

ImpactIonizationBandState ImpactIonizationStateTable::state(std::size_t idx_k, std::size_t idx_band) const {
    return ImpactIonizationBandState{.k_index    = idx_k,
                                     .band_index = idx_band,
                                     .k_SI       = m_states->get_vertex_position(idx_k),
                                     .energy_eV  = m_states->get_eigenvalue_eV(idx_k, idx_band)};
}

const Eigen::VectorXd& ImpactIonizationStateTable::eigenvalues_at_k(std::size_t idx_k) const {
    return m_states->get_eigenvalues_at_k(idx_k);
}

const Eigen::MatrixXcd& ImpactIonizationStateTable::eigenvectors_at_k(std::size_t idx_k) const {
    return m_states->get_eigenvectors_at_k(idx_k);
}

Eigen::VectorXcd ImpactIonizationStateTable::coefficients(std::size_t idx_k, std::size_t idx_band) const {
    return m_states->get_eigenvectors_at_k(idx_k).col(static_cast<Eigen::Index>(idx_band));
}

}  // namespace uepm::mesh_bz
