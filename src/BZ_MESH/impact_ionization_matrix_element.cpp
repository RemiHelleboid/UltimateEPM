/**
 * @file impact_ionization_matrix_element.cpp
 * @brief Plane-wave two-body matrix elements for ab initio impact ionization.
 */

#include "impact_ionization_matrix_element.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <unordered_map>

namespace uepm::mesh_bz {

namespace {

struct MomentumKey {
    long long x = 0;
    long long y = 0;
    long long z = 0;

    bool operator==(const MomentumKey& other) const noexcept { return x == other.x && y == other.y && z == other.z; }

    MomentumKey opposite() const noexcept { return MomentumKey{.x = -x, .y = -y, .z = -z}; }
};

struct MomentumKeyHash {
    std::size_t operator()(const MomentumKey& key) const noexcept {
        std::size_t seed = 0;
        const auto  mix  = [&seed](long long value) {
            const std::size_t hashed = std::hash<long long>{}(value);
            seed ^= hashed + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
        };
        mix(key.x);
        mix(key.y);
        mix(key.z);
        return seed;
    }
};

struct TransitionDensity {
    vector3   q_SI{};
    complex_d amplitude{0.0, 0.0};
};

using TransitionDensityMap = std::unordered_map<MomentumKey, TransitionDensity, MomentumKeyHash>;

void require_finite_vector(const vector3& vector, const char* name) {
    if (!std::isfinite(vector.x()) || !std::isfinite(vector.y()) || !std::isfinite(vector.z())) {
        throw std::invalid_argument(std::string(name) + " must be finite.");
    }
}

void validate_state(const ImpactIonizationPlaneWaveState& state, std::size_t basis_size, const char* name) {
    require_finite_vector(state.k_SI, name);
    if (!std::isfinite(state.energy_eV)) {
        throw std::invalid_argument(std::string(name) + " energy must be finite.");
    }
    if (static_cast<std::size_t>(state.coefficients.size()) != basis_size) {
        throw std::invalid_argument(std::string(name) + " coefficient count does not match basis size.");
    }
    for (Eigen::Index idx = 0; idx < state.coefficients.size(); ++idx) {
        const auto value = state.coefficients[idx];
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
            throw std::invalid_argument(std::string(name) + " coefficients must be finite.");
        }
    }
}

MomentumKey quantize_momentum(const vector3& q_SI, double tolerance_SI) {
    return MomentumKey{.x = std::llround(q_SI.x() / tolerance_SI),
                       .y = std::llround(q_SI.y() / tolerance_SI),
                       .z = std::llround(q_SI.z() / tolerance_SI)};
}

TransitionDensityMap build_transition_density(const ImpactIonizationPlaneWaveState& initial,
                                              const ImpactIonizationPlaneWaveState& final,
                                              const std::vector<vector3>&           basis_vectors_SI,
                                              double                                momentum_tolerance_SI) {
    TransitionDensityMap density;
    const std::size_t    basis_size = basis_vectors_SI.size();
    density.reserve(basis_size * basis_size);

    for (std::size_t idx_initial_g = 0; idx_initial_g < basis_size; ++idx_initial_g) {
        const auto initial_amplitude = initial.coefficients[static_cast<Eigen::Index>(idx_initial_g)];
        if (std::norm(initial_amplitude) == 0.0) {
            continue;
        }
        for (std::size_t idx_final_g = 0; idx_final_g < basis_size; ++idx_final_g) {
            const auto final_amplitude = final.coefficients[static_cast<Eigen::Index>(idx_final_g)];
            if (std::norm(final_amplitude) == 0.0) {
                continue;
            }

            const vector3 q_SI =
                initial.k_SI + basis_vectors_SI[idx_initial_g] - final.k_SI - basis_vectors_SI[idx_final_g];
            const MomentumKey key   = quantize_momentum(q_SI, momentum_tolerance_SI);
            auto&             entry = density[key];
            if (entry.amplitude == complex_d{0.0, 0.0}) {
                entry.q_SI = q_SI;
            }
            entry.amplitude += std::conj(final_amplitude) * initial_amplitude;
        }
    }

    return density;
}

void validate_config(const ImpactIonizationMatrixElementConfig& config) {
    if (!(config.momentum_tolerance_SI > 0.0) || !std::isfinite(config.momentum_tolerance_SI)) {
        throw std::invalid_argument("Impact-ionization momentum tolerance must be positive and finite.");
    }
    if (!(config.normalization_volume_m3 > 0.0) || !std::isfinite(config.normalization_volume_m3)) {
        throw std::invalid_argument("Impact-ionization normalization volume must be positive and finite.");
    }
}

}  // namespace

complex_d compute_screened_two_body_matrix_element(const ImpactIonizationPlaneWaveState&      initial_a,
                                                   const ImpactIonizationPlaneWaveState&      initial_b,
                                                   const ImpactIonizationPlaneWaveState&      final_a,
                                                   const ImpactIonizationPlaneWaveState&      final_b,
                                                   const std::vector<vector3>&                basis_vectors_SI,
                                                   const impact_screened_interaction&         screened_interaction,
                                                   const ImpactIonizationMatrixElementConfig& config) {
    validate_config(config);
    if (!screened_interaction) {
        throw std::invalid_argument("Impact-ionization screened interaction callback is empty.");
    }
    if (basis_vectors_SI.empty()) {
        throw std::invalid_argument("Impact-ionization basis vector list is empty.");
    }
    for (const auto& basis_vector : basis_vectors_SI) {
        require_finite_vector(basis_vector, "basis vector");
    }

    const std::size_t basis_size = basis_vectors_SI.size();
    validate_state(initial_a, basis_size, "initial_a");
    validate_state(initial_b, basis_size, "initial_b");
    validate_state(final_a, basis_size, "final_a");
    validate_state(final_b, basis_size, "final_b");

    const auto density_a = build_transition_density(initial_a, final_a, basis_vectors_SI, config.momentum_tolerance_SI);
    const auto density_b = build_transition_density(initial_b, final_b, basis_vectors_SI, config.momentum_tolerance_SI);
    const double energy_transfer_eV = std::abs(initial_a.energy_eV - final_a.energy_eV);

    complex_d matrix_element{0.0, 0.0};
    for (const auto& [key, density_entry_a] : density_a) {
        const auto density_b_it = density_b.find(key.opposite());
        if (density_b_it == density_b.end()) {
            continue;
        }
        const complex_d interaction = screened_interaction(density_entry_a.q_SI, energy_transfer_eV);
        if (!std::isfinite(interaction.real()) || !std::isfinite(interaction.imag())) {
            throw std::runtime_error("Impact-ionization screened interaction returned a non-finite value.");
        }
        matrix_element += density_entry_a.amplitude * density_b_it->second.amplitude * interaction;
    }

    return matrix_element / config.normalization_volume_m3;
}

std::array<complex_d, 2> compute_direct_exchange_impact_ionization_matrix_element(
    const ImpactIonizationPlaneWaveState&      initial_hot_electron,
    const ImpactIonizationPlaneWaveState&      initial_valence_electron,
    const ImpactIonizationPlaneWaveState&      final_electron_1,
    const ImpactIonizationPlaneWaveState&      final_electron_2,
    const std::vector<vector3>&                basis_vectors_SI,
    const impact_screened_interaction&         screened_interaction,
    const ImpactIonizationMatrixElementConfig& config) {
    const complex_d direct   = compute_screened_two_body_matrix_element(initial_hot_electron,
                                                                      initial_valence_electron,
                                                                      final_electron_1,
                                                                      final_electron_2,
                                                                      basis_vectors_SI,
                                                                      screened_interaction,
                                                                      config);
    const complex_d exchange = compute_screened_two_body_matrix_element(initial_hot_electron,
                                                                        initial_valence_electron,
                                                                        final_electron_2,
                                                                        final_electron_1,
                                                                        basis_vectors_SI,
                                                                        screened_interaction,
                                                                        config);
    return {direct, exchange};
}

double antisymmetrized_impact_ionization_strength(complex_d direct, complex_d exchange) {
    return std::norm(direct - exchange);
}

}  // namespace uepm::mesh_bz
