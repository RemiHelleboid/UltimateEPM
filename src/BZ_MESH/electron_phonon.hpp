/**
 * @file electron_phonon.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2025-10-04
 *
 *
 */

#pragma once

#include <fmt/core.h>
#include <fmt/format.h>

#include <Eigen/Sparse>  // needed for member type in class; keep here
#include <algorithm>     // max_element
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <random>
#include <stdexcept>
#include <vector>

#include "bz_states.hpp"
#include "elph_common.hpp"
#include "elph_deformation_potential.hpp"
#include "elph_dispersion.hpp"
#include "elph_overlap_integral.hpp"
#include "physical_constants.hpp"

namespace uepm::mesh_bz {

// ---------- Eigen aliases ----------
using EigenSparseMatrix = Eigen::SparseMatrix<double>;
using EigenTriplet      = Eigen::Triplet<double>;

struct PGamma {
    std::vector<double> m_energies_eV;
    std::vector<double> m_P_Gamma_values;

    double interpolate_P_Gamma(double energy_eV) const;
};

struct SelectedFinalState {
    std::size_t           idx_final_band;
    std::size_t           idx_final_tetra;
    uepm::mesh_bz::Tetra* ptr_final_tetra;
    vector3               k_final;     // Randomly selected final k-point on the final energy surface
    double                E_final_eV;  // Energy of the final state
};

struct IterativeBTEOptions {
    double      m_temperature_K  = 300.0;
    std::size_t m_max_iterations = 100;
    double      m_tolerance      = 1e-3;

    /**  Energy window around band edges to consider states for the Iterative BTE.
     *   This is used to limit the number of states considered in the Iterative BTE.
     */
    double m_energy_window_eV = 0.3;

    std::size_t              m_nb_bands = 2;      // number of bands to consider in the Iterative BTE
    std::vector<std::size_t> m_indices_vertices;  // vertices considered in the Iterative BTE
};

class ElectronPhonon : public BZ_States {
 private:
    double m_temperature_K         = 300.0;
    double m_rho_kg_m3             = 2.329e3;
    double m_radius_wigner_seitz_m = 0.0;

    MeshParticleType m_elph_particle_type = MeshParticleType::conduction;
    std::size_t      m_nb_bands_elph      = 0;

    bool m_parallelize_over_k = true;

    HoleOverlapIntParams m_hole_overlap_int_params;
    DeformationPotential m_ac_defpot_e, m_op_defpot_e;
    DeformationPotential m_ac_defpot_h, m_op_defpot_h;

    // 4 dispersions: (ac/op) × (L/T); index via md_index()
    std::array<PhononDispersion, 4> m_phonon_dispersion;

    // Transport rates. m_phonon_rates_transport[band][k1]
    std::vector<std::vector<double>> m_phonon_rates_transport;  // precomputed 1/τ_tr(E) on uniform grid

    // Precomputed rates on mesh: [vertex][band][8]
    std::vector<std::vector<Rate8>> m_list_phonon_scattering_rates;
    std::vector<double>             m_count_weight_tetra_per_vertex;

    PGamma m_P_Gamma_data;

    const double m_fit_optical  = 1.0;
    const double m_fit_acoustic = 1.0;

 public:
    explicit ElectronPhonon(const uepm::pseudopotential::Material& material) : BZ_States(material) {}

    void   load_phonon_parameters(const std::string& filename);
    void   plot_phonon_dispersion(const std::string& filename) const;
    double get_max_phonon_energy() const;
    void   export_phonon_dispersion(const std::string& filename) const;

    void             set_parallelize_over_k(bool b) noexcept { m_parallelize_over_k = b; }
    void             set_temperature(double T) noexcept { m_temperature_K = T; }
    void             set_density(double rho) noexcept { m_rho_kg_m3 = rho; }
    void             set_particle_type(MeshParticleType type) noexcept { m_elph_particle_type = type; }
    MeshParticleType get_particle_type() const noexcept { return m_elph_particle_type; }
    void             set_nb_bands_elph(std::size_t nb) noexcept { m_nb_bands_elph = nb; }
    std::size_t      get_nb_bands_elph() const noexcept { return m_nb_bands_elph; }

    Rate8 compute_electron_phonon_transition_rates_pair(std::size_t idx_n1,
                                                        std::size_t idx_k1,
                                                        std::size_t idx_n2,
                                                        std::size_t idx_tetra_final);

    RateValues compute_electron_phonon_rate(std::size_t idx_n1, std::size_t idx_k1);
    RateValues compute_hole_phonon_rate(std::size_t idx_n1, std::size_t idx_k1);

    double scale_q_norm(double q_norm) const;
    void   compute_electron_phonon_rates_over_mesh(double energy_max = 100.0, bool irreducible_wedge_only = false);
    void   add_electron_phonon_rates_to_mesh(const std::string& initial_filename, const std::string& final_filename);
    void   clean_all_elph_data();

    SelectedFinalState select_electron_phonon_final_state(std::size_t     idx_band_initial,
                                                          const vector3&  k_initial,
                                                          PhononMode      mode,
                                                          PhononDirection direction,
                                                          PhononEvent     event,
                                                          std::mt19937&   rng) const;
    SelectedFinalState select_electron_phonon_final_state(std::size_t    idx_band_initial,
                                                          const vector3& k_initial,
                                                          int            idx_phonon_branch,
                                                          std::mt19937&  rng) const;

    void export_rate_values(const std::string& filename) const;
    void compute_plot_electron_phonon_rates_vs_energy_over_mesh(double max_energy, double energy_step, const std::string& filename);

    void          read_phonon_scattering_rates_from_file(const std::filesystem::path& path);
    Rate8         interpolate_phonon_scattering_rate_at_location(const vector3& location, const std::size_t& idx_band) const;
    inline double sum_modes(const Rate8& r) const noexcept;
    double        compute_P_Gamma() const;

    Eigen::Matrix3d compute_electron_MRTA_mobility_tensor(double fermi_level_eV, double temperature_K, bool conduction_only = true);
    double          compute_electron_MRTA_mobility_isotropic(double fermi_level_eV, double temperature_K, bool conduction_only = true);
    double          mean_electron_energy_equilibrium(double fermi_level_eV, double temperature_K, bool excess_above_cbm = false) const;
    void            test_elph() const;
};

}  // namespace uepm::mesh_bz
