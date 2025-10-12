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

#include <Eigen/Sparse>  // needed for member type in class; keep here
#include <algorithm>     // max_element
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <random>
#include <stdexcept>
#include <vector>

#include "bz_states.hpp"  // for vector3, BZ_States, material
#include "elph_common.hpp"
#include "elph_deformation_potential.hpp"
#include "elph_dispersion.hpp"
#include "elph_overlap_integral.hpp"
#include "physical_constants.hpp"

namespace uepm::mesh_bz {

// ---------- Eigen aliases ----------
using EigenSparseMatrix = Eigen::SparseMatrix<double>;
using EigenTriplet      = Eigen::Triplet<double>;

struct Rates_nk_npkp_ctor {
    PhononMode      mode      = PhononMode::none;
    PhononDirection direction = PhononDirection::none;
    PhononEvent     event     = PhononEvent::none;
    /**
     * @brief Preliminary constructor for the rate matrix (n,k) → (n',k') for a given (m,d,e) triplet.
    // (n,k) → (n',k') transition rate matrix. /!\ k' are the barycenters of final tetrahedra, not the mesh vertices!
    // So the size is (nk, n'k') where n'k' are the number of final tetrahedra × number of bands.
     *
     */
    EigenSparseMatrix matrix;  // (nk, n'k')
};

struct ItterativeBTEOptions {
    double      m_temperature     = 300.0;
    std::size_t m_max_itterations = 100;
    double      m_tolerance       = 1e-3;

    /**  Energy window around band edges to consider states for the Itterative BTE.
     *   This is used to limit the number of states considered in the Itterative BTE.
     */
    double m_energy_window_eV = 0.3;

    std::size_t              m_nb_bands = 2;      // number of bands to consider in the Itterative BTE
    std::vector<std::size_t> m_indices_vertices;  // vertices considered in the Itterative BTE
};

class ElectronPhonon : public BZ_States {
 private:
    double m_temperature           = 300.0;
    double m_rho                   = 2.329e3;
    double m_radius_wigner_seitz_m = 0.0;

    int  m_nb_threads         = 1;
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

    // (n,k) → (n',k') transition rate matrices. /!\ k' are the barycenters of final tetrahedra, not the mesh vertices!
    // So the size is (nk, n'k') where n'k' are the number of final tetrahedra × number of bands.
    // One matrix per (m,d,e) triplet.
    std::vector<Rates_nk_npkp_ctor> m_rates_nk_npkp;

 public:
    explicit ElectronPhonon(const uepm::pseudopotential::Material& material) : BZ_States(material) {}

    void   load_phonon_parameters(const std::string& filename);
    void   plot_phonon_dispersion(const std::string& filename) const;
    double get_max_phonon_energy() const;

    Rate8      compute_transition_rates_pair(int idx_n1, std::size_t idx_k1, int idx_n2, std::size_t idx_tetra_final, bool push_nk_npkp);
    RateValues compute_electron_phonon_rate(int idx_n1, std::size_t idx_k1, bool populate_nk_npkp = false);
    RateValues compute_hole_phonon_rate(int idx_n1, std::size_t idx_k1);

    void set_nb_threads(int nb_threads) {
        if (nb_threads < 1) {
            throw std::invalid_argument("nb_threads must be >= 1");
        }
        m_nb_threads = nb_threads;
    }
    void set_parallelize_over_k(bool b) noexcept { m_parallelize_over_k = b; }

    void set_temperature(double T) noexcept { m_temperature = T; }
    void set_density(double rho) noexcept { m_rho = rho; }

    void compute_electron_phonon_rates_over_mesh(double energy_max             = 100.0,
                                                 bool   irreducible_wedge_only = false,
                                                 bool   populate_nk_npkp       = false);
    void add_electron_phonon_rates_to_mesh(const std::string& initial_filename, const std::string& final_filename);
    void compute_electron_phonon_rates_over_mesh_nk_npkp(bool irreducible_wedge_only = false);

    std::pair<int, std::size_t> select_final_state(std::size_t     idx_band_initial,
                                                   const vector3&  k_initial,
                                                   PhononMode      mode,
                                                   PhononDirection direction,
                                                   PhononEvent     event,
                                                   std::mt19937&   rng) const;

    void export_rate_values(const std::string& filename) const;

    void compute_plot_electron_phonon_rates_vs_energy_over_mesh(int                nb_bands,
                                                                double             max_energy,
                                                                double             energy_step,
                                                                const std::string& filename,
                                                                bool               irreducible_wedge_only = false);

    void          read_phonon_scattering_rates_from_file(const std::filesystem::path& path);
    Rate8         interpolate_phonon_scattering_rate_at_location(const vector3& location, const std::size_t& idx_band) const;
    inline double sum_modes(const Rate8& r) const noexcept;
    double        compute_P_Gamma() const;

    void compute_RTA_mobility();
};

}  // namespace uepm::mesh_bz
