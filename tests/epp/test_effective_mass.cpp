#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include <doctest/doctest.h>

#include "Vector3D.h"
#include "effective_mass.hpp"
#include "physical_constants.hpp"

namespace {

double curvature_from_mass(double mass_m0) {
    return (uepm::constants::h_bar * uepm::constants::h_bar) /
           (uepm::constants::eV_to_J * uepm::constants::m_e * mass_m0);
}

std::vector<uepm::pseudopotential::valley_fit_sample> make_synthetic_kane_samples(
    const Vector3D<double>&       k0,
    double                        lattice_constant_m,
    const std::array<double, 3>&  masses_m0,
    double                        alpha_eV_inv,
    const std::vector<double>&    radii,
    uepm::pseudopotential::band_edge_kind edge_kind) {
    const double c0          = curvature_from_mass(masses_m0[0]);
    const double c1          = curvature_from_mass(masses_m0[1]);
    const double c2          = curvature_from_mass(masses_m0[2]);
    const double reduced_to_k = 2.0 * uepm::constants::pi / lattice_constant_m;
    const double sign         = edge_kind == uepm::pseudopotential::band_edge_kind::minimum ? 1.0 : -1.0;

    std::vector<uepm::pseudopotential::valley_fit_sample> samples;
    samples.push_back({k0, 1.25});
    for (double h : radii) {
        for (int ix = -1; ix <= 1; ++ix) {
            for (int iy = -1; iy <= 1; ++iy) {
                for (int iz = -1; iz <= 1; ++iz) {
                    if (ix == 0 && iy == 0 && iz == 0) {
                        continue;
                    }
                    const Vector3D<double> k(k0.X + h * static_cast<double>(ix),
                                             k0.Y + h * static_cast<double>(iy),
                                             k0.Z + h * static_cast<double>(iz));
                    const double qx    = (k.X - k0.X) * reduced_to_k;
                    const double qy    = (k.Y - k0.Y) * reduced_to_k;
                    const double qz    = (k.Z - k0.Z) * reduced_to_k;
                    const double gamma = 0.5 * (c0 * qx * qx + c1 * qy * qy + c2 * qz * qz);
                    const double energy =
                        alpha_eV_inv > 0.0 ? (std::sqrt(1.0 + 4.0 * alpha_eV_inv * gamma) - 1.0) /
                                                  (2.0 * alpha_eV_inv)
                                            : gamma;
                    samples.push_back({k, 1.25 + sign * energy});
                }
            }
        }
    }
    return samples;
}

std::vector<uepm::pseudopotential::valley_fit_sample> make_synthetic_kane_samples(
    const Vector3D<double>&       k0,
    double                        lattice_constant_m,
    const std::array<double, 3>&  masses_m0,
    double                        alpha_eV_inv,
    uepm::pseudopotential::band_edge_kind edge_kind) {
    return make_synthetic_kane_samples(k0, lattice_constant_m, masses_m0, alpha_eV_inv, {0.006, 0.012, 0.018}, edge_kind);
}

}  // namespace

TEST_CASE("effective-mass fit recovers anisotropic Kane parameters around a minimum") {
    const Vector3D<double>      k0(0.85, 0.0, 0.0);
    const double                lattice_constant_m = 5.431e-10;
    const std::array<double, 3> masses{0.19, 0.32, 0.91};
    const double                alpha = 0.45;
    const auto                  samples =
        make_synthetic_kane_samples(k0, lattice_constant_m, masses, alpha, uepm::pseudopotential::band_edge_kind::minimum);

    const auto result = uepm::pseudopotential::fit_effective_mass_and_nonparabolicity(
        samples, k0, 1.25, lattice_constant_m, uepm::pseudopotential::band_edge_kind::minimum);

    auto fitted_masses = result.principal_masses_m0;
    auto expected      = masses;
    std::sort(fitted_masses.begin(), fitted_masses.end());
    std::sort(expected.begin(), expected.end());

    CHECK(fitted_masses[0] == doctest::Approx(expected[0]).epsilon(2.0e-2));
    CHECK(fitted_masses[1] == doctest::Approx(expected[1]).epsilon(2.0e-2));
    CHECK(fitted_masses[2] == doctest::Approx(expected[2]).epsilon(2.0e-2));
    CHECK(result.non_parabolicity_eV_inv == doctest::Approx(alpha).epsilon(2.0e-2));
    CHECK(result.rms_error_meV < 0.1);
}

TEST_CASE("effective-mass fit supports valence-like maxima") {
    const Vector3D<double>      k0(0.0, 0.0, 0.0);
    const double                lattice_constant_m = 5.431e-10;
    const std::array<double, 3> masses{0.23, 0.23, 0.49};
    const double                alpha = 0.2;
    const auto                  samples =
        make_synthetic_kane_samples(k0, lattice_constant_m, masses, alpha, uepm::pseudopotential::band_edge_kind::maximum);

    const auto result = uepm::pseudopotential::fit_effective_mass_and_nonparabolicity(
        samples, k0, 1.25, lattice_constant_m, uepm::pseudopotential::band_edge_kind::maximum);

    auto fitted_masses = result.principal_masses_m0;
    auto expected      = masses;
    std::sort(fitted_masses.begin(), fitted_masses.end());
    std::sort(expected.begin(), expected.end());

    CHECK(fitted_masses[0] == doctest::Approx(expected[0]).epsilon(2.0e-2));
    CHECK(fitted_masses[1] == doctest::Approx(expected[1]).epsilon(2.0e-2));
    CHECK(fitted_masses[2] == doctest::Approx(expected[2]).epsilon(2.0e-2));
    CHECK(result.non_parabolicity_eV_inv == doctest::Approx(alpha).epsilon(2.0e-2));
    CHECK(result.rms_error_meV < 0.1);
}

TEST_CASE("split fit keeps near-edge masses while fitting alpha over a wider energy window") {
    const Vector3D<double>      k0(0.85, 0.0, 0.0);
    const double                lattice_constant_m = 5.431e-10;
    const std::array<double, 3> masses{0.19, 0.32, 0.91};
    const double                alpha = 0.55;
    const auto                  mass_samples = make_synthetic_kane_samples(
        k0, lattice_constant_m, masses, alpha, {0.001, 0.002, 0.003}, uepm::pseudopotential::band_edge_kind::minimum);
    const auto alpha_samples = make_synthetic_kane_samples(
        k0, lattice_constant_m, masses, alpha, {0.006, 0.012, 0.018}, uepm::pseudopotential::band_edge_kind::minimum);

    const auto result = uepm::pseudopotential::fit_effective_mass_then_nonparabolicity(
        mass_samples,
        alpha_samples,
        k0,
        1.25,
        lattice_constant_m,
        uepm::pseudopotential::band_edge_kind::minimum,
        0.4,
        true);

    auto fitted_masses = result.principal_masses_m0;
    auto expected      = masses;
    std::sort(fitted_masses.begin(), fitted_masses.end());
    std::sort(expected.begin(), expected.end());

    CHECK(fitted_masses[0] == doctest::Approx(expected[0]).epsilon(2.0e-2));
    CHECK(fitted_masses[1] == doctest::Approx(expected[1]).epsilon(2.0e-2));
    CHECK(fitted_masses[2] == doctest::Approx(expected[2]).epsilon(2.0e-2));
    CHECK(result.non_parabolicity_eV_inv == doctest::Approx(alpha).epsilon(2.0e-2));
    CHECK(result.alpha_sample_count > 0);
}

TEST_CASE("fixed-mass alpha fit can clamp unphysical negative alpha") {
    const Vector3D<double>      k0(0.85, 0.0, 0.0);
    const double                lattice_constant_m = 5.431e-10;
    const std::array<double, 3> masses{0.19, 0.32, 0.91};
    const auto                  samples =
        make_synthetic_kane_samples(k0, lattice_constant_m, masses, 0.0, uepm::pseudopotential::band_edge_kind::minimum);

    auto mass_fit = uepm::pseudopotential::fit_effective_mass_tensor(
        samples, k0, 1.25, lattice_constant_m, uepm::pseudopotential::band_edge_kind::minimum);
    mass_fit.hessian_eV_m2[0][0] *= 0.5;
    mass_fit.hessian_eV_m2[1][1] *= 0.5;
    mass_fit.hessian_eV_m2[2][2] *= 0.5;

    const auto unclamped =
        uepm::pseudopotential::fit_nonparabolicity_with_fixed_mass(samples, mass_fit, lattice_constant_m, 0.4, false);
    const auto clamped =
        uepm::pseudopotential::fit_nonparabolicity_with_fixed_mass(samples, mass_fit, lattice_constant_m, 0.4, true);

    CHECK(unclamped.non_parabolicity_eV_inv < 0.0);
    CHECK(clamped.non_parabolicity_eV_inv == doctest::Approx(0.0));
}
