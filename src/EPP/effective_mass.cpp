#include "effective_mass.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "physical_constants.hpp"

namespace uepm::pseudopotential {
namespace {

constexpr double kinetic_tolerance_eV = 1.0e-12;

double edge_sign(band_edge_kind edge_kind) { return edge_kind == band_edge_kind::minimum ? 1.0 : -1.0; }

Eigen::Vector3d to_eigen(const Vector3D<double>& value) { return Eigen::Vector3d{value.X, value.Y, value.Z}; }

std::array<std::array<double, 3>, 3> to_array(const Eigen::Matrix3d& matrix) {
    std::array<std::array<double, 3>, 3> result{};
    for (Eigen::Index row = 0; row < 3; ++row) {
        for (Eigen::Index col = 0; col < 3; ++col) {
            result[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)] = matrix(row, col);
        }
    }
    return result;
}

Eigen::Matrix3d to_matrix(const std::array<std::array<double, 3>, 3>& array) {
    Eigen::Matrix3d matrix;
    for (Eigen::Index row = 0; row < 3; ++row) {
        for (Eigen::Index col = 0; col < 3; ++col) {
            matrix(row, col) = array[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)];
        }
    }
    return matrix;
}

Eigen::MatrixXd build_quadratic_design(const std::vector<valley_fit_sample>& samples,
                                       const Vector3D<double>&               k0_reduced,
                                       double                                lattice_constant_m) {
    const double          reduced_to_k_m = 2.0 * uepm::constants::pi / lattice_constant_m;
    const Eigen::Vector3d k0             = to_eigen(k0_reduced);
    Eigen::MatrixXd       design(samples.size(), 6);

    for (std::size_t i = 0; i < samples.size(); ++i) {
        const Eigen::Vector3d dk_reduced = to_eigen(samples[i].k_reduced) - k0;
        const Eigen::Vector3d q          = dk_reduced * reduced_to_k_m;
        const double          qx         = q.x();
        const double          qy         = q.y();
        const double          qz         = q.z();

        design(static_cast<Eigen::Index>(i), 0) = 0.5 * qx * qx;
        design(static_cast<Eigen::Index>(i), 1) = 0.5 * qy * qy;
        design(static_cast<Eigen::Index>(i), 2) = 0.5 * qz * qz;
        design(static_cast<Eigen::Index>(i), 3) = qx * qy;
        design(static_cast<Eigen::Index>(i), 4) = qx * qz;
        design(static_cast<Eigen::Index>(i), 5) = qy * qz;
    }
    return design;
}

Eigen::VectorXd kinetic_energies(const std::vector<valley_fit_sample>& samples,
                                 double                                edge_energy_eV,
                                 band_edge_kind                        edge_kind) {
    const double    sign = edge_sign(edge_kind);
    Eigen::VectorXd kinetic(samples.size());
    for (std::size_t i = 0; i < samples.size(); ++i) {
        kinetic(static_cast<Eigen::Index>(i)) = sign * (samples[i].energy_eV - edge_energy_eV);
    }
    return kinetic;
}

Eigen::Matrix3d fit_hessian_from_gamma_target(const Eigen::MatrixXd& design, const Eigen::VectorXd& gamma_target) {
    const Eigen::VectorXd coeff = design.colPivHouseholderQr().solve(gamma_target);
    Eigen::Matrix3d       hessian;
    hessian << coeff(0), coeff(3), coeff(4), coeff(3), coeff(1), coeff(5), coeff(4), coeff(5), coeff(2);
    return hessian;
}

std::array<double, 3> masses_from_hessian(const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>& solver) {
    std::array<double, 3> masses{};
    for (Eigen::Index i = 0; i < 3; ++i) {
        const double curvature = solver.eigenvalues()(i);
        if (curvature > 0.0 && std::isfinite(curvature)) {
            const double mass_kg =
                (uepm::constants::h_bar * uepm::constants::h_bar) / (uepm::constants::eV_to_J * curvature);
            masses[static_cast<std::size_t>(i)] = mass_kg / uepm::constants::m_e;
        } else {
            masses[static_cast<std::size_t>(i)] = std::numeric_limits<double>::quiet_NaN();
        }
    }
    return masses;
}

}  // namespace

effective_mass_fit_result fit_effective_mass_and_nonparabolicity(const std::vector<valley_fit_sample>& samples,
                                                                 const Vector3D<double>&               k0_reduced,
                                                                 double                                edge_energy_eV,
                                                                 double         lattice_constant_m,
                                                                 band_edge_kind edge_kind) {
    if (samples.size() < 6) {
        throw std::invalid_argument("effective-mass fit needs at least six samples");
    }
    if (!(lattice_constant_m > 0.0) || !std::isfinite(lattice_constant_m)) {
        throw std::invalid_argument("lattice constant must be positive");
    }

    const double          reduced_to_k_m = 2.0 * uepm::constants::pi / lattice_constant_m;
    const Eigen::Vector3d k0             = to_eigen(k0_reduced);
    Eigen::MatrixXd       design         = build_quadratic_design(samples, k0_reduced, lattice_constant_m);
    Eigen::VectorXd       kinetic        = kinetic_energies(samples, edge_energy_eV, edge_kind);

    double          alpha = 0.0;
    Eigen::Matrix3d hessian;
    for (int iteration = 0; iteration < 100; ++iteration) {
        Eigen::VectorXd gamma_target(samples.size());
        for (std::size_t i = 0; i < samples.size(); ++i) {
            const double energy                        = kinetic(static_cast<Eigen::Index>(i));
            gamma_target(static_cast<Eigen::Index>(i)) = energy * (1.0 + alpha * energy);
        }

        hessian = fit_hessian_from_gamma_target(design, gamma_target);

        double alpha_numerator   = 0.0;
        double alpha_denominator = 0.0;
        for (std::size_t i = 0; i < samples.size(); ++i) {
            const Eigen::Vector3d dk_reduced = to_eigen(samples[i].k_reduced) - k0;
            const Eigen::Vector3d q          = dk_reduced * reduced_to_k_m;
            const double          energy     = kinetic(static_cast<Eigen::Index>(i));
            if (energy <= kinetic_tolerance_eV) {
                continue;
            }
            const double gamma = 0.5 * q.transpose() * hessian * q;
            alpha_numerator += energy * energy * (gamma - energy);
            alpha_denominator += energy * energy * energy * energy;
        }
        const double next_alpha = alpha_denominator > 0.0 ? alpha_numerator / alpha_denominator : 0.0;
        if (std::abs(next_alpha - alpha) < 1.0e-10) {
            alpha = next_alpha;
            break;
        }
        alpha = next_alpha;
    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(hessian);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("effective-mass fit failed to diagonalize the Hessian");
    }

    const auto masses = masses_from_hessian(solver);

    double residual2 = 0.0;
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const Eigen::Vector3d dk_reduced = to_eigen(samples[i].k_reduced) - k0;
        const Eigen::Vector3d q          = dk_reduced * reduced_to_k_m;
        const double          energy     = kinetic(static_cast<Eigen::Index>(i));
        const double          gamma      = 0.5 * q.transpose() * hessian * q;
        const double          residual   = gamma - energy * (1.0 + alpha * energy);
        residual2 += residual * residual;
    }

    effective_mass_fit_result result;
    result.k0_reduced              = k0_reduced;
    result.edge_energy_eV          = edge_energy_eV;
    result.edge_kind               = edge_kind;
    result.non_parabolicity_eV_inv = alpha;
    result.rms_error_meV           = 1000.0 * std::sqrt(residual2 / static_cast<double>(samples.size()));
    result.mass_rms_error_meV      = result.rms_error_meV;
    result.alpha_rms_error_meV     = result.rms_error_meV;
    result.sample_count            = samples.size();
    result.mass_sample_count       = samples.size();
    result.alpha_sample_count      = samples.size();
    result.hessian_eV_m2           = to_array(hessian);
    result.principal_axes          = to_array(solver.eigenvectors());
    result.principal_masses_m0     = masses;
    return result;
}

effective_mass_fit_result fit_effective_mass_tensor(const std::vector<valley_fit_sample>& samples,
                                                    const Vector3D<double>&               k0_reduced,
                                                    double                                edge_energy_eV,
                                                    double                                lattice_constant_m,
                                                    band_edge_kind                        edge_kind) {
    if (samples.size() < 6) {
        throw std::invalid_argument("effective-mass tensor fit needs at least six samples");
    }
    if (!(lattice_constant_m > 0.0) || !std::isfinite(lattice_constant_m)) {
        throw std::invalid_argument("lattice constant must be positive");
    }

    const Eigen::MatrixXd design  = build_quadratic_design(samples, k0_reduced, lattice_constant_m);
    const Eigen::VectorXd kinetic = kinetic_energies(samples, edge_energy_eV, edge_kind);
    const Eigen::VectorXd coeff   = design.colPivHouseholderQr().solve(kinetic);
    Eigen::Matrix3d       hessian;
    hessian << coeff(0), coeff(3), coeff(4), coeff(3), coeff(1), coeff(5), coeff(4), coeff(5), coeff(2);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(hessian);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("effective-mass tensor fit failed to diagonalize the Hessian");
    }

    double residual2 = 0.0;
    for (Eigen::Index i = 0; i < design.rows(); ++i) {
        const double residual = (design.row(i) * coeff)(0) - kinetic(i);
        residual2 += residual * residual;
    }

    effective_mass_fit_result result;
    result.k0_reduced          = k0_reduced;
    result.edge_energy_eV      = edge_energy_eV;
    result.edge_kind           = edge_kind;
    result.rms_error_meV       = 1000.0 * std::sqrt(residual2 / static_cast<double>(samples.size()));
    result.mass_rms_error_meV  = result.rms_error_meV;
    result.sample_count        = samples.size();
    result.mass_sample_count   = samples.size();
    result.hessian_eV_m2       = to_array(hessian);
    result.principal_axes      = to_array(solver.eigenvectors());
    result.principal_masses_m0 = masses_from_hessian(solver);
    return result;
}

effective_mass_fit_result fit_nonparabolicity_with_fixed_mass(const std::vector<valley_fit_sample>& samples,
                                                              const effective_mass_fit_result&      mass_fit,
                                                              double                                lattice_constant_m,
                                                              double max_kinetic_energy_eV,
                                                              bool   clamp_nonnegative) {
    if (!(max_kinetic_energy_eV > 0.0) || !std::isfinite(max_kinetic_energy_eV)) {
        throw std::invalid_argument("maximum kinetic energy must be positive");
    }
    if (!(lattice_constant_m > 0.0) || !std::isfinite(lattice_constant_m)) {
        throw std::invalid_argument("lattice constant must be positive");
    }

    const double          sign           = edge_sign(mass_fit.edge_kind);
    const double          reduced_to_k_m = 2.0 * uepm::constants::pi / lattice_constant_m;
    const Eigen::Vector3d k0             = to_eigen(mass_fit.k0_reduced);
    const Eigen::Matrix3d hessian        = to_matrix(mass_fit.hessian_eV_m2);

    double      alpha_numerator   = 0.0;
    double      alpha_denominator = 0.0;
    std::size_t used_samples      = 0;
    for (const auto& sample : samples) {
        const double energy = sign * (sample.energy_eV - mass_fit.edge_energy_eV);
        if (energy <= kinetic_tolerance_eV || energy > max_kinetic_energy_eV) {
            continue;
        }
        const Eigen::Vector3d q     = (to_eigen(sample.k_reduced) - k0) * reduced_to_k_m;
        const double          gamma = 0.5 * q.transpose() * hessian * q;
        alpha_numerator += energy * energy * (gamma - energy);
        alpha_denominator += energy * energy * energy * energy;
        ++used_samples;
    }
    if (used_samples == 0 || alpha_denominator <= 0.0) {
        throw std::runtime_error("no alpha-fit samples inside the requested energy window");
    }

    double alpha = alpha_numerator / alpha_denominator;
    if (clamp_nonnegative && alpha < 0.0) {
        alpha = 0.0;
    }

    double residual2 = 0.0;
    for (const auto& sample : samples) {
        const double energy = sign * (sample.energy_eV - mass_fit.edge_energy_eV);
        if (energy <= kinetic_tolerance_eV || energy > max_kinetic_energy_eV) {
            continue;
        }
        const Eigen::Vector3d q        = (to_eigen(sample.k_reduced) - k0) * reduced_to_k_m;
        const double          gamma    = 0.5 * q.transpose() * hessian * q;
        const double          residual = gamma - energy * (1.0 + alpha * energy);
        residual2 += residual * residual;
    }

    auto result                    = mass_fit;
    result.non_parabolicity_eV_inv = alpha;
    result.rms_error_meV           = 1000.0 * std::sqrt(residual2 / static_cast<double>(used_samples));
    result.alpha_rms_error_meV     = result.rms_error_meV;
    result.sample_count            = used_samples;
    result.alpha_sample_count      = used_samples;
    return result;
}

effective_mass_fit_result fit_effective_mass_then_nonparabolicity(const std::vector<valley_fit_sample>& mass_samples,
                                                                  const std::vector<valley_fit_sample>& alpha_samples,
                                                                  const Vector3D<double>&               k0_reduced,
                                                                  double                                edge_energy_eV,
                                                                  double         lattice_constant_m,
                                                                  band_edge_kind edge_kind,
                                                                  double         max_kinetic_energy_eV,
                                                                  bool           clamp_nonnegative) {
    const auto mass_fit =
        fit_effective_mass_tensor(mass_samples, k0_reduced, edge_energy_eV, lattice_constant_m, edge_kind);
    return fit_nonparabolicity_with_fixed_mass(alpha_samples,
                                               mass_fit,
                                               lattice_constant_m,
                                               max_kinetic_energy_eV,
                                               clamp_nonnegative);
}

std::string to_string(band_edge_kind edge_kind) { return edge_kind == band_edge_kind::minimum ? "minimum" : "maximum"; }

band_edge_kind band_edge_kind_from_string(const std::string& value) {
    if (value == "min" || value == "minimum" || value == "electron") {
        return band_edge_kind::minimum;
    }
    if (value == "max" || value == "maximum" || value == "hole") {
        return band_edge_kind::maximum;
    }
    throw std::invalid_argument("edge kind must be one of: min, minimum, electron, max, maximum, hole");
}

}  // namespace uepm::pseudopotential
