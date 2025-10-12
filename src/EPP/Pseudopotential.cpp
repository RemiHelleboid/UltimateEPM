#include "Pseudopotential.h"

#define _USE_MATH_DEFINES 1
#include <array>
#include <cmath>
#include <complex>
#include <numbers>

#include "physical_constants.hpp"

namespace uepm::pseudopotential {

Pseudopotential::Pseudopotential(double V3S, double V4S, double V8S, double V11S, double V3A, double V4A, double V8A, double V11A)
    : m_V3S(V3S),
      m_V4S(V4S),
      m_V8S(V8S),
      m_V11S(V11S),
      m_V3A(V3A),
      m_V4A(V4A),
      m_V8A(V8A),
      m_V11A(V11A) {}

#include <array>
#include <cmath>
#include <complex>
#include <numbers>

std::complex<double> Pseudopotential::GetValue(const Vector3D<int>& G, const Vector3D<double>& tau, double lattice_constant) const {
    // Optional: enforce fcc selection rule (all-even or all-odd Miller indices)
    const bool same_parity = ((G.X & 1) == (G.Y & 1)) && ((G.Y & 1) == (G.Z & 1));
    if (!same_parity) return {0.0, 0.0};

    const int    G2   = G * G;  // integer dot product h^2 + k^2 + l^2
    const double kfac = 2.0 * std::numbers::pi_v<double> / lattice_constant;
    const double Gtau = kfac * (tau * G);  // dimensionless phase

    double VS = 0.0, VA = 0.0;
    switch (G2) {
        case 3:
            VS = m_V3S;
            VA = m_V3A;
            break;
        case 4:
            VS = m_V4S;
            VA = m_V4A;
            break;
        case 8:
            VS = m_V8S;
            VA = m_V8A;
            break;
        case 11:
            VS = m_V11S;
            VA = m_V11A;
            break;
        default: 
            break;
    }

    // V(G) = VS cos(G·τ) + i VA sin(G·τ)
    const double c = std::cos(Gtau);
    const double s = std::sin(Gtau);
    return {VS * c, VA * s};
}

void Pseudopotential::print_parameters() const {
    std::cout << "V3S = " << m_V3S << std::endl;
    std::cout << "V4S = " << m_V4S << std::endl;
    std::cout << "V8S = " << m_V8S << std::endl;
    std::cout << "V11S = " << m_V11S << std::endl;
    std::cout << "V3A = " << m_V3A << std::endl;
    std::cout << "V4A = " << m_V4A << std::endl;
    std::cout << "V8A = " << m_V8A << std::endl;
    std::cout << "V11A = " << m_V11A << std::endl;
}

}  // namespace uepm::pseudopotential