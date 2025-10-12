#pragma once

#include <complex>

#include "Vector3D.h"

namespace uepm::pseudopotential {

class Pseudopotential {
 public:
    Pseudopotential(double V3S  = 0,
                    double V4S  = 0,
                    double V8S  = 0,
                    double V11S = 0,
                    double V3A  = 0,
                    double V4A  = 0,
                    double V8A  = 0,
                    double V11A = 0);

 protected:
    double m_V3S;
    double m_V4S;
    double m_V8S;
    double m_V11S;

    double m_V3A;
    double m_V4A;
    double m_V8A;
    double m_V11A;

 public:
    std::complex<double> GetValue(const Vector3D<int>& G, const Vector3D<double>& tau, double lattice_constant) const;
    void                 print_parameters() const;
};

}  // namespace uepm::pseudopotential