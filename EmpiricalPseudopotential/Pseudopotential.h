#pragma once

#include <complex>
#include "Vector3D.h"

namespace EmpiricalPseudopotential
{


	class Pseudopotential
	{
	public:
		Pseudopotential(double V3S = 0, double V8S = 0, double V11S = 0, double V3A = 0, double V4A = 0, double V11A = 0);

	protected:
		double m_V3S;
		double m_V8S;
		double m_V11S;

		double m_V3A;
		double m_V4A;
		double m_V11A;

	public:
		std::complex<double> GetValue(const Vector3D<int>& G, const Vector3D<double>& tau = Vector3D<double>(1. / 8., 1. / 8., 1. / 8.)) const;
	};

}