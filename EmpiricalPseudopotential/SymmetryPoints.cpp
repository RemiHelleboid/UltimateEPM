#include "SymmetryPoints.h"

namespace EmpiricalPseudopotential
{


	SymmetryPoints::SymmetryPoints()
		:
		path{ {"K", "W", "X", "G", "L", "W"} }
	{
		symmetryPoints["L"] = SymmetryPoint("L", Vector3D<double>(0.5, 0.5, 0.5));
		symmetryPoints["G"] = SymmetryPoint("G", Vector3D<double>(0., 0., 0.));
		symmetryPoints["X"] = SymmetryPoint("X", Vector3D<double>(1., 0., 0.));
		symmetryPoints["W"] = SymmetryPoint("W", Vector3D<double>(1., 0.5, 0.));
		symmetryPoints["K"] = SymmetryPoint("K", Vector3D<double>(0.75, 0.75, 0.));
		symmetryPoints["U"] = SymmetryPoint("U", Vector3D<double>(1., 0.25, 0.25));
	}


	SymmetryPoints::~SymmetryPoints()
	{
	}


	std::vector<Vector3D<double>> SymmetryPoints::GeneratePoints(int nrPoints, std::vector<unsigned int>& symmetryPointsPositions)
	{
		std::vector<Vector3D<double>> result;

		symmetryPointsPositions.clear();
		symmetryPointsPositions.reserve(path.size());

		if (nrPoints <= path.size() * 2 + 1) return result;

		result.reserve(nrPoints);

		double length = 0;

		for (int i = 1; i < path.size(); ++i)
		{
			const Vector3D<double> dif = symmetryPoints[path[i]].position - symmetryPoints[path[i - 1]].position;
			length += dif.Length();
		}

		const double stepSize = length / (nrPoints - 1);

		for (int i = 1; i < path.size(); ++i)
		{
			const Vector3D<double> startPos = symmetryPoints[path[i - 1]].position;
			const Vector3D<double> dif = symmetryPoints[path[i]].position - startPos;
			const double difLength = dif.Length();

			Vector3D<double> stepVec = dif / difLength * stepSize;

			if (1 == i) symmetryPointsPositions.push_back(0);
			else symmetryPointsPositions.push_back(static_cast<unsigned int>(result.size() + 1));

			for (Vector3D<double> pos = startPos; (pos - startPos).Length() < difLength; pos += stepVec)
				result.push_back(pos);
		}


		return std::move(result);
	}

}