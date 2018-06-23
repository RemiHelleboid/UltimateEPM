#pragma once

#include <vector>
#include <atomic>

#include "Vector3D.h"

#include "Material.h"
#include "SymmetryPoints.h"

namespace EmpiricalPseudopotential
{

	class BandStructure
	{
	public:
		BandStructure();

		Materials materials;
		SymmetryPoints symmetryPoints;

		std::vector<std::vector<double>> results;
		std::vector<unsigned int> symmetryPointsPositions;


		void Initialize(unsigned int nrPoints = 600,  unsigned int nearestNeighborsNumber = 10);
		std::vector<std::vector<double>> Compute(const Material& material, unsigned int startPoint, unsigned int endPoint, unsigned int nrLevels, std::atomic_bool& terminate);

		double AdjustValues();

		unsigned int GetPointsNumber() const { return static_cast<unsigned int>(kpoints.size()); }
	private:
		std::vector<Vector3D<int>> basisVectors;

		std::vector<Vector3D<double>> kpoints;

		unsigned int nearestNeighbors;
		std::vector<unsigned int> G2;

		static bool FindBandgap(const std::vector<std::vector<double>>& results, double& maxValValence, double& minValConduction);
		bool GenerateBasisVectors(unsigned int nearestNeighborsNumber);
	};

}

