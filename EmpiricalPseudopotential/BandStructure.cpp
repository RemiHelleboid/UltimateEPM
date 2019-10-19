#include "BandStructure.h"
#include "Hamiltonian.h"

namespace EmpiricalPseudopotential
{

	BandStructure::BandStructure()
		: nearestNeighbors(9),
		G2{ 0, 3, 4, 8, 11, 12, 16, 19, 20, 24 }
	{
		basisVectors.reserve(127);
	}


	bool BandStructure::GenerateBasisVectors(unsigned int nearestNeighborsNumber)
	{
		if (nearestNeighborsNumber < 2 || nearestNeighborsNumber > 10) return false;

		nearestNeighbors = nearestNeighborsNumber - 1;
		basisVectors.clear();


		const int size = static_cast<int>(ceil(sqrt(static_cast<double>(G2[nearestNeighbors]))));
		
		/*
		for (int i = -size; i <= size; ++i)
			for (int j = -size; j <= size; ++j)
				for (int k = -size; k <= size; ++k)
				{
					const Vector3D<int> vect(i, j, k);
					const double vectSquared = vect * vect;

					for (unsigned int nearestNeighbor = 0; nearestNeighbor <= nearestNeighbors; ++nearestNeighbor)
						if (vectSquared == G2[nearestNeighbor])
						{
							basisVectors.push_back(vect);
							break;
						}
				}
		*/


		// the other way is commented above
		// the following should be easier to understand

		// the Bravais lattice is a fcc lattice (with 2 atoms - the anion and cation - in the unit cell)
		// the reciprocal lattice is a bcc lattice 

		// the basis vectors for the reciprocal space
		
		const Vector3D<int> b1(-1, 1, 1), b2(1, -1, 1), b3(1, 1, -1);

		// you can also get them from the Bravais lattice vectors
		// like this:

		// Bravais lattice vectors for fcc:
		//const Vector3D<double> a1(0, 0.5, 0.5), a2(0.5, 0, 0.5), a3(0.5, 0.5, 0);

		// reciprocal lattice:
		//const Vector3D<double> b1 = a2 % a3 / (a1 * (a2 % a3)); 		
		//const Vector3D<double> b2 = a3 % a1 / (a2 * (a3 % a1));
		//const Vector3D<double> b3 = a1 % a2 / (a3 * (a1 % a2));
		

		for (int i = -size; i <= size; ++i)
			for (int j = -size; j <= size; ++j)
				for (int k = -size; k <= size; ++k)
				{
					const Vector3D<int> vect = b1 * i + b2 * j + b3 * k; // reciprocal lattice vector
					const double vectSquared =  vect * vect;

					if (vectSquared <= G2[nearestNeighbors]) // if it's under the cutoff length, add it
						basisVectors.push_back(vect);
				}
				
		return true;
	}


	void BandStructure::Initialize(std::vector<std::string> path, unsigned int nrPoints,  unsigned int nearestNeighborsNumber)
	{
		kpoints.clear();
		results.clear();

		kpoints.reserve(nrPoints);
		results.reserve(nrPoints);

		m_path.swap(path);

		GenerateBasisVectors(nearestNeighborsNumber);

		kpoints = symmetryPoints.GeneratePoints(m_path, nrPoints, symmetryPointsPositions);
	}


	std::vector<std::vector<double>> BandStructure::Compute(const Material& material, unsigned int startPoint, unsigned int endPoint, unsigned int nrLevels, std::atomic_bool& terminate)
	{	
		std::vector<std::vector<double>> res;

		Hamiltonian hamiltonian(material, basisVectors);

		for (unsigned int i = startPoint; i < endPoint && !terminate; ++i)
		{
			hamiltonian.SetMatrix(kpoints[i]);
			hamiltonian.Diagonalize();

			const Eigen::VectorXd& eigenvals = hamiltonian.eigenvalues();

			res.emplace_back();
			res.back().reserve(nrLevels);

			for (unsigned int level = 0; level < nrLevels && level < eigenvals.rows(); ++level)
				res.back().push_back(eigenvals(level));
		}

		return std::move(res);
	}


	double BandStructure::AdjustValues()
	{
		double maxValValence;
		double minValConduction;

		double bandgap = 0;

		if (FindBandgap(results, maxValValence, minValConduction))
			bandgap = minValConduction - maxValValence;

		// adjust values to a guessed zero
		for (auto& p : results)
			for (auto& v : p)
			{
				v -= maxValValence;
				
				// computation is done with atomic units
				// results are in Hartree, here they are converted to eV
				v *= 27.211385; 
			}

		return bandgap * 27.211385;
	}


	bool BandStructure::FindBandgap(const std::vector<std::vector<double>>& results, double& maxValValence, double& minValConduction)
	{
		maxValValence = DBL_MIN;
		if (results.empty() || results.front().size() < 2) return false;


		const unsigned int nrLevels = static_cast<unsigned int>(results.front().size());
		double fallbackMaxVal = 0;

		for (unsigned int levelLow = 2; levelLow < nrLevels - 1; ++levelLow)
		{
			maxValValence = DBL_MIN;
			minValConduction = DBL_MAX;

			for (auto& p : results)
			{
				const double valLow = p[levelLow];
				const double valHigh = p[levelLow + 1ULL];

				maxValValence = std::max(maxValValence, valLow);
				minValConduction = std::min(minValConduction, valHigh);
			}

			if (3 == levelLow) fallbackMaxVal = maxValValence;

			if (maxValValence + 0.35 < minValConduction)
				return true;
		}

		maxValValence = fallbackMaxVal;

		return false;
	}

}