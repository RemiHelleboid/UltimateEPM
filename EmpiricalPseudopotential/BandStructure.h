#pragma once

#include <atomic>
#include <vector>

#include "Material.h"
#include "SymmetryPoints.h"
#include "Vector3D.h"

namespace EmpiricalPseudopotential {

class BandStructure {
 public:
    BandStructure();

    Materials      materials;
    SymmetryPoints symmetryPoints;

    std::vector<std::vector<double>> m_results;
    std::vector<unsigned int>        symmetryPointsPositions;

    void Initialize(std::vector<std::string> path, unsigned int nrPoints = 600, unsigned int nearestNeighborsNumber = 10);
    std::vector<std::vector<double>> Compute(const Material&   material,
                                             unsigned int      startPoint,
                                             unsigned int      endPoint,
                                             unsigned int      nrLevels,
                                             std::atomic_bool& terminate);

    std::vector<std::vector<double>> Compute(const Material&   material,
                                             unsigned int      startPoint,
                                             unsigned int      endPoint,
                                             unsigned int      nrLevels);

    double AdjustValues();

    unsigned int GetPointsNumber() const { return static_cast<unsigned int>(kpoints.size()); }


    const std::vector<std::string>& GetPath() const { return m_path; }

   void print_results() const;
   void export_result_in_file(const std::string& filename) const;

 private:
    std::vector<std::string>      m_path;
    std::vector<Vector3D<int>>    basisVectors;
    std::vector<Vector3D<double>> kpoints;

    static bool FindBandgap(const std::vector<std::vector<double>>& results, double& maxValValence, double& minValConduction);
    bool        GenerateBasisVectors(unsigned int nearestNeighborsNumber);
};

void export_vector_bands_result_in_file(const std::string& filename, std::vector<std::vector<double>>) ;

}  // namespace EmpiricalPseudopotential
