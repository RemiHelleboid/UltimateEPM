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

    void Initialize(const Material&           material,
                    std::size_t               nb_bands,
                    std::vector<std::string>& path,
                    unsigned int              nrPoints,
                    unsigned int              nearestNeighborsNumber);

    //  void Initialize(std::vector<Vector3D<double>> list_k_points, unsigned int nearestNeighborsNumber = 10);
    std::vector<std::vector<double>> Compute();
    std::vector<std::vector<double>> Compute_parralel(int nb_threads);

    double AdjustValues();

    unsigned int GetPointsNumber() const { return static_cast<unsigned int>(m_kpoints.size()); }

    const std::vector<std::string>& GetPath() const { return m_path; }

    void print_results() const;
    void export_result_in_file(const std::string& filename) const;

    std::string path_band_filename() const;

 private:
    std::vector<std::string> m_path;
    unsigned int             m_nb_points;

    Material     m_material;
    unsigned int m_nb_bands;
    unsigned int m_nearestNeighborsNumber;

    std::vector<Vector3D<int>>    basisVectors;
    std::vector<Vector3D<double>> m_kpoints;

    static bool FindBandgap(const std::vector<std::vector<double>>& results, double& maxValValence, double& minValConduction);
    bool        GenerateBasisVectors(unsigned int nearestNeighborsNumber);
};

void export_vector_bands_result_in_file(const std::string& filename, std::vector<std::vector<double>>);

}  // namespace EmpiricalPseudopotential
