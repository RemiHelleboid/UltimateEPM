#include "BandStructure.h"

#include <cfloat>
#include <fstream>
#include <iostream>

#include "Hamiltonian.h"

namespace EmpiricalPseudopotential {

BandStructure::BandStructure() { basisVectors.reserve(127); }

bool BandStructure::GenerateBasisVectors(unsigned int nearestNeighborsNumber) {
    static const std::vector<unsigned int> G2{
        0,   3,   4,   8,   11,  12,  16,  19,  20,  24,  27,  32,  35,  36,  40,  43,  44,  48,  51,  52,  56,  59,  67,  68,
        75,  76,  80,  83,  84,  88,  96,  99,  104, 107, 108, 115, 116, 120, 123, 128, 131, 132, 136, 139, 140, 144, 147, 152,
        155, 160, 163, 164, 168, 171, 172, 176, 179, 180, 184, 187, 192, 195, 196, 200, 203, 204, 208, 211, 212, 216, 219, 224,
        227, 228, 232, 236, 243, 244, 248, 251, 259, 260, 264, 267, 268, 272, 275, 276, 280, 283, 291, 296, 299, 300, 304, 307,
        308, 312, 315, 320, 323, 324, 331, 332, 339, 355, 356, 360, 363, 371, 376, 384, 387, 395, 420, 451};

    if (nearestNeighborsNumber < 2 || nearestNeighborsNumber > G2.size()) return false;
    const unsigned int nearestNeighbors = nearestNeighborsNumber - 1;
    basisVectors.clear();
    const int           size = static_cast<int>(ceil(sqrt(static_cast<double>(G2[nearestNeighbors]))));
    const Vector3D<int> b1(-1, 1, 1), b2(1, -1, 1), b3(1, 1, -1);

    for (int i = -size; i <= size; ++i)
        for (int j = -size; j <= size; ++j)
            for (int k = -size; k <= size; ++k) {
                const Vector3D<int> vect        = b1 * i + b2 * j + b3 * k;  // reciprocal lattice vector
                const double        vectSquared = vect * vect;

                if (vectSquared <= G2[nearestNeighbors])  // if it's under the cutoff length, add it
                    basisVectors.push_back(vect);
            }

    return true;
}

void BandStructure::Initialize(const Material&           material,
                               std::size_t               nb_bands,
                               std::vector<std::string>& path,
                               unsigned int              nbPoints,
                               unsigned int              nearestNeighborsNumber) {
    m_material               = material;
    m_nb_bands               = nb_bands;
    m_path                   = path;
    m_nb_points               = nbPoints;
    m_nearestNeighborsNumber = nearestNeighborsNumber;
    m_kpoints.clear();
    m_results.clear();
    m_kpoints.reserve(m_nb_points);
    m_results.reserve(m_nb_points);

    if (!GenerateBasisVectors(nearestNeighborsNumber)) {
        throw std::runtime_error("BandStructure::Initialize: GenerateBasisVectors failed");
    }
    m_kpoints = symmetryPoints.GeneratePoints(m_path, m_nb_points, symmetryPointsPositions);
}

std::vector<std::vector<double>> BandStructure::Compute() {
    std::cout << "Computing band structure..." << std::endl;
    std::vector<std::vector<double>> res;
    m_results.clear();

    Hamiltonian hamiltonian(m_material, basisVectors);

    for (unsigned int i = 0; i < m_nb_points; ++i) {
        hamiltonian.SetMatrix(m_kpoints[i]);
        hamiltonian.Diagonalize();

        const Eigen::VectorXd& eigenvals = hamiltonian.eigenvalues();

        m_results.emplace_back();
        m_results.back().reserve(m_nb_bands);
        for (unsigned int level = 0; level < m_nb_bands && level < eigenvals.rows(); ++level) {
            m_results.back().push_back(eigenvals(level));
        }
    }
    std::cout << "Done!" << std::endl;
    return std::move(res);
}

double BandStructure::AdjustValues() {
    constexpr double hartree_to_eV = 27.2113845;
    double           maxValValence;
    double           minValConduction;

    double bandgap = 0;

    if (FindBandgap(m_results, maxValValence, minValConduction)) bandgap = minValConduction - maxValValence;

    // adjust values to a guessed zero
    for (auto& p : m_results)
        for (auto& v : p) {
            v -= maxValValence;

            // computation is done with atomic units
            // results are in Hartree, here they are converted to eV
            v *= 27.211385;
        }

    return bandgap * 27.211385;
}

bool BandStructure::FindBandgap(const std::vector<std::vector<double>>& results, double& maxValValence, double& minValConduction) {
    maxValValence = DBL_MIN;
    if (results.empty() || results.front().size() < 2) return false;

    const unsigned int nrLevels       = static_cast<unsigned int>(results.front().size());
    double             fallbackMaxVal = 0;

    for (unsigned int levelLow = 2; levelLow < nrLevels - 1; ++levelLow) {
        maxValValence    = DBL_MIN;
        minValConduction = DBL_MAX;

        for (auto& p : results) {
            const double valLow  = p[levelLow];
            const double valHigh = p[levelLow + 1ULL];

            maxValValence    = std::max(maxValValence, valLow);
            minValConduction = std::min(minValConduction, valHigh);
        }

        if (3 == levelLow) fallbackMaxVal = maxValValence;

        if (maxValValence + 0.35 < minValConduction) return true;
    }

    maxValValence = fallbackMaxVal;

    return false;
}

void BandStructure::print_results() const {
    for (auto& p : m_results) {
        for (auto& v : p)
            std::cout << v << " ";
        std::cout << std::endl;
    }
}

void BandStructure::export_result_in_file(const std::string& filename) const {
    std::ofstream file(filename);

    for (auto& p : m_results) {
        for (auto& v : p)
            file << v << " ";
        file << std::endl;
    }
}

std::string BandStructure::path_band_filename() const {
    std::string path_string;
    if (m_path.empty()) {
        path_string = "";
    } else {
        for (auto& point : m_path) {
            path_string += point;
        }
    }
    std::string filename = "EEP_" + m_material.name + "_nb_bands_" + std::to_string(m_results.front().size()) + "_path_" + path_string +
                           "_size_basis_" + std::to_string(basisVectors.size());
    return filename;
}

void export_vector_bands_result_in_file(const std::string& filename, std::vector<std::vector<double>> results) {
    std::ofstream file(filename);

    for (auto& p : results) {
        for (auto& v : p)
            file << v << ",";
        file << std::endl;
    }
}

}  // namespace EmpiricalPseudopotential