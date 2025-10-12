#include "SymmetryPoints.h"

namespace uepm::pseudopotential {

SymmetryPoints::SymmetryPoints() {
    symmetryPoints["L"] = SymmetryPoint("L", Vector3D<double>(0.5, 0.5, 0.5));
    symmetryPoints["G"] = SymmetryPoint("G", Vector3D<double>(0., 0., 0.));
    symmetryPoints["X"] = SymmetryPoint("X", Vector3D<double>(1., 0., 0.));
    symmetryPoints["W"] = SymmetryPoint("W", Vector3D<double>(1., 0.5, 0.));
    symmetryPoints["K"] = SymmetryPoint("K", Vector3D<double>(0.75, 0.75, 0.));
    symmetryPoints["U"] = SymmetryPoint("U", Vector3D<double>(1., 0.25, 0.25));
}

/**
 * @brief Generates the list of k points for the given path of "symetry points".
 *
 * The path will be for example {"G", "X", "W", "L", "G", "K", "X"}
 * The number of points will be the total of points in the output.
 * The symmetryPointsPositions will be the positions (indexes) of the symmetry points in the output.
 *
 * The points are spaced in such way that the number of points between two symmetry points is "proportional" to the distance between them.
 *
 * @param path
 * @param nrPoints
 * @param symmetryPointsPositions
 * @return std::vector<Vector3D<double>>
 */
std::vector<Vector3D<double>> SymmetryPoints::GeneratePoints(const std::vector<std::string>& path,
                                                             unsigned int                    nrPoints,
                                                             std::vector<unsigned int>&      symmetryPointsPositions) {
    std::vector<Vector3D<double>> result;

    symmetryPointsPositions.clear();
    symmetryPointsPositions.reserve(path.size());

    if (nrPoints <= path.size() * 2 + 1) {
        std::cout << "Error: nrPoints must be greater than twice the number of symmetry points in the path + 1" << std::endl;
        return result;
    };

    result.reserve(nrPoints);

    // Calculate the length of each path segment, and add it to the total length.
    double length = 0.0;
    for (std::size_t index_HS_point = 1; index_HS_point < path.size(); ++index_HS_point) {
        const Vector3D<double> dif = symmetryPoints[path[index_HS_point]].position - symmetryPoints[path[index_HS_point - 1]].position;
        length += dif.Length();
    }

    const double stepSize = length / (nrPoints - 1.);

    for (std::size_t indexHSpoint = 1; indexHSpoint < path.size(); ++indexHSpoint) {
        const Vector3D<double> startPos  = symmetryPoints[path[indexHSpoint - 1ULL]].position;
        const Vector3D<double> dif       = symmetryPoints[path[indexHSpoint]].position - startPos;
        const double           difLength = dif.Length();

        Vector3D<double> stepVec = dif / difLength * stepSize;

        if (indexHSpoint == 1) {
            symmetryPointsPositions.push_back(0.0);
        } else {
            symmetryPointsPositions.push_back(static_cast<unsigned int>(result.size() + 1));
        }

        // Fill the path segment with k points spaced by stepVec.
        for (Vector3D<double> pos = startPos; (pos - startPos).Length() < difLength; pos += stepVec) {
            result.push_back(pos);
        }
    }

    return result;
}

}  // namespace uepm::pseudopotential