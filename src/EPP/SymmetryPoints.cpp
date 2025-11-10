/**
 * @file SymmetryPoints.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2025-11-01
 * 
 * 
 */

#include "SymmetryPoints.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

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
    using Vec = Vector3D<double>;
    std::vector<Vec> result;
    symmetryPointsPositions.clear();

    // --- Validation ---
    const std::size_t M = path.size();
    if (M < 2) {
        throw std::invalid_argument("Path must contain at least two symmetry points.");
    }
    if (nrPoints < M) {
        throw std::invalid_argument("nrPoints must be at least the number of nodes in the path.");
    }

    std::vector<Vec> nodes;
    nodes.reserve(M);
    for (const auto& key : path) {
        auto it = symmetryPoints.find(key);
        if (it == symmetryPoints.end()) {
            throw std::invalid_argument("Unknown symmetry point label: " + key);
        }
        nodes.push_back(it->second.position);
    }

    std::vector<double> segLen(M - 1, 0.0);
    for (std::size_t i = 1; i < M; ++i) {
        segLen[i - 1] = (nodes[i] - nodes[i - 1]).Length();
    }
    const double totalLen = std::accumulate(segLen.begin(), segLen.end(), 0.0);

    const unsigned R = nrPoints - static_cast<unsigned>(M);

    if (totalLen == 0.0) {
        result.reserve(nrPoints);
        for (std::size_t j = 0; j < M; ++j) {
            symmetryPointsPositions.push_back(static_cast<unsigned>(result.size()));
            result.push_back(nodes[j]);
        }
        // Add remaining duplicates of the last node
        for (unsigned i = 0; i < R; ++i) {
            result.push_back(nodes.back());
        }
        return result;
    }

    std::vector<unsigned>    m(M - 1, 0);
    std::vector<double>      want(M - 1, 0.0);
    std::vector<std::size_t> order(M - 1);
    for (std::size_t s = 0; s < M - 1; ++s) {
        want[s]  = (segLen[s] / totalLen) * static_cast<double>(R);
        m[s]     = static_cast<unsigned>(std::floor(want[s]));
        order[s] = s;
    }
    unsigned assigned = std::accumulate(m.begin(), m.end(), 0u);
    unsigned rem      = R - assigned;

    // Give remaining points to the segments with largest fractional parts
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        double fa = want[a] - std::floor(want[a]);
        double fb = want[b] - std::floor(want[b]);
        return fa > fb;  // descending fractional part
    });
    for (unsigned i = 0; i < rem; ++i) {
        ++m[order[i]];
    }

    // Build result ensuring every node is exact and appears once.
    result.reserve(nrPoints);
    symmetryPointsPositions.reserve(M);

    // First node
    symmetryPointsPositions.push_back(static_cast<unsigned>(result.size()));
    result.push_back(nodes.front());

    for (std::size_t s = 0; s < M - 1; ++s) {
        const Vec&     A  = nodes[s];
        const Vec&     B  = nodes[s + 1];
        const unsigned ms = m[s];

        if (ms > 0) {
            const Vec    d     = B - A;
            const double denom = static_cast<double>(ms + 1);
            for (unsigned i = 1; i <= ms; ++i) {
                double t = static_cast<double>(i) / denom;
                result.push_back(A + d * t);
                // std::cout << "Added point at t=" << t <<  " : " << result.back().X << " " << result.back().Y << " " << result.back().Z << std::endl;
            }
        }

        symmetryPointsPositions.push_back(static_cast<unsigned>(result.size()));
        result.push_back(B);
    }
    // for (auto &V : result) {
    //     std::cout << "K-point: " << V.X << " " << V.Y << " " << V.Z << std::endl;
    // }


    return result;
}

}  // namespace uepm::pseudopotential