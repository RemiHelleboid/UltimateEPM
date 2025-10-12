// main.cpp
// Build: see CMakeLists.txt below
// Translates the provided Python into modern C++ using the Gmsh C++ API.
// - Creates initial 3D points (Gamma, L, X, K, W, U)
// - Tetrahedralizes with gmsh::algorithm::tetrahedralize
// - Optionally refines the mesh N times
// - Extracts node coordinates
// - Applies the 48 cubic symmetry operations (axis permutations x reflections)
// - Creates a new discrete volume, tetrahedralizes the symmetrized point cloud
// - Writes the mesh to MSH and optionally opens the GUI
//
// Usage: ./delaunay_sym [nb_levels=0] [outname=delaunay_mesh_refined.msh] [--nogui]

#include <gmsh.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {
using Mat3 = std::array<std::array<int, 3>, 3>;
using Vec3 = std::array<double, 3>;

struct VecHash {
    std::size_t operator()(Vec3 const &v) const noexcept {
        // Simple FNV-1a over rounded doubles cast to int64
        auto h   = static_cast<std::size_t>(1469598103934665603ull);
        auto mix = [&](long long x) {
            h ^= static_cast<std::size_t>(x);
            h *= static_cast<std::size_t>(1099511628211ull);
        };
        mix(static_cast<long long>(v[0] * 1e8));
        mix(static_cast<long long>(v[1] * 1e8));
        mix(static_cast<long long>(v[2] * 1e8));
        return h;
    }
};

struct VecEq {
    bool operator()(Vec3 const &a, Vec3 const &b) const noexcept {
        return std::llround(a[0] * 1e8) == std::llround(b[0] * 1e8) && std::llround(a[1] * 1e8) == std::llround(b[1] * 1e8) &&
               std::llround(a[2] * 1e8) == std::llround(b[2] * 1e8);
    }
};

inline Vec3 mul(const Mat3 &M, const Vec3 &p) {
    Vec3 r{0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            r[i] += static_cast<double>(M[i][j]) * p[j];
        }
    }
    return r;
}

std::vector<Mat3> permutation_matrices() {
    std::vector<Mat3>  mats;
    std::array<int, 3> idx{0, 1, 2};
    std::sort(idx.begin(), idx.end());
    do {
        Mat3 P{{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
        for (int i = 0; i < 3; ++i) {
            P[i][idx[i]] = 1;
        }
        mats.push_back(P);
    } while (std::next_permutation(idx.begin(), idx.end()));
    return mats;  // 6
}

std::vector<Mat3> reflection_matrices() {
    std::vector<Mat3> mats;
    for (int sx : {-1, 1}) {
        for (int sy : {-1, 1}) {
            for (int sz : {-1, 1}) {
                mats.push_back(Mat3{{{sx, 0, 0}, {0, sy, 0}, {0, 0, sz}}});
            }
        }
    }
    return mats;  // 8
}

std::vector<Mat3> symmetry_ops_full() {
    std::vector<Mat3> ops;
    auto              perms = permutation_matrices();
    auto              refls = reflection_matrices();
    ops.reserve(48);
    for (const auto &P : perms) {
        for (const auto &R : refls) {
            Mat3 M{};
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    int s = 0;
                    for (int k = 0; k < 3; ++k) {
                        s += R[i][k] * P[k][j];
                    }
                    M[i][j] = s;
                }
            }
            ops.push_back(M);
        }
    }
    // Remove duplicates (there should be exactly 48 unique ops)
    auto key = [](const Mat3 &M) {
        std::array<int, 9> flat{};
        int                t = 0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                flat[t++] = M[i][j];
            }
        }
        return flat;
    };
    std::sort(ops.begin(), ops.end(), [&](const Mat3 &A, const Mat3 &B) { return key(A) < key(B); });
    ops.erase(std::unique(ops.begin(), ops.end(), [&](const Mat3 &A, const Mat3 &B) { return key(A) == key(B); }), ops.end());
    if (ops.size() != 48) {
        std::cerr << "[warn] symmetry op count = " << ops.size() << " (expected 48)\n";
    }
    return ops;
}

std::vector<double> flatten_xyz(const std::vector<Vec3> &pts) {
    std::vector<double> coords;
    coords.reserve(pts.size() * 3);
    for (auto &p : pts) {
        coords.push_back(p[0]);
        coords.push_back(p[1]);
        coords.push_back(p[2]);
    }
    return coords;
}

std::vector<std::size_t> iota_tags(std::size_t n, std::size_t start = 1) {
    std::vector<std::size_t> v(n);
    for (std::size_t i = 0; i < n; ++i) {
        v[i] = start + i;
    }
    return v;
}

}  // namespace

int main(int argc, char **argv) try {
    int         nb_levels = 0;
    std::string outname   = "delaunay_mesh_refined.msh";
    bool        nogui     = false;
    if (argc > 1) {
        nb_levels = std::atoi(argv[1]);
    }
    if (argc > 2) {
        outname = argv[2];
    }
    if (argc > 3 && std::string_view(argv[3]) == std::string_view{"--nogui"}) {
        nogui = true;
    }

    // Initial six points (Gamma, L, X, K, W, U)
    std::vector<Vec3> initPts = {
        {0.0, 0.0, 0.0},    // Gamma
        {0.5, 0.5, 0.5},    // L
        {1.0, 0.0, 0.0},    // X
        {0.75, 0.75, 0.0},  // K
        {1.0, 0.5, 0.0},    // W
        {1.0, 0.25, 0.25}   // U
    };

    // -------- Stage 1: tetrahedralize these points and (optionally) refine --------
    gmsh::initialize();
    gmsh::option::setNumber("Mesh.Algorithm3D", 10);  // Delaunay (new)

    // Build a discrete volume from the initial points using gmsh::algorithm::tetrahedralize
    std::vector<double>      coords0 = flatten_xyz(initPts);
    std::vector<std::size_t> tets0;
    gmsh::algorithm::tetrahedralize(coords0, tets0);

    // Create a discrete volume entity and add nodes/elements
    const int                dim       = 3;
    const int                volTag    = gmsh::model::addDiscreteEntity(dim);
    std::vector<std::size_t> nodeTags0 = iota_tags(initPts.size(), 1);
    gmsh::model::mesh::addNodes(dim, volTag, nodeTags0, coords0);
    // Element type 4 = 4-node tetrahedron, node connectivity concatenated per element
    gmsh::model::mesh::addElementsByType(volTag, /*type*/ 4, /*elt tags*/ {}, tets0);

    for (int i = 0; i < nb_levels; ++i) {
        gmsh::model::mesh::refine();
    }

    // Gather all nodes of the refined mesh
    std::vector<std::size_t> nodeTags;
    std::vector<double>      nodeCoords;
    std::vector<double>      _param;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, _param);
    std::cout << "Nodes: " << nodeTags.size() << "\n";

    // -------- Stage 2: apply the 48 symmetry ops and unique the points --------
    std::vector<Vec3> nodes;
    nodes.reserve(nodeTags.size());
    for (size_t i = 0; i + 2 < nodeCoords.size(); i += 3) {
        nodes.push_back({nodeCoords[i + 0], nodeCoords[i + 1], nodeCoords[i + 2]});
    }

    auto                                     ops = symmetry_ops_full();
    std::unordered_set<Vec3, VecHash, VecEq> uniq;
    uniq.reserve(nodes.size() * ops.size());
    for (const auto &M : ops) {
        for (const auto &p : nodes) {
            auto q = mul(M, p);
            // No bounding box clipping as per the commented Python code
            uniq.insert({q[0], q[1], q[2]});
        }
    }

    std::vector<Vec3> symPts;
    symPts.reserve(uniq.size());
    for (const auto &p : uniq) {
        symPts.push_back(p);
    }

    // -------- Stage 3: create a new discrete volume from the symmetrized points --------
    gmsh::model::add("delaunay_mesh_refined");
    const int vol2 = gmsh::model::addDiscreteEntity(3);

    // Flatten and tetrahedralize the expanded point cloud
    std::vector<double>      coords = flatten_xyz(symPts);
    std::vector<std::size_t> tets;
    gmsh::algorithm::tetrahedralize(coords, tets);

    // Node tags 1..N
    std::vector<std::size_t> nodeTags2 = iota_tags(symPts.size(), 1);
    gmsh::model::mesh::addNodes(3, vol2, nodeTags2, coords);
    gmsh::model::mesh::addElementsByType(vol2, 4, {}, tets);

    // Save as MSH (set version if you prefer 2.2 like the Python script)
    // gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
    gmsh::write(outname);

    if (!nogui) {
        // Launch the GUI; close manually
        gmsh::fltk::run();
    }

    gmsh::finalize();
    return 0;
} catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    try {
        gmsh::finalize();
    } catch (...) {
    }
    return 1;
}
