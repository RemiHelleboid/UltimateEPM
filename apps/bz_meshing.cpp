/**
 * @file bz_meshing.cpp
 * @brief FCC IBZ/BZ mesher with symmetry-exact full-BZ construction.
 *
 * Key properties:
 * - Meshes the IBZ wedge with Gmsh.
 * - Expands the IBZ mesh to the full BZ by symmetry-copying BOTH nodes and tetrahedra.
 * - Keeps the first N full-BZ nodes equal to the IBZ node ordering, so k-star files remain compatible.
 * - Avoids re-tetrahedralizing the symmetric point cloud.
 * - Avoids post-expansion mesh optimization, which would otherwise break exact symmetry.
 */

#include <fmt/format.h>
#include <gmsh.h>
#include <tclap/CmdLine.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

using Mat3 = std::array<std::array<int, 3>, 3>;
using Vec3 = std::array<double, 3>;
using Tet4 = std::array<std::size_t, 4>;

constexpr double kCoordQuant       = 1e12;
constexpr double kDegenerateTetTol = 1e-18;

double canonical_zero(double x) noexcept { return (std::abs(x) < 1e-15) ? 0.0 : x; }

long long quantize(double x) noexcept { return std::llround(canonical_zero(x) * kCoordQuant); }

struct VecHash {
    std::size_t operator()(const Vec3& v) const noexcept {
        std::size_t h   = static_cast<std::size_t>(1469598103934665603ull);
        auto        mix = [&](long long x) {
            h ^= static_cast<std::size_t>(x);
            h *= static_cast<std::size_t>(1099511628211ull);
        };
        mix(quantize(v[0]));
        mix(quantize(v[1]));
        mix(quantize(v[2]));
        return h;
    }
};

struct VecEq {
    bool operator()(const Vec3& a, const Vec3& b) const noexcept {
        return quantize(a[0]) == quantize(b[0]) && quantize(a[1]) == quantize(b[1]) && quantize(a[2]) == quantize(b[2]);
    }
};

struct TetKeyHash {
    std::size_t operator()(const Tet4& t) const noexcept {
        std::size_t h   = static_cast<std::size_t>(1469598103934665603ull);
        auto        mix = [&](std::size_t x) {
            h ^= x;
            h *= static_cast<std::size_t>(1099511628211ull);
        };
        mix(t[0]);
        mix(t[1]);
        mix(t[2]);
        mix(t[3]);
        return h;
    }
};

Vec3 operator+(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }

Vec3 operator-(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }

Vec3 cross(const Vec3& a, const Vec3& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }

double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

double signed_six_volume(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) { return dot(b - a, cross(c - a, d - a)); }

Vec3 mul(const Mat3& M, const Vec3& p) {
    Vec3 r{0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            r[i] += static_cast<double>(M[i][j]) * p[j];
        }
        r[i] = canonical_zero(r[i]);
    }
    return r;
}

bool is_identity(const Mat3& M) {
    static constexpr Mat3 I{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    return M == I;
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
    return mats;
}

std::vector<Mat3> reflection_matrices() {
    std::vector<Mat3> mats;
    mats.reserve(8);
    for (int sx : {-1, 1}) {
        for (int sy : {-1, 1}) {
            for (int sz : {-1, 1}) {
                mats.push_back(Mat3{{{sx, 0, 0}, {0, sy, 0}, {0, 0, sz}}});
            }
        }
    }
    return mats;
}

std::vector<Mat3> symmetry_ops_full() {
    std::vector<Mat3> ops;
    const auto        perms = permutation_matrices();
    const auto        refls = reflection_matrices();
    ops.reserve(48);

    for (const auto& P : perms) {
        for (const auto& R : refls) {
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

    auto key = [](const Mat3& M) {
        std::array<int, 9> flat{};
        int                t = 0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                flat[t++] = M[i][j];
            }
        }
        return flat;
    };

    std::sort(ops.begin(), ops.end(), [&](const Mat3& A, const Mat3& B) { return key(A) < key(B); });
    ops.erase(std::unique(ops.begin(), ops.end(), [&](const Mat3& A, const Mat3& B) { return key(A) == key(B); }), ops.end());

    if (ops.size() != 48) {
        throw std::runtime_error("symmetry_ops_full: expected 48 cubic operations");
    }

    auto itI = std::find_if(ops.begin(), ops.end(), [](const Mat3& M) { return is_identity(M); });
    if (itI == ops.end()) {
        throw std::runtime_error("symmetry_ops_full: identity operation not found");
    }
    std::rotate(ops.begin(), itI, std::next(itI));

    return ops;
}

std::vector<double> flatten_xyz(const std::vector<Vec3>& pts) {
    std::vector<double> coords;
    coords.reserve(pts.size() * 3);
    for (const auto& p : pts) {
        coords.push_back(p[0]);
        coords.push_back(p[1]);
        coords.push_back(p[2]);
    }
    return coords;
}

std::vector<std::size_t> iota_tags(std::size_t n, std::size_t start = 1) {
    std::vector<std::size_t> v(n);
    std::iota(v.begin(), v.end(), start);
    return v;
}

Tet4 sorted_key(Tet4 t) {
    std::sort(t.begin(), t.end());
    return t;
}

bool has_duplicate_vertices(const Tet4& t) {
    Tet4 s = sorted_key(t);
    return (s[0] == s[1]) || (s[1] == s[2]) || (s[2] == s[3]);
}

Tet4 orient_positive(const Tet4& t, const std::vector<Vec3>& pts) {
    Tet4         out = t;
    const double v6  = signed_six_volume(pts[out[0]], pts[out[1]], pts[out[2]], pts[out[3]]);

    if (std::abs(v6) < kDegenerateTetTol) {
        throw std::runtime_error("orient_positive: degenerate tetrahedron encountered");
    }
    if (v6 < 0.0) {
        std::swap(out[0], out[1]);
    }
    return out;
}

struct EdgeTags {
    int GX{-1};
    int LK{-1};
    int LU{-1};
    int GL{-1};
};

struct PointTags {
    int L{-1};
    int Gamma{-1};
    int X{-1};
};

struct BuildResult {
    int       volumeTag{-1};
    EdgeTags  edges{};
    PointTags points{};
};

BuildResult build_ibz_wedge(double h, double gammaMesh) {
    using namespace gmsh::model::occ;

    const int Gamma = addPoint(0.0, 0.0, 0.0, h * gammaMesh);
    const int L     = addPoint(0.5, 0.5, 0.5, h * gammaMesh);
    const int X     = addPoint(1.0, 0.0, 0.0, h * gammaMesh);
    const int K     = addPoint(0.75, 0.75, 0.0, h);
    const int W     = addPoint(1.0, 0.5, 0.0, h);
    const int U     = addPoint(1.0, 0.25, 0.25, h);

    const int lLK    = addLine(K, L);
    const int lLU    = addLine(L, U);
    const int lUW    = addLine(U, W);
    const int lWK    = addLine(W, K);
    const int clKLUW = addCurveLoop({lLK, lLU, lUW, lWK});
    const int sKLUW  = addPlaneSurface({clKLUW});

    const int lWX   = addLine(W, X);
    const int lXU   = addLine(X, U);
    const int clUWX = addCurveLoop({lUW, lWX, lXU});
    const int sUWX  = addPlaneSurface({clUWX});

    const int lGL   = addLine(Gamma, L);
    const int lKG   = addLine(K, Gamma);
    const int clGLK = addCurveLoop({lGL, -lLK, lKG});
    const int sGLK  = addPlaneSurface({clGLK});

    const int lXG    = addLine(X, Gamma);
    const int clGLUX = addCurveLoop({lGL, lLU, -lXU, lXG});
    const int sGLUX  = addPlaneSurface({clGLUX});

    const int clGKWX = addCurveLoop({-lKG, -lWK, lWX, lXG});
    const int sGKWX  = addPlaneSurface({clGKWX});

    const int sl  = addSurfaceLoop({sKLUW, sUWX, sGLK, sGLUX, sGKWX});
    const int vol = addVolume({sl});

    BuildResult out;
    out.volumeTag = vol;
    out.edges     = EdgeTags{lXG, lLK, lLU, lGL};
    out.points    = PointTags{L, Gamma, X};
    return out;
}

void apply_uniform_mesh(double h) {
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeMin", h);
    gmsh::option::setNumber("Mesh.MeshSizeMax", h);

    using namespace gmsh::model::mesh::field;
    const int          f = add("MathEval");
    std::ostringstream oss;
    oss << std::setprecision(17) << h;
    setString(f, "F", oss.str());
    setAsBackgroundMesh(f);
}

void set_distance_to_point(int fieldId, int pointTag) {
    gmsh::model::mesh::field::setNumbers(fieldId, "PointsList", std::vector<double>{static_cast<double>(pointTag)});
}

struct MeshKnobs {
    double h{0.01};
    double meshGamma{1.0};

    double delta_t0{0.85};
    double delta_axial{0.12};
    double delta_radial{0.03};
    double h_delta{0.003};

    bool   enable_tube{true};
    double tube_size_min_factor{0.6};
    double tube_rmin{0.02};
    double tube_rmax{0.06};

    bool   enable_L{true};
    double L_radius{0.06};
    double h_L{0.004};

    bool   enable_L_tubes{true};
    double L_tube_size_min_factor{0.7};
    double L_tube_rmin{0.02};
    double L_tube_rmax{0.05};
};

int apply_background_fields(const BuildResult& g, const MeshKnobs& k) {
    using namespace gmsh::model::mesh::field;

    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 1);

    std::vector<int> fieldsToMin;

    const int fEll = add("MathEval");
    {
        const double t0 = k.delta_t0;
        const double ax = std::max(1e-12, k.delta_axial);
        const double rr = std::max(1e-12, k.delta_radial);

        std::ostringstream oss;
        oss << "sqrt(((x-" << std::setprecision(17) << t0 << ")/" << ax << ")^2 + "
            << "(y/" << rr << ")^2 + (z/" << rr << ")^2)";
        setString(fEll, "F", oss.str());
    }

    const int fThrDelta = add("Threshold");
    setNumber(fThrDelta, "InField", fEll);
    setNumber(fThrDelta, "SizeMin", k.h_delta);
    setNumber(fThrDelta, "SizeMax", k.h);
    setNumber(fThrDelta, "DistMin", 0.5);
    setNumber(fThrDelta, "DistMax", 1.0);
    fieldsToMin.push_back(fThrDelta);

    if (k.enable_tube && g.edges.GX > 0) {
        const int fDistEdge = add("Distance");
        setNumbers(fDistEdge, "EdgesList", std::vector<double>{static_cast<double>(std::abs(g.edges.GX))});

        const int fThrEdge = add("Threshold");
        setNumber(fThrEdge, "InField", fDistEdge);
        setNumber(fThrEdge, "SizeMin", std::max(1e-12, k.tube_size_min_factor * k.h));
        setNumber(fThrEdge, "SizeMax", k.h);
        setNumber(fThrEdge, "DistMin", k.tube_rmin);
        setNumber(fThrEdge, "DistMax", k.tube_rmax);
        fieldsToMin.push_back(fThrEdge);
    }

    if (k.enable_L && g.points.L > 0) {
        const int fDistL = add("Distance");
        set_distance_to_point(fDistL, g.points.L);

        const int fThrL = add("Threshold");
        setNumber(fThrL, "InField", fDistL);
        setNumber(fThrL, "SizeMin", k.h_L);
        setNumber(fThrL, "SizeMax", k.h);
        setNumber(fThrL, "DistMin", std::max(1e-9, 0.5 * k.L_radius));
        setNumber(fThrL, "DistMax", k.L_radius);
        fieldsToMin.push_back(fThrL);

        if (k.enable_L_tubes) {
            const int fDistLEdges = add("Distance");
            setNumbers(fDistLEdges,
                       "EdgesList",
                       std::vector<double>{static_cast<double>(std::abs(g.edges.LK)),
                                           static_cast<double>(std::abs(g.edges.LU)),
                                           static_cast<double>(std::abs(g.edges.GL))});

            const int fThrLEdges = add("Threshold");
            setNumber(fThrLEdges, "InField", fDistLEdges);
            setNumber(fThrLEdges, "SizeMin", std::max(1e-12, k.L_tube_size_min_factor * k.h));
            setNumber(fThrLEdges, "SizeMax", k.h);
            setNumber(fThrLEdges, "DistMin", k.L_tube_rmin);
            setNumber(fThrLEdges, "DistMax", k.L_tube_rmax);
            fieldsToMin.push_back(fThrLEdges);
        }
    }

    const int fConst = add("MathEval");
    {
        std::ostringstream oss;
        oss << std::setprecision(17) << k.h;
        setString(fConst, "F", oss.str());
    }
    fieldsToMin.push_back(fConst);

    const int fMin = add("Min");
    setNumbers(fMin, "FieldsList", std::vector<double>(fieldsToMin.begin(), fieldsToMin.end()));
    setAsBackgroundMesh(fMin);

    return fMin;
}

enum class Mode { Conduction, Valence };

MeshKnobs preset_knobs(Mode mode, int level) {
    struct Pack {
        double h;
        double dax;
        double drad;
        double hD;
        double Lrad;
        double hL;
        double tubeFac;
        double LtubeFac;
    };

    static const Pack ladder[5] = {{0.050, 0.20, 0.10, 0.0050, 0.080, 0.0100, 0.80, 0.90},
                                   {0.100, 0.20, 0.10, 0.0025, 0.070, 0.0050, 0.70, 0.80},
                                   {0.050, 0.30, 0.20, 0.0050, 0.060, 0.0080, 0.60, 0.70},
                                   {0.010, 0.16, 0.025, 0.0022, 0.050, 0.0032, 0.55, 0.60},
                                   {0.005, 0.18, 0.020, 0.0016, 0.040, 0.0026, 0.50, 0.50}};

    const Pack& p = ladder[std::clamp(level, 0, 4)];

    MeshKnobs k;
    k.h                      = p.h;
    k.meshGamma              = 1.0;
    k.delta_axial            = p.dax;
    k.delta_radial           = p.drad;
    k.h_delta                = p.hD;
    k.L_radius               = p.Lrad;
    k.h_L                    = p.hL;
    k.tube_size_min_factor   = p.tubeFac;
    k.L_tube_size_min_factor = p.LtubeFac;
    k.tube_rmin              = 0.02;
    k.tube_rmax              = 0.06;
    k.L_tube_rmin            = 0.02;
    k.L_tube_rmax            = 0.05;

    if (mode == Mode::Conduction) {
        k.delta_t0       = 0.85;
        k.enable_tube    = true;
        k.enable_L       = true;
        k.enable_L_tubes = true;
        k.h_L *= 1.20;
        k.L_radius *= 0.90;
    } else {
        k.delta_t0       = 0.80;
        k.enable_tube    = true;
        k.enable_L       = true;
        k.enable_L_tubes = true;
        k.h_delta *= 1.25;
        k.h_L *= 0.75;
        k.L_radius *= 1.15;
        k.tube_size_min_factor = std::min(0.95, k.tube_size_min_factor + 0.10);
    }

    return k;
}

template <typename TArg, typename TValue>
void apply_override(TValue& dst, const TArg& arg) {
    if (arg.isSet()) {
        dst = arg.getValue();
    }
}

void apply_toggle(bool& dst, const TCLAP::SwitchArg& onArg, const TCLAP::SwitchArg& offArg, const std::string& name) {
    if (onArg.getValue() && offArg.getValue()) {
        throw std::runtime_error("conflicting options for " + name);
    }
    if (onArg.getValue()) {
        dst = true;
    }
    if (offArg.getValue()) {
        dst = false;
    }
}

struct IbzMesh {
    std::vector<Vec3> nodes;
    std::vector<Tet4> tets;
};

IbzMesh extract_ibz_mesh(int volumeTag) {
    gmsh::model::mesh::removeDuplicateNodes();
    gmsh::model::mesh::removeDuplicateElements();

    std::vector<std::size_t> nodeTags;
    std::vector<double>      nodeCoords;
    std::vector<double>      nodeParams;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams);

    if (nodeCoords.size() != 3 * nodeTags.size()) {
        throw std::runtime_error("extract_ibz_mesh: inconsistent node coordinate array");
    }

    IbzMesh mesh;
    mesh.nodes.reserve(nodeTags.size());

    std::unordered_map<std::size_t, std::size_t> localFromTag;
    localFromTag.reserve(nodeTags.size());

    for (std::size_t i = 0; i < nodeTags.size(); ++i) {
        localFromTag.emplace(nodeTags[i], i);
        mesh.nodes.push_back(
            {canonical_zero(nodeCoords[3 * i + 0]), canonical_zero(nodeCoords[3 * i + 1]), canonical_zero(nodeCoords[3 * i + 2])});
    }

    std::vector<int>                      elemTypes;
    std::vector<std::vector<std::size_t>> elemTags;
    std::vector<std::vector<std::size_t>> elemNodeTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, 3, volumeTag);

    bool foundTet4 = false;
    for (std::size_t it = 0; it < elemTypes.size(); ++it) {
        if (elemTypes[it] != 4) {
            continue;
        }
        foundTet4 = true;

        const auto& conn = elemNodeTags[it];
        if (conn.size() % 4 != 0) {
            throw std::runtime_error("extract_ibz_mesh: tetra connectivity is not a multiple of 4");
        }

        const std::size_t nTet = conn.size() / 4;
        mesh.tets.reserve(nTet);

        for (std::size_t e = 0; e < nTet; ++e) {
            Tet4 t{};
            for (int k = 0; k < 4; ++k) {
                const std::size_t tag = conn[4 * e + k];
                auto              jt  = localFromTag.find(tag);
                if (jt == localFromTag.end()) {
                    throw std::runtime_error("extract_ibz_mesh: tetra references unknown node tag");
                }
                t[k] = jt->second;
            }

            if (has_duplicate_vertices(t)) {
                continue;
            }

            const double v6 = signed_six_volume(mesh.nodes[t[0]], mesh.nodes[t[1]], mesh.nodes[t[2]], mesh.nodes[t[3]]);
            if (std::abs(v6) < kDegenerateTetTol) {
                continue;
            }

            mesh.tets.push_back(orient_positive(t, mesh.nodes));
        }
    }

    if (!foundTet4) {
        throw std::runtime_error("extract_ibz_mesh: no 4-node tetrahedra found in IBZ mesh");
    }

    if (mesh.nodes.empty() || mesh.tets.empty()) {
        throw std::runtime_error("extract_ibz_mesh: empty IBZ mesh");
    }

    return mesh;
}

struct ExpandedMesh {
    std::vector<Vec3>                     nodes;
    std::vector<Tet4>                     tets;
    std::vector<std::vector<std::size_t>> orbitIds;
};

ExpandedMesh expand_mesh_by_symmetry(const IbzMesh& ibz, const std::vector<Mat3>& ops) {
    ExpandedMesh out;

    const std::size_t nIbz = ibz.nodes.size();
    out.nodes.reserve(nIbz * ops.size());
    out.orbitIds.resize(nIbz);

    std::unordered_map<Vec3, std::size_t, VecHash, VecEq> nodeIndex;
    nodeIndex.reserve(nIbz * ops.size());

    auto get_or_insert_node = [&](const Vec3& p) -> std::size_t {
        auto it = nodeIndex.find(p);
        if (it != nodeIndex.end()) {
            return it->second;
        }
        const std::size_t id = out.nodes.size();
        nodeIndex.emplace(p, id);
        out.nodes.push_back(p);
        return id;
    };

    for (std::size_t i = 0; i < nIbz; ++i) {
        get_or_insert_node(ibz.nodes[i]);
    }

    for (std::size_t i = 0; i < nIbz; ++i) {
        std::unordered_set<std::size_t> seen;
        seen.reserve(ops.size());

        for (const auto& M : ops) {
            const std::size_t id = get_or_insert_node(mul(M, ibz.nodes[i]));
            if (seen.insert(id).second) {
                out.orbitIds[i].push_back(id);
            }
        }

        if (out.orbitIds[i].empty()) {
            throw std::runtime_error("expand_mesh_by_symmetry: empty orbit");
        }
    }

    std::unordered_set<Tet4, TetKeyHash> seenTets;
    seenTets.reserve(ibz.tets.size() * ops.size());

    out.tets.reserve(ibz.tets.size() * ops.size());

    for (const auto& M : ops) {
        for (const auto& tIbz : ibz.tets) {
            Tet4 tFull{};
            for (int k = 0; k < 4; ++k) {
                tFull[k] = get_or_insert_node(mul(M, ibz.nodes[tIbz[k]]));
            }

            if (has_duplicate_vertices(tFull)) {
                continue;
            }

            Tet4 key = sorted_key(tFull);
            if (!seenTets.insert(key).second) {
                continue;
            }

            const double v6 = signed_six_volume(out.nodes[tFull[0]], out.nodes[tFull[1]], out.nodes[tFull[2]], out.nodes[tFull[3]]);
            if (std::abs(v6) < kDegenerateTetTol) {
                continue;
            }

            out.tets.push_back(orient_positive(tFull, out.nodes));
        }
    }

    if (out.nodes.size() < nIbz) {
        throw std::runtime_error("expand_mesh_by_symmetry: internal node ordering failure");
    }

    return out;
}

void export_kstar(const std::string& filename, const std::vector<std::vector<std::size_t>>& orbitIds) {
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("could not open " + filename + " for writing");
    }

    out << "# iw m ids...\n";
    for (std::size_t iw = 0; iw < orbitIds.size(); ++iw) {
        out << iw << " " << orbitIds[iw].size();
        for (std::size_t id : orbitIds[iw]) {
            out << " " << id;
        }
        out << "\n";
    }
}

double mesh_volume(const std::vector<Vec3>& pts, const std::vector<Tet4>& tets) {
    double sum = 0.0;
    for (const auto& t : tets) {
        sum += std::abs(signed_six_volume(pts[t[0]], pts[t[1]], pts[t[2]], pts[t[3]])) / 6.0;
    }
    return sum;
}

void write_full_bz_discrete_mesh(const ExpandedMesh& fullMesh, const std::string& modelName) {
    gmsh::model::add(modelName);
    gmsh::model::setCurrent(modelName);

    const int volTag = gmsh::model::addDiscreteEntity(3);

    const std::vector<double>      coords   = flatten_xyz(fullMesh.nodes);
    const std::vector<std::size_t> nodeTags = iota_tags(fullMesh.nodes.size(), 1);

    gmsh::model::mesh::addNodes(3, volTag, nodeTags, coords);

    std::vector<std::size_t> elemTags = iota_tags(fullMesh.tets.size(), 1);
    std::vector<std::size_t> conn;
    conn.reserve(4 * fullMesh.tets.size());
    for (const auto& t : fullMesh.tets) {
        conn.push_back(t[0] + 1);
        conn.push_back(t[1] + 1);
        conn.push_back(t[2] + 1);
        conn.push_back(t[3] + 1);
    }

    gmsh::model::mesh::addElementsByType(volTag, 4, elemTags, conn);

    gmsh::model::addPhysicalGroup(3, std::vector<int>{volTag}, 1);
    gmsh::model::setPhysicalName(3, 1, "Full_BZ");
}

}  // namespace

int main(int argc, char** argv) try {
    TCLAP::CmdLine cmd("FCC IBZ/BZ mesher with symmetry-exact full-BZ construction", ' ', "2.0");

    TCLAP::ValueArg<std::string> modeArg("m", "mode", "Band mode: conduction|valence", false, "conduction", "string");
    TCLAP::ValueArg<int>         levelArg("l", "level", "Preset refinement level [0..4]", false, 2, "int");
    TCLAP::SwitchArg             uniformArg("", "uniform", "Use strictly uniform mesh", cmd, false);

    TCLAP::ValueArg<double> hArg("", "mesh", "Global base size h", false, 0.0, "float");
    TCLAP::ValueArg<double> gammaArg("", "mesh-gamma", "Refinement factor at Gamma point", false, 0.0, "float");

    TCLAP::ValueArg<double> dt0Arg("", "delta-t0", "Delta position along Gamma->X in [0,1]", false, 0.0, "float");
    TCLAP::ValueArg<double> daxArg("", "delta-axial", "Delta semi-axis along Gamma->X", false, 0.0, "float");
    TCLAP::ValueArg<double> dradArg("", "delta-radial", "Delta semi-axis transverse to Gamma->X", false, 0.0, "float");
    TCLAP::ValueArg<double> dhArg("", "delta-h", "Target size at Delta core", false, 0.0, "float");

    TCLAP::SwitchArg        tubeOn("", "tube", "Force-enable Gamma-X tube", cmd, false);
    TCLAP::SwitchArg        tubeOff("", "no-tube", "Disable Gamma-X tube", cmd, false);
    TCLAP::ValueArg<double> tubeFacArg("", "tube-size-min-factor", "Min size factor on Gamma-X tube core (x h)", false, 0.0, "float");
    TCLAP::ValueArg<double> tubeRminArg("", "tube-rmin", "Gamma-X tube inner radius", false, 0.0, "float");
    TCLAP::ValueArg<double> tubeRmaxArg("", "tube-rmax", "Gamma-X tube outer radius", false, 0.0, "float");

    TCLAP::SwitchArg        LOn("", "L", "Force-enable L refinement", cmd, false);
    TCLAP::SwitchArg        LOff("", "no-L", "Disable L refinement", cmd, false);
    TCLAP::ValueArg<double> LradArg("", "L-radius", "Refinement radius around L", false, 0.0, "float");
    TCLAP::ValueArg<double> LhArg("", "L-h", "Target size at L core", false, 0.0, "float");

    TCLAP::SwitchArg        LTubesOn("", "L-tubes", "Force-enable L-star tubes", cmd, false);
    TCLAP::SwitchArg        LTubesOff("", "no-L-tubes", "Disable L-star tubes", cmd, false);
    TCLAP::ValueArg<double> LTubesFacArg("", "L-tube-size-min-factor", "Min size factor on L-star tube cores (x h)", false, 0.0, "float");
    TCLAP::ValueArg<double> LTubesRminArg("", "L-tube-rmin", "L-star tube inner radius", false, 0.0, "float");
    TCLAP::ValueArg<double> LTubesRmaxArg("", "L-tube-rmax", "L-star tube outer radius", false, 0.0, "float");

    TCLAP::ValueArg<std::string> outArg("o", "outfile", "Output full-BZ mesh filename (.msh)", false, "bz.msh", "string");
    TCLAP::SwitchArg             noGuiArg("", "nogui", "Do not open the GUI", cmd, false);

    cmd.add(modeArg);
    cmd.add(levelArg);

    cmd.add(hArg);
    cmd.add(gammaArg);

    cmd.add(dt0Arg);
    cmd.add(daxArg);
    cmd.add(dradArg);
    cmd.add(dhArg);

    cmd.add(tubeFacArg);
    cmd.add(tubeRminArg);
    cmd.add(tubeRmaxArg);

    cmd.add(LradArg);
    cmd.add(LhArg);

    cmd.add(LTubesFacArg);
    cmd.add(LTubesRminArg);
    cmd.add(LTubesRmaxArg);

    cmd.add(outArg);
    cmd.parse(argc, argv);

    Mode mode;
    if (modeArg.getValue() == "conduction") {
        mode = Mode::Conduction;
    } else if (modeArg.getValue() == "valence") {
        mode = Mode::Valence;
    } else {
        throw std::runtime_error("invalid mode: expected 'conduction' or 'valence'");
    }

    const int level = std::clamp(levelArg.getValue(), 0, 4);

    MeshKnobs k = preset_knobs(mode, level);

    apply_override(k.h, hArg);
    apply_override(k.meshGamma, gammaArg);

    apply_override(k.delta_t0, dt0Arg);
    apply_override(k.delta_axial, daxArg);
    apply_override(k.delta_radial, dradArg);
    apply_override(k.h_delta, dhArg);

    apply_override(k.tube_size_min_factor, tubeFacArg);
    apply_override(k.tube_rmin, tubeRminArg);
    apply_override(k.tube_rmax, tubeRmaxArg);

    apply_override(k.L_radius, LradArg);
    apply_override(k.h_L, LhArg);

    apply_override(k.L_tube_size_min_factor, LTubesFacArg);
    apply_override(k.L_tube_rmin, LTubesRminArg);
    apply_override(k.L_tube_rmax, LTubesRmaxArg);

    apply_toggle(k.enable_tube, tubeOn, tubeOff, "Gamma-X tube");
    apply_toggle(k.enable_L, LOn, LOff, "L refinement");
    apply_toggle(k.enable_L_tubes, LTubesOn, LTubesOff, "L-star tubes");

    if (!k.enable_L) {
        k.enable_L_tubes = false;
    }

    if (k.h <= 0.0) {
        throw std::runtime_error("mesh size must be > 0");
    }
    if (k.meshGamma <= 0.0) {
        throw std::runtime_error("mesh-gamma must be > 0");
    }
    if (k.delta_axial <= 0.0 || k.delta_radial <= 0.0 || k.h_delta <= 0.0) {
        throw std::runtime_error("Delta refinement parameters must be > 0");
    }
    if (k.enable_L && (k.L_radius <= 0.0 || k.h_L <= 0.0)) {
        throw std::runtime_error("L refinement parameters must be > 0");
    }

    gmsh::initialize();

    gmsh::option::setNumber("Mesh.Binary", 1);
    gmsh::option::setNumber("Mesh.Algorithm3D", 4);
    gmsh::option::setNumber("General.Verbosity", 10);

    const std::string ibzModelName = "IBZ_Wedge_Model";
    gmsh::model::add(ibzModelName);
    gmsh::model::setCurrent(ibzModelName);

    const BuildResult geom = build_ibz_wedge(k.h, k.meshGamma);
    gmsh::model::occ::synchronize();

    if (uniformArg.getValue()) {
        apply_uniform_mesh(k.h);
        std::cout << "Meshing mode: uniform\n";
    } else {
        apply_background_fields(geom, k);
        std::cout << "Meshing mode: adaptive\n";
    }

    gmsh::model::addPhysicalGroup(3, std::vector<int>{geom.volumeTag}, 1);
    gmsh::model::setPhysicalName(3, 1, "IBZ_Wedge");

    gmsh::model::mesh::generate(3);
    gmsh::model::mesh::removeDuplicateNodes();
    gmsh::model::mesh::removeDuplicateElements();

    const IbzMesh ibz = extract_ibz_mesh(geom.volumeTag);

    std::cout << "IBZ nodes: " << ibz.nodes.size() << "\n";
    std::cout << "IBZ tetrahedra: " << ibz.tets.size() << "\n";

    const std::vector<Mat3> ops  = symmetry_ops_full();
    const ExpandedMesh      full = expand_mesh_by_symmetry(ibz, ops);

    std::cout << "Full-BZ nodes: " << full.nodes.size() << "\n";
    std::cout << "Full-BZ tetrahedra: " << full.tets.size() << "\n";

    const double volIBZ = mesh_volume(ibz.nodes, ibz.tets);
    const double volBZ  = mesh_volume(full.nodes, full.tets);
    std::cout << std::setprecision(16);
    std::cout << "IBZ volume = " << volIBZ << "\n";
    std::cout << "Full-BZ volume = " << volBZ << "\n";
    if (volIBZ > 0.0) {
        std::cout << "Volume ratio Full/IBZ = " << (volBZ / volIBZ) << "\n";
    }

    const std::string outMesh       = outArg.getValue();
    const std::string stem          = std::filesystem::path(outMesh).stem().string();
    const std::string kstarFilename = fmt::format("{}_kstar_ibz_to_bz.txt", stem);
    export_kstar(kstarFilename, full.orbitIds);

    const std::string fullModelName = "Full_BZ_Model";
    write_full_bz_discrete_mesh(full, fullModelName);

    if (std::filesystem::exists(outMesh)) {
        std::filesystem::remove(outMesh);
    }
    gmsh::write(outMesh);

    std::cout << "Wrote mesh: " << outMesh << "\n";
    std::cout << "Wrote k-star map: " << kstarFilename << "\n";

    if (!noGuiArg.getValue()) {
        gmsh::fltk::run();
    }

    gmsh::finalize();
    return 0;

} catch (const TCLAP::ArgException& e) {
    std::cerr << "TCLAP error: " << e.error() << " for arg " << e.argId() << "\n";
    try {
        gmsh::finalize();
    } catch (...) {
    }
    return 2;
} catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    try {
        gmsh::finalize();
    } catch (...) {
    }
    return 1;
}