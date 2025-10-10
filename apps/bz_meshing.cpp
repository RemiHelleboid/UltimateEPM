/**
 * @file bz_meshing.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2025-10-10
 * 
 * 
 */

#include <gmsh.h>
#include <tclap/CmdLine.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

using Mat3 = std::array<std::array<int, 3>, 3>;
using Vec3 = std::array<double, 3>;

// --------------------------- Small utils ---------------------------

struct VecHash {
    std::size_t operator()(Vec3 const &v) const noexcept {
        // FNV-1a on rounded doubles
        auto h   = static_cast<std::size_t>(1469598103934665603ull);
        auto mix = [&](long long x) {
            h ^= static_cast<std::size_t>(x);
            h *= static_cast<std::size_t>(1099511628211ull);
        };
        mix(static_cast<long long>(std::llround(v[0] * 1e8)));
        mix(static_cast<long long>(std::llround(v[1] * 1e8)));
        mix(static_cast<long long>(std::llround(v[2] * 1e8)));
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
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            r[i] += static_cast<double>(M[i][j]) * p[j];
    return r;
}

std::vector<Mat3> permutation_matrices() {
    std::vector<Mat3>  mats;
    std::array<int, 3> idx{0, 1, 2};
    std::sort(idx.begin(), idx.end());
    do {
        Mat3 P{{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
        for (int i = 0; i < 3; ++i)
            P[i][idx[i]] = 1;
        mats.push_back(P);
    } while (std::next_permutation(idx.begin(), idx.end()));
    return mats;  // 6
}

std::vector<Mat3> reflection_matrices() {
    std::vector<Mat3> mats;
    for (int sx : {-1, 1})
        for (int sy : {-1, 1})
            for (int sz : {-1, 1})
                mats.push_back(Mat3{{{sx, 0, 0}, {0, sy, 0}, {0, 0, sz}}});
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
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    int s = 0;
                    for (int k = 0; k < 3; ++k)
                        s += R[i][k] * P[k][j];
                    M[i][j] = s;
                }
            ops.push_back(M);
        }
    }
    // unique them
    auto key = [](const Mat3 &M) {
        std::array<int, 9> flat{};
        int                t = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                flat[t++] = M[i][j];
        return flat;
    };
    std::sort(ops.begin(), ops.end(), [&](const Mat3 &A, const Mat3 &B) { return key(A) < key(B); });
    ops.erase(std::unique(ops.begin(), ops.end(), [&](const Mat3 &A, const Mat3 &B) { return key(A) == key(B); }), ops.end());
    if (ops.size() != 48) std::cerr << "[warn] symmetry op count = " << ops.size() << " (expected 48)\n";
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
    for (std::size_t i = 0; i < n; ++i)
        v[i] = start + i;
    return v;
}

// --------------------------- Geometry (IBZ wedge) ---------------------------

struct EdgeTags {
    int GX{-1}, LK{-1}, LU{-1}, GL{-1};
};
struct PointTags {
    int L{-1}, Gamma{-1}, X{-1};
};

struct BuildResult {
    int       volumeTag{-1};
    EdgeTags  edges{};
    PointTags points{};
};

BuildResult build_ibz_wedge(double h, double gammaMesh) {
    using namespace gmsh::model::occ;
    // Points
    int Gamma = addPoint(0.0, 0.0, 0.0, h);
    int L     = addPoint(0.5, 0.5, 0.5, h * gammaMesh);
    int X     = addPoint(1.0, 0.0, 0.0, h * gammaMesh);
    int K     = addPoint(0.75, 0.75, 0.0, h);
    int W     = addPoint(1.0, 0.5, 0.0, h);
    int U     = addPoint(1.0, 0.25, 0.25, h);

    // KLUW (quad)
    int l1      = addLine(K, L);  // LK
    int l2      = addLine(L, U);  // LU
    int l3      = addLine(U, W);  // UW
    int l4      = addLine(W, K);  // WK
    int cl_quad = addCurveLoop({l1, l2, l3, l4});
    int s10     = addPlaneSurface({cl_quad});

    // UWX (tri)
    int l6      = addLine(W, X);  // WX
    int l7      = addLine(X, U);  // XU
    int cl_tri1 = addCurveLoop({l3, l6, l7});
    int s11     = addPlaneSurface({cl_tri1});

    // ΓLK (tri)
    int l8      = addLine(Gamma, L);  // ΓL
    int l10     = addLine(K, Gamma);  // KΓ
    int cl_tri2 = addCurveLoop({l8, -l1, l10});
    int s12     = addPlaneSurface({cl_tri2});

    // ΓLUX (quad)
    int l12      = addLine(X, Gamma);  // XΓ (Γ–X)
    int cl_quad2 = addCurveLoop({l8, l2, -l7, l12});
    int s13      = addPlaneSurface({cl_quad2});

    // ΓKWX (quad)
    int cl_quad3 = addCurveLoop({-l10, -l4, l6, l12});
    int s14      = addPlaneSurface({cl_quad3});

    int sl  = addSurfaceLoop({s10, s11, s12, s13, s14});
    int vol = addVolume({sl});

    BuildResult r;
    r.volumeTag = vol;
    r.edges     = EdgeTags{l12, l1, l2, l8};  // GX, LK, LU, GL
    r.points    = PointTags{L, Gamma, X};
    return r;
}

// --------------------------- Size fields ---------------------------

void set_distance_to_point_robust(int fieldId, int pointTag) {
    // Prefer binding to geometry vertex list:
    try {
        gmsh::model::mesh::field::setNumbers(fieldId, "PointsList", std::vector<double>{(double)pointTag});
        return;
    } catch (...) {
        // Fallback: bind via NodesList: create/get node at that vertex
        std::vector<std::size_t> nodeTags;
        std::vector<double>      nodeCoords, dummy;
        gmsh::model::mesh::getNodes(nodeTags, nodeCoords, dummy, 0, pointTag);
        if (nodeTags.empty()) {
            try {
                gmsh::model::mesh::embed(0, {(double)pointTag}, 3, 1);
            } catch (...) {
            }
            gmsh::model::mesh::getNodes(nodeTags, nodeCoords, dummy, 0, pointTag);
        }
        if (!nodeTags.empty()) {
            gmsh::model::mesh::field::setNumbers(fieldId, "NodesList", std::vector<double>{(double)nodeTags[0]});
        }
    }
}

struct MeshKnobs {
    // global
    double h{0.01};
    double meshGamma{1.0};
    // Δ cigar
    double delta_t0{0.85};
    double delta_axial{0.12};
    double delta_radial{0.03};
    double h_delta{0.003};
    bool   enable_tube{true};
    double tube_size_min_factor{0.6};
    double tube_rmin{0.02};
    double tube_rmax{0.06};
    // L
    bool   enable_L{true};
    double L_radius{0.06};
    double h_L{0.004};
    bool   enable_L_tubes{true};
    double L_tube_size_min_factor{0.7};
    double L_tube_rmin{0.02};
    double L_tube_rmax{0.05};
};

int apply_background_fields(const BuildResult &g, const MeshKnobs &k) {
    using namespace gmsh::model::mesh::field;

    // Make background field authoritative-ish (tweak to taste)
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 1);

    std::vector<int> fieldsToMin;

    // Δ cigar via MathEval ellipsoidal distance (normalized)
    int f_ell = add("MathEval");
    {
        double t0 = k.delta_t0;
        double ax = std::max(1e-9, k.delta_axial);
        double rr = std::max(1e-9, k.delta_radial);
        std::ostringstream oss;
        // r = sqrt(((x-t0)/ax)^2 + (y/rr)^2 + (z/rr)^2)
        oss << "sqrt(((x-" << std::setprecision(17) << t0 << ")/" << ax << ")^2 + "
            << "(y/" << rr << ")^2 + (z/" << rr << ")^2)";
        setString(f_ell, "F", oss.str());
    }

    int f_thr_delta = add("Threshold");
    setNumber(f_thr_delta, "InField", f_ell);
    setNumber(f_thr_delta, "SizeMin", k.h_delta);
    setNumber(f_thr_delta, "SizeMax", k.h);
    setNumber(f_thr_delta, "DistMin", 0.5);
    setNumber(f_thr_delta, "DistMax", 1.0);
    fieldsToMin.push_back(f_thr_delta);

    // Optional: gentle tube along Γ–X
    if (k.enable_tube && g.edges.GX > 0) {
        int f_dist_edge = add("Distance");
        setNumbers(f_dist_edge, "EdgesList", std::vector<double>{(double)std::abs(g.edges.GX)});

        int f_thr_edge = add("Threshold");
        setNumber(f_thr_edge, "InField", f_dist_edge);
        setNumber(f_thr_edge, "SizeMin", std::max(1e-9, k.tube_size_min_factor * k.h));
        setNumber(f_thr_edge, "SizeMax", k.h);
        setNumber(f_thr_edge, "DistMin", k.tube_rmin);
        setNumber(f_thr_edge, "DistMax", k.tube_rmax);
        fieldsToMin.push_back(f_thr_edge);
    }

    // L bubble + optional star tubes
    if (k.enable_L && g.points.L > 0) {
        int f_dist_L = add("Distance");
        set_distance_to_point_robust(f_dist_L, g.points.L);

        int f_thr_L = add("Threshold");
        setNumber(f_thr_L, "InField", f_dist_L);
        setNumber(f_thr_L, "SizeMin", k.h_L);
        setNumber(f_thr_L, "SizeMax", k.h);
        setNumber(f_thr_L, "DistMin", std::max(1e-6, 0.5 * k.L_radius));
        setNumber(f_thr_L, "DistMax", k.L_radius);
        fieldsToMin.push_back(f_thr_L);

        if (k.enable_L_tubes) {
            int f_dist_L_edges = add("Distance");
            setNumbers(f_dist_L_edges,
                       "EdgesList",
                       std::vector<double>{(double)std::abs(g.edges.LK), (double)std::abs(g.edges.LU), (double)std::abs(g.edges.GL)});

            int f_thr_L_edges = add("Threshold");
            setNumber(f_thr_L_edges, "InField", f_dist_L_edges);
            setNumber(f_thr_L_edges, "SizeMin", std::max(1e-9, k.L_tube_size_min_factor * k.h));
            setNumber(f_thr_L_edges, "SizeMax", k.h);
            setNumber(f_thr_L_edges, "DistMin", k.L_tube_rmin);
            setNumber(f_thr_L_edges, "DistMax", k.L_tube_rmax);
            fieldsToMin.push_back(f_thr_L_edges);
        }
    }

    // Global cap (version-proof)
    int f_const = add("MathEval");
    {
        std::ostringstream oss;
        oss << std::setprecision(17) << k.h;
        setString(f_const, "F", oss.str());
    }
    fieldsToMin.push_back(f_const);

    int f_min = add("Min");
    setNumbers(f_min, "FieldsList", std::vector<double>(fieldsToMin.begin(), fieldsToMin.end()));
    setAsBackgroundMesh(f_min);

    return f_min;
}

// --------------------------- Presets (mode + level) ---------------------------

enum class Mode { Conduction, Valence };

MeshKnobs apply_presets(Mode mode, int level, const MeshKnobs &base, bool no_L = true) {
    MeshKnobs k = base;  // start from user/base

    // Level 0 = coarsest, Level 4 = finest.
    struct Pack {
        double h, dax, drad, hD, Lrad, hL, tubeFac, LtubeFac;
    };

    static const Pack ladder[5] = {
        {0.1, 0.20, 0.1, 0.0050, 0.08, 0.0060, 0.8, 0.9},    // L0
        {0.1, 0.2, 0.1, 0.0025, 0.07, 0.0050, 0.7, 0.8},    // L1
        {0.05, 0.3, 0.2, 0.0025, 0.06, 0.0040, 0.6, 0.7},    // L2
        {0.01, 0.16, 0.025, 0.0022, 0.05, 0.0032, 0.55, 0.6},  // L3
        {0.005, 0.18, 0.020, 0.0016, 0.04, 0.0026, 0.5, 0.5}    // L4
    };
    const Pack &p            = ladder[std::clamp(level, 0, 4)];
    k.h                      = p.h;
    k.delta_axial            = p.dax;
    k.delta_radial           = p.drad;
    k.h_delta                = p.hD;
    k.L_radius               = p.Lrad;
    k.h_L                    = p.hL;
    k.tube_size_min_factor   = p.tubeFac;
    k.L_tube_size_min_factor = p.LtubeFac;


    // Mode nudges:
    if (mode == Mode::Conduction) {
        // emphasize Δ (Γ->X); slightly de-emphasize L
        k.delta_t0    = 0.85;
        k.enable_tube = (not no_L);
        k.enable_L    = (not no_L);
        k.h_L *= 1.2;  // a bit coarser at L
        k.L_radius *= 0.9;
    } else {
        // Valence: emphasize L; Δ still there but milder
        k.delta_t0 = 0.80;  // can shift the cigar slightly toward Γ
        k.h_delta *= 1.3;   // slightly coarser at Δ
        k.enable_L = false;
        k.h_L *= 0.8;  // finer at L
        k.L_radius *= 1.1;
    }
    k.enable_L_tubes = (not no_L);
    return k;
}

// --------------------------- Kmap file ---------------------------
void export_kmap(const std::string &filename, const std::vector<std::vector<std::size_t>>& maps_k_iwedge) {
    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "[error] could not open " << filename << " for writing\n";
        return;
    }
    ofs << "# Irreducible wedge to full BZ k-point map\n";
    ofs << "# Format: k_index_iwedge k_indices_BZ\n";
    for (std::size_t iw = 0; iw < maps_k_iwedge.size(); ++iw) {
        ofs << iw;
        for (std::size_t k : maps_k_iwedge[iw]) {
            ofs << " " << k;
        }
        ofs << "\n";
    }
    ofs.close();
    std::cout << "[info] wrote kmap to " << filename << "\n";
}

// --------------------------- Program ---------------------------

}  // namespace

int main(int argc, char **argv) try {
    // ----------- CLI -----------
    TCLAP::CmdLine cmd("FCC IBZ/BZ mesher (C++), with Δ-cigar and L refinement", ' ', "1.0");

    TCLAP::ValueArg<std::string> modeArg("m", "mode", "Band mode: conduction|valence", false, "conduction", "string");
    TCLAP::ValueArg<int>         lvlArg("l", "level", "Preset refinement level [0..4]", false, 2, "int");

    // Global mesh + gamma scaling
    TCLAP::ValueArg<double> hArg("", "mesh", "Global base size h", false, 0.05, "float");
    TCLAP::ValueArg<double> gammaArg("", "mesh-gamma", "Refinement factor at Γ point", false, 1.0, "float");

    // Δ cigar knobs
    TCLAP::ValueArg<double> dt0Arg("", "delta-t0", "Δ position along Γ->X in [0,1]", false, 0.85, "float");
    TCLAP::ValueArg<double> daxArg("", "delta-axial", "Δ semi-axis along Γ->X", false, 0.12, "float");
    TCLAP::ValueArg<double> dradArg("", "delta-radial", "Δ semi-axis across (y,z)", false, 0.03, "float");
    TCLAP::ValueArg<double> dhArg("", "delta-h", "Target size at Δ core", false, 0.003, "float");
    TCLAP::SwitchArg        tubeOn("", "tube", "Enable Γ–X tube", cmd, true);
    TCLAP::ValueArg<double> tubeFacArg("", "tube-size-min-factor", "Min size factor on Γ–X tube core (×h)", false, 0.6, "float");
    TCLAP::ValueArg<double> tubeRminArg("", "tube-rmin", "Γ–X tube inner radius", false, 0.02, "float");
    TCLAP::ValueArg<double> tubeRmaxArg("", "tube-rmax", "Γ–X tube outer radius", false, 0.06, "float");

    // L point
    TCLAP::SwitchArg        LOn("", "L", "Enable L refinement", cmd, true);
    TCLAP::ValueArg<double> LradArg("", "L-radius", "Refinement radius around L", false, 0.06, "float");
    TCLAP::ValueArg<double> LhArg("", "L-h", "Target size at L core", false, 0.004, "float");
    TCLAP::SwitchArg        LtubesOn("", "L-tubes", "Enable tubes along L–K, L–U, Γ–L", cmd, true);
    TCLAP::ValueArg<double> LtubeFacArg("", "L-tube-size-min-factor", "Min size factor on L-star tube cores (×h)", false, 0.7, "float");
    TCLAP::ValueArg<double> LtubeRminArg("", "L-tube-rmin", "L tubes inner radius", false, 0.02, "float");
    TCLAP::ValueArg<double> LtubeRmaxArg("", "L-tube-rmax", "L tubes outer radius", false, 0.05, "float");

    // Output, GUI
    TCLAP::ValueArg<std::string> outArg("o", "outfile", "Output mesh filename (.msh)", false, "ibz.msh", "string");
    TCLAP::SwitchArg             noGui("", "nogui", "Do not open GUI", cmd, false);

    cmd.add(modeArg);
    cmd.add(lvlArg);
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
    cmd.add(LtubeFacArg);
    cmd.add(LtubeRminArg);
    cmd.add(LtubeRmaxArg);
    cmd.add(outArg);

    cmd.parse(argc, argv);

    // ----------- Presets + overrides -----------
    Mode mode  = (modeArg.getValue() == "valence") ? Mode::Valence : Mode::Conduction;
    int  level = std::clamp(lvlArg.getValue(), 0, 4);

    MeshKnobs base;
    base.h                    = hArg.getValue();
    base.meshGamma            = gammaArg.getValue();
    base.delta_t0             = dt0Arg.getValue();
    base.delta_axial          = daxArg.getValue();
    base.delta_radial         = dradArg.getValue();
    base.h_delta              = dhArg.getValue();
    base.enable_tube          = tubeOn.getValue();
    base.tube_size_min_factor = tubeFacArg.getValue();
    base.tube_rmin            = tubeRminArg.getValue();
    base.tube_rmax            = tubeRmaxArg.getValue();

    base.enable_L               = LOn.getValue();
    base.L_radius               = LradArg.getValue();
    base.h_L                    = LhArg.getValue();
    base.enable_L_tubes         = LtubesOn.getValue();
    base.L_tube_size_min_factor = LtubeFacArg.getValue();
    base.L_tube_rmin            = LtubeRminArg.getValue();
    base.L_tube_rmax            = LtubeRmaxArg.getValue();

    MeshKnobs k = apply_presets(mode, level, base);

    // ----------- Gmsh -----------
    gmsh::initialize();
    gmsh::model::add("BZ_from_IBZ_cpp");
    gmsh::option::setNumber("Mesh.Algorithm3D", 1);  // Delaunay 3D

    // Build geometry
    auto g = build_ibz_wedge(k.h, k.meshGamma);
    gmsh::model::occ::synchronize();

    // Fields
    apply_background_fields(g, k);

    // Physical tagging (optional)
    gmsh::model::addPhysicalGroup(3, std::vector<int>{g.volumeTag}, 1);
    gmsh::model::setPhysicalName(3, 1, "IBZ_Wedge");

    // Mesh
    gmsh::model::mesh::generate(3);

    // Export nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double>      nodeCoords, parametric;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, parametric);
    std::cout << "Nodes in IBZ mesh: " << nodeTags.size() << "\n";

    // Collect into Vec3 list
    std::vector<Vec3> nodes;
    nodes.reserve(nodeTags.size());
    for (size_t i = 0; i + 2 < nodeCoords.size(); i += 3)
        nodes.push_back({nodeCoords[i], nodeCoords[i + 1], nodeCoords[i + 2]});

    // Symmetry expansion
    std::cout << "Expanding nodes to full BZ using symmetry...\n";
    std::cout << "Expanding nodes to full BZ using symmetry...\n";
    auto ops = symmetry_ops_full();

    const std::size_t                     nb_nodes = nodes.size();
    std::vector<std::vector<std::size_t>> sym_map(nb_nodes);
    for (auto &v : sym_map)
        v.reserve(ops.size());

    // Index map: point -> stable ID
    std::unordered_map<Vec3, std::size_t, VecHash, VecEq> index;
    index.reserve(nb_nodes * ops.size());

    // Output point cloud (stable order)
    std::vector<Vec3> symPts;
    symPts.reserve(nb_nodes * ops.size());

    auto get_id = [&](const Vec3 &q) -> std::size_t {
        auto it = index.find(q);
        if (it != index.end()) return it->second;
        std::size_t id = symPts.size();
        index.emplace(q, id);
        symPts.push_back(q);
        return id;
    };

    for (std::size_t i = 0; i < nb_nodes; ++i) {
        const auto &p = nodes[i];
        for (const auto &M : ops) {
            Vec3        q  = mul(M, p);
            std::size_t id = get_id(q);
            sym_map[i].push_back(id);
        }
    }

    // ... now 'symPts' is your deduplicated full-BZ point cloud in stable order,
    // and 'sym_map[i]' lists the IDs of the 48 symmetry images of original node i.

    std::cout << "Nodes in full BZ mesh (after symmetry expansion): " << sym_map.size() << "\n";
    export_kmap("kmap_ibz_to_bz.txt", sym_map);

    // std::vector<Vec3> symPts;
    // symPts.reserve(uniq.size());
    // for (const auto &p : uniq)
    //     symPts.push_back(p);

    // Build discrete BZ from point cloud
    gmsh::model::add("bz_from_ibz_full_symmetry_cpp");
    int                      vol2   = gmsh::model::addDiscreteEntity(3);
    std::vector<double>      coords = flatten_xyz(symPts);
    std::vector<std::size_t> tets;
    gmsh::algorithm::tetrahedralize(coords, tets);

    std::vector<std::size_t> nodeTags2 = iota_tags(symPts.size(), 1);
    gmsh::model::mesh::addNodes(3, vol2, nodeTags2, coords);
    gmsh::model::mesh::addElementsByType(vol2, 4, std::vector<std::size_t>{}, tets);

    // Write additional outputs
    std::string out  = outArg.getValue();
    gmsh::write(out);

    if (!noGui.getValue()) gmsh::fltk::run();
    gmsh::finalize();
    return 0;

} catch (const TCLAP::ArgException &e) {
    std::cerr << "TCLAP error: " << e.error() << " for arg " << e.argId() << "\n";
    try {
        gmsh::finalize();
    } catch (...) {
    }
    return 2;
} catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    try {
        gmsh::finalize();
    } catch (...) {
    }
    return 1;
}
