import sys
import argparse
import itertools
import numpy as np
import gmsh


# --------------------------- symmetry helpers ---------------------------

def _permutation_matrices():
    mats = []
    for perm in itertools.permutations(range(3)):
        P = np.zeros((3, 3), dtype=int)
        for i, j in enumerate(perm):
            P[i, j] = 1
        mats.append(P)
    return mats  # 6

def _reflection_matrices():
    mats = []
    for sx, sy, sz in itertools.product([-1, 1], repeat=3):
        mats.append(np.diag([sx, sy, sz]))
    return mats  # 8

def _symmetry_ops_full():
    ops, seen = [], set()
    for P in _permutation_matrices():
        for R in _reflection_matrices():
            M = R @ P
            key = tuple(M.flatten())
            if key not in seen:
                seen.add(key)
                ops.append(M)
    assert len(ops) == 48
    return ops

def _symmetry_ops_octant():
    ops = _permutation_matrices()
    assert len(ops) == 6
    return ops

def _is_identity(M: np.ndarray) -> bool:
    return np.all(M == np.eye(3, dtype=int))

def _apply_affine_to_copy(vol_tag: int, M: np.ndarray):
    new = gmsh.model.occ.copy([(3, vol_tag)])
    A = [
        M[0,0], M[0,1], M[0,2], 0,
        M[1,0], M[1,1], M[1,2], 0,
        M[2,0], M[2,1], M[2,2], 0,
        0,      0,      0,      1
    ]
    gmsh.model.occ.affineTransform(new, A)
    return new

def _get_all_vols():
    return gmsh.model.getEntities(3)  # list[(3, tag)]


# --------------------------- IBZ geometry ---------------------------

def build_ibz_wedge(h=0.01, gamma_mesh=1.0):
    # Points
    Gamma = gmsh.model.occ.addPoint(0.0, 0.0, 0.0, h)  # Γ
    L     = gmsh.model.occ.addPoint(0.5, 0.5, 0.5, h * gamma_mesh)  # L
    X     = gmsh.model.occ.addPoint(1.0, 0.0, 0.0, h * gamma_mesh)  # X
    K     = gmsh.model.occ.addPoint(0.75, 0.75, 0.0, h)  # K
    W     = gmsh.model.occ.addPoint(1.0, 0.5, 0.0, h)    # W
    U     = gmsh.model.occ.addPoint(1.0, 0.25, 0.25, h)  # U

    # KLUW (quad)
    l1 = gmsh.model.occ.addLine(K, L)   # LK
    l2 = gmsh.model.occ.addLine(L, U)   # LU
    l3 = gmsh.model.occ.addLine(U, W)   # UW
    l4 = gmsh.model.occ.addLine(W, K)   # WK
    s10 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])])

    # UWX (tri)
    l6 = gmsh.model.occ.addLine(W, X)   # WX
    l7 = gmsh.model.occ.addLine(X, U)   # XU
    s11 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l3, l6, l7])])

    # ΓLK (tri)
    l8  = gmsh.model.occ.addLine(Gamma, L)  # ΓL
    l10 = gmsh.model.occ.addLine(K, Gamma)  # KΓ
    s12 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l8, -l1, l10])])

    # ΓLUX (quad)
    l12 = gmsh.model.occ.addLine(X, Gamma)  # XΓ (Γ–X edge)
    s13 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l8, l2, -l7, l12])])

    # ΓKWX (quad)
    s14 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([-l10, -l4, l6, l12])])

    # Volume
    v21 = gmsh.model.occ.addVolume([gmsh.model.occ.addSurfaceLoop([s10, s11, s12, s13, s14])])

    # Important edge/point tags for refinement
    edges = {
        "GX": l12,
        "LK": l1,
        "LU": l2,
        "GL": l8,
    }
    points = {
        "L": L,
        "Gamma": Gamma,
        "X": X,
    }
    return v21, (s10, s11, s12, s13, s14), edges, points


# --------------------------- meshing ---------------------------

def build_mode_bz(h, outfile, gamma_mesh=1.0,
                  # Δ-valley refinement (ellipsoidal/cigar around t0 on Γ->X)
                  delta_t0=0.85, delta_axial=0.12, delta_radial=0.03, h_delta=0.003,
                  # tube along Γ–X to improve grading (optional)
                  tube_size_min_factor=0.6, tube_rmin=0.02, tube_rmax=0.06, enable_tube=True,
                  # L-point refinement
                  enable_L=True, L_radius=0.06, h_L=0.004,
                  # star tubes from L along LK, LU, and ΓL
                  enable_L_tubes=True, L_tube_size_min_factor=0.7, L_tube_rmin=0.02, L_tube_rmax=0.05):
    v21, faces, edges, points = build_ibz_wedge(h, gamma_mesh=gamma_mesh)
    gmsh.model.occ.synchronize()

    # Make background field authoritative-ish
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)

    fields_to_min = []

    # ---------- Δ-focused refinement (ELLIPSOID = cigar) ----------
    # Ellipsoidal "distance": r = sqrt(((x-t0)/ax)^2 + (y/rr)^2 + (z/rr)^2)
    # Threshold on r: fine for r <= 0.5, smooth to baseline by r = 1
    f_ell = gmsh.model.mesh.field.add("MathEval")
    t0 = float(delta_t0)
    ax = max(1e-9, float(delta_axial))
    rr = max(1e-9, float(delta_radial))
    expr = f"sqrt(((x-{t0})/{ax})^2 + (y/{rr})^2 + (z/{rr})^2)"
    gmsh.model.mesh.field.setString(f_ell, "F", expr)

    f_thr_delta = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_thr_delta, "InField", f_ell)
    gmsh.model.mesh.field.setNumber(f_thr_delta, "SizeMin", h_delta)  # fine at core
    gmsh.model.mesh.field.setNumber(f_thr_delta, "SizeMax", h)        # baseline outside
    gmsh.model.mesh.field.setNumber(f_thr_delta, "DistMin", 0.5)      # inner core (normalized)
    gmsh.model.mesh.field.setNumber(f_thr_delta, "DistMax", 1.0)      # outer boundary
    fields_to_min.append(f_thr_delta)

    # Optional: gentle tube along Γ–X (helps alignment/gradation)
    if enable_tube and "GX" in edges:
        f_dist_edge = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(f_dist_edge, "EdgesList", [abs(edges["GX"])])

        f_thr_edge = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(f_thr_edge, "InField", f_dist_edge)
        gmsh.model.mesh.field.setNumber(f_thr_edge, "SizeMin", max(1e-9, tube_size_min_factor * h))
        gmsh.model.mesh.field.setNumber(f_thr_edge, "SizeMax", h)
        gmsh.model.mesh.field.setNumber(f_thr_edge, "DistMin", tube_rmin)
        gmsh.model.mesh.field.setNumber(f_thr_edge, "DistMax", tube_rmax)
        fields_to_min.append(f_thr_edge)

    # ---------- L-point refinement (bubble at L) ----------
    if enable_L:
        f_dist_L = gmsh.model.mesh.field.add("Distance")
        # L is a geometry vertex; PointsList is typically supported
        try:
            gmsh.model.mesh.field.setNumbers(f_dist_L, "PointsList", [points["L"]])
        except Exception:
            # Fallback via node on the vertex
            node_tags, _, _ = gmsh.model.mesh.getNodes(0, points["L"])
            if not node_tags:
                try:
                    gmsh.model.mesh.embed(0, [points["L"]], 3, 1)
                except Exception:
                    pass
                node_tags, _, _ = gmsh.model.mesh.getNodes(0, points["L"])
            if node_tags:
                gmsh.model.mesh.field.setNumbers(f_dist_L, "NodesList", [int(node_tags[0])])

        f_thr_L = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(f_thr_L, "InField", f_dist_L)
        gmsh.model.mesh.field.setNumber(f_thr_L, "SizeMin", h_L)
        gmsh.model.mesh.field.setNumber(f_thr_L, "SizeMax", h)
        gmsh.model.mesh.field.setNumber(f_thr_L, "DistMin", max(1e-6, 0.5 * L_radius))
        gmsh.model.mesh.field.setNumber(f_thr_L, "DistMax", L_radius)
        fields_to_min.append(f_thr_L)

        if enable_L_tubes:
            # Star tubes along LK, LU, and ΓL
            edge_list = [abs(edges["LK"]), abs(edges["LU"]), abs(edges["GL"])]
            f_dist_L_edges = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(f_dist_L_edges, "EdgesList", edge_list)

            f_thr_L_edges = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(f_thr_L_edges, "InField", f_dist_L_edges)
            gmsh.model.mesh.field.setNumber(f_thr_L_edges, "SizeMin", max(1e-9, L_tube_size_min_factor * h))
            gmsh.model.mesh.field.setNumber(f_thr_L_edges, "SizeMax", h)
            gmsh.model.mesh.field.setNumber(f_thr_L_edges, "DistMin", L_tube_rmin)
            gmsh.model.mesh.field.setNumber(f_thr_L_edges, "DistMax", L_tube_rmax)
            fields_to_min.append(f_thr_L_edges)

    # Global cap via MathEval (version-proof)
    f_const = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(f_const, "F", str(h))
    fields_to_min.append(f_const)

    # Combine and apply as background mesh
    f_min = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(f_min, "FieldsList", fields_to_min)
    gmsh.model.mesh.field.setAsBackgroundMesh(f_min)

    # ---------- tagging + mesh ----------
    gmsh.model.addPhysicalGroup(3, [v21], 1)
    gmsh.model.setPhysicalName(3, 1, "IBZ_Wedge")

    gmsh.model.mesh.generate(3)
    gmsh.write(outfile or "ibz.msh")

    # ---------- export nodes + symmetry expansion ----------
    nodes = gmsh.model.mesh.getNodes()
    print(f"Number of nodes: {len(nodes[0])}")
    list_points = np.array(nodes[1]).reshape(-1, 3)
    np.savetxt("ibz_nodes.csv", list_points, delimiter=",")

    AllTransforms = _symmetry_ops_full()
    print(f"Number of symmetry operations: {len(AllTransforms)}")
    List_new_nodes = []
    for op in AllTransforms:
        for p in list_points:
            p_new = op @ p
            List_new_nodes.append(tuple(np.round(p_new, 8)))
    List_new_nodes = list(set(List_new_nodes))
    print(f"Number of unique new nodes: {len(List_new_nodes)}")
    np.savetxt("ibz_nodes_full_symmetry.csv", np.array(List_new_nodes),
               delimiter=",", header="X,Y,Z", comments='')

    gmsh.model.add("bz_from_ibz_full_symmetry")

    vol2 = gmsh.model.addDiscreteEntity(3)
    N = len(List_new_nodes)
    list_coords = []
    for node in List_new_nodes:
        list_coords.extend([node[0], node[1], node[2]])
    points_arr = np.array(list_coords, dtype=float)
    tets = gmsh.algorithm.tetrahedralize(points_arr)
    gmsh.model.mesh.addNodes(3, vol2, range(1, N+1), points_arr)
    gmsh.model.mesh.addElementsByType(vol2, 4, [], tets)

    gmsh.write(outfile or "bz_full_symmetry.msh")
    gmsh.write("bz_full_symmetry.vtk")


# --------------------------- main ---------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description="FCC IBZ / BZ builder with Δ-cigar + L refinement.")
    ap.add_argument("--mesh", type=float, default=0.01)
    ap.add_argument("--outfile", type=str, default="")
    ap.add_argument("--nopopup", action="store_true")
    ap.add_argument("--bbox-tol", type=float, default=1e-12)
    ap.add_argument("--mesh-gamma", type=float, default=1.0,
                    help="refinement factor at Γ point (default 1.0 = uniform).")

    # Δ valley (along Γ->X), ellipsoidal
    ap.add_argument("--delta-t0", type=float, default=0.85, help="Δ position along Γ->X in [0,1].")
    ap.add_argument("--delta-axial", type=float, default=0.12, help="Semi-axis length along Γ->X.")
    ap.add_argument("--delta-radial", type=float, default=0.03, help="Semi-axis radius across (y,z).")
    ap.add_argument("--delta-h", type=float, default=0.003, help="Target size at Δ core.")

    # Optional tube along Γ–X
    ap.add_argument("--tube", action="store_true", default=True,
                    help="Enable gentle refinement tube along Γ–X.")
    ap.add_argument("--tube-size-min-factor", type=float, default=0.6,
                    help="Min size factor (× --mesh) on Γ–X tube core.")
    ap.add_argument("--tube-rmin", type=float, default=0.02, help="Tube inner radius.")
    ap.add_argument("--tube-rmax", type=float, default=0.06, help="Tube outer radius.")

    # L point
    ap.add_argument("--L", dest="enable_L", action="store_true", default=True,
                    help="Enable refinement around L.")
    ap.add_argument("--L-radius", type=float, default=0.06, help="Refinement radius around L.")
    ap.add_argument("--L-h", type=float, default=0.004, help="Target size at L core.")
    ap.add_argument("--L-tubes", dest="enable_L_tubes", action="store_true", default=True,
                    help="Enable tubes along L–K, L–U, and Γ–L.")
    ap.add_argument("--L-tube-size-min-factor", type=float, default=0.7,
                    help="Min size factor (× --mesh) on L-star tube cores.")
    ap.add_argument("--L-tube-rmin", type=float, default=0.02, help="L tubes inner radius.")
    ap.add_argument("--L-tube-rmax", type=float, default=0.05, help="L tubes outer radius.")

    args = ap.parse_args(argv)

    gmsh.initialize()
    gmsh.model.add("BZ_from_IBZ")
    # 3D mesh algorithm (1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    try:
        build_mode_bz(
            args.mesh, args.outfile, args.mesh_gamma,
            delta_t0=args.delta_t0,
            delta_axial=args.delta_axial, delta_radial=args.delta_radial,
            h_delta=args.delta_h,
            tube_size_min_factor=args.tube_size_min_factor,
            tube_rmin=args.tube_rmin, tube_rmax=args.tube_rmax,
            enable_tube=args.tube,
            enable_L=args.enable_L, L_radius=args.L_radius, h_L=args.L_h,
            enable_L_tubes=args.enable_L_tubes,
            L_tube_size_min_factor=args.L_tube_size_min_factor,
            L_tube_rmin=args.L_tube_rmin, L_tube_rmax=args.L_tube_rmax
        )
    finally:
        gmsh.finalize()
    return 0


if __name__ == "__main__":
    sys.exit(main())
