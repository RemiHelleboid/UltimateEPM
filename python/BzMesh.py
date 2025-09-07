#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FCC IBZ → BZ builder with strict union option.

Modes:
  --mode {iw,octant,full}   : single IBZ; 6 wedges in +octant; 48 wedges (full BZ)
  --mesh <h>                : target size (default 0.01)
  --outfile <path>          : output .msh (defaults per mode)
  --nopopup                 : skip GUI
  --union                   : boolean-union all wedges into ONE watertight volume
  --octant-filter           : (octant mode) build 48, keep only +octant via bbox
  --bbox-tol <tol>          : bbox tol for octant-filter (default 1e-12)

The IBZ is constructed EXACTLY like your .geo (same points/lines/loops/surfaces).
"""

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


# --------------------------- EXACT IBZ geometry (matches your .geo) ---------------------------

def build_ibz_wedge(h=0.01, gamma_mesh=1.0):
    # Points
    Gamma = gmsh.model.occ.addPoint(0.0, 0.0, 0.0, h * gamma_mesh) # Γ
    L = gmsh.model.occ.addPoint(0.5, 0.5, 0.5, h)     # L
    X = gmsh.model.occ.addPoint(1.0, 0.0, 0.0, h)     # X
    K = gmsh.model.occ.addPoint(0.75, 0.75, 0.0, h)   # K
    W = gmsh.model.occ.addPoint(1.0, 0.5, 0.0, h)     # W
    U = gmsh.model.occ.addPoint(1.0, 0.25, 0.25, h)   # U

    # KLUW (quad)
    l1 = gmsh.model.occ.addLine(K, L)
    l2 = gmsh.model.occ.addLine(L, U)
    l3 = gmsh.model.occ.addLine(U, W)
    l4 = gmsh.model.occ.addLine(W, K)
    s10 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])])

    # UWX (tri)
    l6 = gmsh.model.occ.addLine(W, X)
    l7 = gmsh.model.occ.addLine(X, U)
    s11 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l3, l6, l7])])

    # ΓLK (tri)
    l8  = gmsh.model.occ.addLine(Gamma, L)
    l10 = gmsh.model.occ.addLine(K, Gamma)
    s12 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l8, -l1, l10])])

    # ΓLUX (quad)
    l12 = gmsh.model.occ.addLine(X, Gamma)
    s13 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l8, l2, -l7, l12])])

    # ΓKWX (quad)
    s14 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([-l10, -l4, l6, l12])])

    # Volume
    v21 = gmsh.model.occ.addVolume([gmsh.model.occ.addSurfaceLoop([s10, s11, s12, s13, s14])])
    return v21, (s10, s11, s12, s13, s14)


# --------------------------- strict union utility ---------------------------

def strict_union_all_volumes(batch=8, max_passes=6):
    """
    Robustly reduce ALL current OCC volumes to ONE by repeated fragment/fuse.
    - Works in batches to avoid heavy single-call unions.
    - Stops early if a single volume remains.
    """
    def fuse_two_sets(A, B):
        if not A or not B:
            return A or B
        out, _ = gmsh.model.occ.fuse(A, B)
        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()
        return out

    for _ in range(max_passes):
        vols = _get_all_vols()
        if len(vols) <= 1:
            break

        # Pre-fragment once to split overlapping shells and align faces
        gmsh.model.occ.fragment(vols[:1], vols[1:])
        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()

        # Batch fuse
        vols = _get_all_vols()
        if len(vols) <= 1:
            break

        # Group into chunks and fuse progressively
        chunks = [vols[i:i+batch] for i in range(0, len(vols), batch)]
        # Fuse within each chunk
        fused_chunks = []
        for ch in chunks:
            if len(ch) == 1:
                fused_chunks.append(ch)
                continue
            out, _ = gmsh.model.occ.fuse(ch[:1], ch[1:])
            gmsh.model.occ.removeAllDuplicates()
            gmsh.model.occ.synchronize()
            fused_chunks.append(out)

        # Fuse chunks together
        current = fused_chunks[0]
        for ch in fused_chunks[1:]:
            current = fuse_two_sets(current, ch)

        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()

        if len(_get_all_vols()) <= 1:
            break

    # Final assert (best effort)
    vols = _get_all_vols()
    if len(vols) > 1:
        # As a last resort, keep the volume with the largest bbox diagonal, remove the rest
        sizes = []
        for (d, t) in vols:
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(d, t)
            diag2 = (xmax-xmin)**2 + (ymax-ymin)**2 + (zmax-zmin)**2
            sizes.append((diag2, (d, t)))
        sizes.sort(reverse=True)
        keep = sizes[0][1]
        remove = [x[1] for x in sizes[1:]]
        gmsh.model.occ.remove(remove, recursive=True)
        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()


# --------------------------- build modes ---------------------------

def build_mode_iw(h, outfile, nopopup, do_union, gamma_mesh=1.0):
    v21, faces = build_ibz_wedge(h, gamma_mesh=gamma_mesh)
    gmsh.model.occ.synchronize()

    if do_union:
        # Single wedge already one volume; no-op but keep API symmetry
        strict_union_all_volumes()

    gmsh.model.addPhysicalGroup(3, [v21], 1)
    gmsh.model.setPhysicalName(3, 1, "IBZ_Wedge")
    # gmsh.model.addPhysicalGroup(2, list(faces), 2)
    # gmsh.model.setPhysicalName(2, 2, "IBZ_Faces")

    gmsh.model.mesh.generate(3)
    gmsh.write(outfile or "ibz.msh")
    if not nopopup:
        gmsh.fltk.run()

def build_mode_octant(h, outfile, nopopup, do_union, octant_filter=False, bbox_tol=1e-12):
    v21, faces = build_ibz_wedge(h)
    gmsh.model.occ.synchronize()

    if octant_filter:
        # Build 48 and keep only +octant via bbox
        for M in _symmetry_ops_full():
            if _is_identity(M):
                continue
            _apply_affine_to_copy(v21, M)
        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()
        kept = []
        for (dim, tag) in _get_all_vols():
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
            if xmin >= -bbox_tol and ymin >= -bbox_tol and zmin >= -bbox_tol:
                kept.append((dim, tag))
        to_remove = [ent for ent in _get_all_vols() if ent not in kept]
        if to_remove:
            gmsh.model.occ.remove(to_remove, recursive=True)
        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()
    else:
        # Only 6 permutations (no reflections)
        for M in _symmetry_ops_octant():
            if _is_identity(M):
                continue
            _apply_affine_to_copy(v21, M)
        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()

    if do_union:
        strict_union_all_volumes()

    vol_tags = [tag for (_, tag) in _get_all_vols()]
    gmsh.model.addPhysicalGroup(3, vol_tags, 1)
    gmsh.model.setPhysicalName(3, 1, "BZ_positive_octant")
    # gmsh.model.addPhysicalGroup(2, list(faces), 2)
    # gmsh.model.setPhysicalName(2, 2, "IBZ_Faces")

    gmsh.model.mesh.generate(3)
    gmsh.write(outfile or ("bz_positive_octant_union.msh" if do_union else "bz_positive_octant_from_ibz.msh"))
    if not nopopup:
        gmsh.fltk.run()

def build_mode_full(h, outfile, nopopup, do_union):
    v21, faces = build_ibz_wedge(h)
    gmsh.model.occ.synchronize()

    for M in _symmetry_ops_full():
        if _is_identity(M):
            continue
        _apply_affine_to_copy(v21, M)

    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    if do_union:
        strict_union_all_volumes()

    vol_tags = [tag for (_, tag) in _get_all_vols()]
    gmsh.model.addPhysicalGroup(3, vol_tags, 1)
    gmsh.model.setPhysicalName(3, 1, "BZ_full")
    # gmsh.model.addPhysicalGroup(2, list(faces), 2)
    # gmsh.model.setPhysicalName(2, 2, "IBZ_Faces")

    gmsh.model.mesh.generate(3)
    gmsh.write(outfile or ("full_bz_union.msh" if do_union else "full_bz_from_ibz.msh"))
    if not nopopup:
        gmsh.fltk.run()


# --------------------------- main ---------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description="FCC IBZ / BZ builder with strict union.")
    ap.add_argument("--mode", choices=["iw", "octant", "full"], default="iw")
    ap.add_argument("--mesh", type=float, default=0.01)
    ap.add_argument("--outfile", type=str, default="")
    ap.add_argument("--nopopup", action="store_true")
    ap.add_argument("--union", action="store_true", help="Boolean-union wedges into ONE volume.")
    ap.add_argument("--octant-filter", action="store_true", default=True,
                    help="(octant) build 48 then keep only +octant via bbox.")
    ap.add_argument("--bbox-tol", type=float, default=1e-12)
    ap.add_argument("--mesh-gamma", type=float, default=1.0,
                    help="refinement factor at Γ point (default 1.0 = uniform).")
    args = ap.parse_args(argv)

    gmsh.initialize()
    gmsh.model.add("BZ_from_IBZ")
    # gmsh.option.setNumber("Mesh.CharacteristicLengthMin", args.mesh)
    # gmsh.option.setNumber("Mesh.CharacteristicLengthMax", args.mesh)

    # 3D mesh algorithm (1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT)
    gmsh.option.setNumber("Mesh.Algorithm3D", 7)

    try:
        if args.mode == "iw":
            build_mode_iw(args.mesh, args.outfile, args.nopopup, args.union, args.mesh_gamma)
        elif args.mode == "octant":
            build_mode_octant(args.mesh, args.outfile, args.nopopup, args.union,
                              octant_filter=args.octant_filter, bbox_tol=args.bbox_tol)
        else:
            build_mode_full(args.mesh, args.outfile, args.nopopup, args.union)
    finally:
        gmsh.finalize()
    return 0

if __name__ == "__main__":
    main()

