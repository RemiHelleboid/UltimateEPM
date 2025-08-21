#!/usr/bin/env python3
import argparse, sys, os
from collections import defaultdict
import math

def quantize_xyz(x, y, z, inv_tol):
    return (round(x*inv_tol), round(y*inv_tol), round(z*inv_tol))

def main():
    p = argparse.ArgumentParser(
        description="Deduplicate mesh (nodes/elements) and preserve NodeData views by remapping them to surviving nodes."
    )
    p.add_argument("-i","--input", required=True, help="Input .msh (mesh + NodeData views)")
    p.add_argument("-o","--output", required=True, help="Output cleaned .msh (mesh + remapped views)")
    p.add_argument("--msh", type=float, default=4.1, help="Output MSH version (4.1 or 2.2). Default 4.1")
    p.add_argument("--tol", type=float, default=1e-12, help="Coordinate tolerance to identify duplicates (in current units)")
    args = p.parse_args()

    try:
        import gmsh
    except Exception:
        print("ERROR: gmsh Python module not found.", file=sys.stderr)
        sys.exit(1)

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Verbosity", 2)

        if not os.path.isfile(args.input):
            print(f"ERROR: '{args.input}' not found.", file=sys.stderr); sys.exit(1)
        print(f"[gmsh] Open: {args.input}")
        gmsh.open(args.input)

        # -------- 1) Snapshot node coords and all views BEFORE dedup ----------
        old_tags, old_coords, _ = gmsh.model.mesh.getNodes()
        n_old = len(old_tags)
        tag2xyz = {}
        for i, t in enumerate(old_tags):
            x, y, z = old_coords[3*i:3*i+3]
            tag2xyz[t] = (x, y, z)

        views = []
        for vtag in gmsh.view.getTags():
            # name (optional)
            try:
                name = gmsh.option.getString(f"View[{gmsh.view.getIndex(vtag)}].Name")
            except Exception:
                name = f"view_{vtag}"

            typ = "NodeData"
            tags, data, time, numComp = [], [], 0.0, 1
            typ, tags, data, time, numComp = gmsh.view.getHomogeneousModelData(vtag, 0)
            if typ != "NodeData" or numComp != 1:
                # Skip non-scalar NodeData views (you can extend if needed)
                print(f"[warn] Skipping view {name} (type={typ}, numComp={numComp})")
                continue
            views.append({"name": name, "tags": tags, "data": data})

        print(f"[snapshot] Views captured: {len(views)} | nodes: {n_old}")

        # -------- 2) Deduplicate geometry ----------
        gmsh.model.mesh.removeDuplicateNodes()
        gmsh.model.mesh.removeDuplicateElements()

        # -------- 3) Build coordinate->newtag map AFTER dedup ----------
        new_tags, new_coords, _ = gmsh.model.mesh.getNodes()
        inv_tol = 1.0 / args.tol if args.tol > 0 else 1e12

        xyzkey2newtag = {}
        for i, t in enumerate(new_tags):
            x, y, z = new_coords[3*i:3*i+3]
            xyzkey2newtag[quantize_xyz(x, y, z, inv_tol)] = t

        # For writing views, we’ll keep original new_tags order
        newtag2index = {t: i for i, t in enumerate(new_tags)}

        # -------- 4) Rebuild views on surviving nodes ----------
        # Clear existing views to avoid duplicates
        for vtag in gmsh.view.getTags():
            gmsh.view.remove(vtag)

        # Get model name for addHomogeneousModelData
        model_name = ""
        try:
            model_name = gmsh.model.getCurrent()
        except Exception:
            pass

        for v in views:
            name = v["name"]
            orig_tags = v["tags"]
            data = v["data"]

            # Accumulate values per surviving new_tag (average if multiple map to same)
            acc = defaultdict(lambda: [0.0, 0])
            missing = 0
            for t, val in zip(orig_tags, data):
                x, y, z = tag2xyz.get(t, (None, None, None))
                if x is None:
                    missing += 1
                    continue
                key = quantize_xyz(x, y, z, inv_tol)
                new_t = xyzkey2newtag.get(key, None)
                if new_t is None:
                    missing += 1
                    continue
                acc[new_t][0] += val
                acc[new_t][1] += 1

            # Create arrays aligned to new_tags order (sparse → dense)
            out_vals = [0.0] * len(new_tags)
            out_tags = list(new_tags)  # store for clarity, though gmsh ignores provided tags for NodeData hom. model

            for new_t, (s, c) in acc.items():
                idx = newtag2index[new_t]
                out_vals[idx] = s / c

            vtag_new = gmsh.view.add(name)
            # Note: we pass new_tags (the current node tag list) and out_vals aligned to it
            gmsh.view.addHomogeneousModelData(vtag_new, 0, model_name, "NodeData", list(new_tags), out_vals)
            # Optional: hide in GUI
            gmsh.option.setNumber(f"View[{gmsh.view.getIndex(vtag_new)}].Visible", 0)

            if missing:
                print(f"[view '{name}'] remapped with {missing} samples that had no surviving node (dedup removed them).")

        # -------- 5) Write cleaned mesh + views ----------
        gmsh.option.setNumber("Mesh.MshFileVersion", args.msh)
        gmsh.write(args.output)
        print(f"[gmsh] Wrote: {args.output}")

        # Quick check
        gmsh.clear()
        gmsh.open(args.output)
        v_out = gmsh.view.getTags()
        print(f"[verify] Views in output: {len(v_out)}")

    finally:
        try:
            gmsh.finalize()
        except Exception:
            pass

if __name__ == "__main__":
    main()
