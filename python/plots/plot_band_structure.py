#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations
import argparse
import json
import csv
from pathlib import Path
from typing import Dict, List, Sequence, Tuple, Any, Optional

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


# ---------- IO ----------

def load_band_data(filename: Path) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    CSV with optional header lines starting with '#', e.g.:
      # Material Si
      # NBands 12
      # Nonlocal Yes
      # Path LGXK
      # Symmetry points: 0 163 351 499

    Returns (bands: (nb, nk)), info dict.
    """
    material: str = "material"
    nb_bands: Optional[int] = None
    nonlocal_corr: Optional[bool] = None
    path_str: Optional[str] = None
    sym_idx: Optional[List[int]] = None

    header_lines: List[str] = []
    with filename.open("r", encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            header_lines.append(line.strip())

    for line in header_lines:
        s = line[1:].lstrip()  # drop '#'
        if s.startswith("Material "):
            material = s.split(" ", 1)[1].strip()
        elif s.startswith("NBands "):
            try:
                nb_bands = int(s.split(" ", 1)[1].strip())
            except Exception:
                nb_bands = None
        elif s.startswith("Nonlocal "):
            nonlocal_corr = s.split(" ", 1)[1].strip().lower() in {"yes", "true", "1"}
        elif s.startswith("Path "):
            path_str = s.split(" ", 1)[1].strip()
        elif s.lower().startswith("symmetry points"):
            after = s.split(":", 1)[1] if ":" in s else s.split(" ", 2)[-1]
            sym_idx = [int(x) for x in after.strip().split() if x.strip().lstrip("-").isdigit()]

    # Load numeric data; skip detected header lines
    data = np.loadtxt(filename, delimiter=",", comments="#", skiprows=len(header_lines)+1)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    info = {
        "material": material,
        "nb_bands": nb_bands,
        "nonlocal": nonlocal_corr,
        "path": path_str,                 # e.g. "LGXK"
        "symmetry_points": sym_idx,       # e.g. [0, 163, 351, 499]
        "header_lines": len(header_lines)
    }
    return data.T, info  # (nb, nk), info


# ---------- Metrics ----------

def estimate_gap_components(bands: np.ndarray, eps: float = 1e-1) -> Dict[str, Any]:
    """
    Simple VBM/CBM finder: values <= eps -> valence, >= eps -> conduction.
    """
    nb, nk = bands.shape
    valence_max = np.full(nk, -np.inf, dtype=float)
    valence_band_idx = np.full(nk, -1, dtype=int)
    conduction_min = np.full(nk, np.inf, dtype=float)
    conduction_band_idx = np.full(nk, -1, dtype=int)

    for b in range(nb):
        E = bands[b]
        mv = E <= eps
        if np.any(mv):
            better = (E > valence_max) & mv
            valence_max[better] = E[better]
            valence_band_idx[better] = b
        mc = E >= eps
        if np.any(mc):
            better = (E < conduction_min) & mc
            conduction_min[better] = E[better]
            conduction_band_idx[better] = b

    vbm_k = int(np.nanargmax(valence_max))
    cbm_k = int(np.nanargmin(conduction_min))
    vbm_E = float(valence_max[vbm_k])
    cbm_E = float(conduction_min[cbm_k])
    fundamental_gap = cbm_E - vbm_E

    Eg_k = conduction_min - valence_max
    min_direct_k = int(np.nanargmin(Eg_k))
    min_direct_gap = float(Eg_k[min_direct_k])

    return {
        "VBM": {"E_eV": float(vbm_E), "band": int(valence_band_idx[vbm_k]), "k_index": vbm_k},
        "CBM": {"E_eV": float(cbm_E), "band": int(conduction_band_idx[cbm_k]), "k_index": cbm_k},
        "fundamental_gap_eV": float(fundamental_gap),
        "is_direct": bool(vbm_k == cbm_k),
        "direct_gap_min_eV": float(min_direct_gap),
        "direct_gap_min_k_index": min_direct_k,
        "valence_max_per_k_eV": valence_max.tolist(),
        "conduction_min_per_k_eV": conduction_min.tolist(),
    }


def bandwidths(bands: np.ndarray) -> List[float]:
    return (bands.max(axis=1) - bands.min(axis=1)).astype(float).tolist()


# ---------- Symmetry dumps ----------

def dump_all_bands_at_symmetry_points(
    bands: np.ndarray,
    labels_tex: Sequence[str],
    sym_indices: Sequence[int],
) -> Tuple[Dict[str, Any], List[Tuple[str, int, int, float]], List[Tuple[str, int, List[float]]]]:
    nk = bands.shape[1]
    ticks = np.clip(np.asarray(sym_indices, dtype=int), 0, nk - 1)

    out_json: Dict[str, Any] = {}
    rows_long: List[Tuple[str, int, int, float]] = []
    rows_wide: List[Tuple[str, int, List[float]]] = []

    for lab, kidx in zip(labels_tex, ticks):
        col = bands[:, kidx]
        vmax = float(col[col <= 0.0].max()) if np.any(col <= 0.0) else float('nan')
        cmin = float(col[col >= 0.0].min()) if np.any(col >= 0.0) else float('nan')
        dgap = (cmin - vmax) if np.isfinite(vmax) and np.isfinite(cmin) else float('nan')
        out_json[str(lab)] = {
            "k_index": int(kidx),
            "bands_eV": col.astype(float).tolist(),
            "valence_max_eV": vmax,
            "conduction_min_eV": cmin,
            "direct_gap_eV": float(dgap),
        }
        for bi, e in enumerate(col):
            rows_long.append((str(lab), int(kidx), int(bi), float(e)))
        rows_wide.append((str(lab), int(kidx), col.astype(float).tolist()))
    return out_json, rows_long, rows_wide


# ---------- Plot ----------

def plot_and_save_figs(
    bands_list: List[np.ndarray],
    labels_tex: Sequence[str],
    material: str,
    out_dir: Path,
    nbbands: int,
    show_plot: bool,
    ylim: Optional[Tuple[float, float]],
    zoom_pad: float,
    sym_indices: Sequence[int],
) -> Tuple[int, float]:
    mpl.rcParams['figure.figsize'] = [3.8, 3.0]
    mpl.rcParams['axes.grid'] = True
    mpl.rcParams['grid.alpha'] = 0.35
    mpl.rcParams['grid.linewidth'] = 0.4

    fig, ax = plt.subplots()
    linestyles = ["-", "--", "-.", ":"]
    nk_ref = bands_list[0].shape[1]
    last_gap = 0.0

    for i, bands_all in enumerate(bands_list):
        nb = min(nbbands, bands_all.shape[0])
        bands = bands_all[:nb, :]
        ls = linestyles[i % len(linestyles)]
        for band in bands:
            ax.plot(band, ls=ls, lw=0.6, color="k")
        comps = estimate_gap_components(bands)
        last_gap = comps["fundamental_gap_eV"]

    ax.set_title(f"Band structure of {material}")
    ax.set_xlabel(r"$\mathbf{k}$")
    ax.set_ylabel("Energy (eV)")
    ax.axhline(0.0, lw=0.6, ls=":", color="k", alpha=0.6)

    ticks = np.clip(np.asarray(sym_indices, dtype=int), 0, nk_ref - 1)
    ax.set_xticks(ticks, labels=labels_tex)
    ax.set_xticks([], minor=True)
    for x in ticks:
        ax.axvline(x, lw=0.4, ls="--", alpha=0.35)

    if ylim is not None:
        ax.set_ylim(*ylim)
    fig.tight_layout()
    fig.savefig(out_dir / f"band_structure_{material}.png", dpi=300)

    ax.set_ylim(-4.0, max(0.0, last_gap) + zoom_pad)
    fig.tight_layout()
    fig.savefig(out_dir / f"band_structure_{material}_zoom.png", dpi=300)

    if show_plot:
        plt.show()

    return nk_ref, last_gap


# ---------- CLI ----------

def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description="Band plotter + metrics + symmetry-point dump (header-driven).")
    ap.add_argument("-f", "--file", dest="band_files", nargs="+", required=True,
                    help="Band-structure CSV files (columns are bands; header allowed).")
    ap.add_argument("-b", "--nbbands", dest="nb_bands", type=int, default=10,
                    help="Number of bands to analyze/plot.")
    ap.add_argument("-o", "--outputdir", dest="output_dir", default="./",
                    help="Directory to save figures and reports.")
    ap.add_argument("--fermi", type=float, default=0.0,
                    help="Energy shift (eV) subtracted from all bands.")
    ap.add_argument("--ylim", nargs=2, type=float, metavar=("YMIN", "YMAX"),
                    help="Y-axis limits for the full plot.")
    ap.add_argument("--zoom-pad", type=float, default=4.0,
                    help="Padding above CBM for the zoom figure.")
    ap.add_argument("--no-plot", dest="show_plot", action="store_false",
                    help="Do not display the plot window.")
    ap.set_defaults(show_plot=True)
    ap.add_argument("--report", default=None,
                    help="Basename for metrics reports (without extension). Default: band_metrics_<material>")
    args = ap.parse_args(argv)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    files = [Path(f) for f in args.band_files]

    # First file provides metadata (material, labels, indices)
    arr0, info0 = load_band_data(files[0])

    # Require header path and symmetry indices; this script is header-driven by design.
    if not info0.get("path") or not info0.get("symmetry_points"):
        raise ValueError(
            "Missing header metadata. Expected lines like:\n"
            "# Path LGXK\n"
            "# Symmetry points: 0 163 351 499"
        )

    label_symbols = list(info0["path"])
    labels_tex = [r"$\Gamma$" if s == "G" else s for s in label_symbols]
    material = info0.get("material", "material")
    sym_indices = info0["symmetry_points"]

    # Load and prep all band files
    bands_list: List[np.ndarray] = []
    arr = arr0
    nb = min(args.nb_bands, arr.shape[0])
    arr = arr[:nb, :]
    if args.fermi != 0.0:
        arr = arr - args.fermi
    bands_list.append(arr)

    for f in files[1:]:
        arr_i, _ = load_band_data(f)
        nb = min(args.nb_bands, arr_i.shape[0])
        arr_i = arr_i[:nb, :]
        if args.fermi != 0.0:
            arr_i = arr_i - args.fermi
        bands_list.append(arr_i)

    nk, _ = plot_and_save_figs(
        bands_list, labels_tex, material, out_dir, args.nb_bands,
        args.show_plot, tuple(args.ylim) if args.ylim else None, args.zoom_pad,
        sym_indices=sym_indices,
    )

    # Metrics from the first file
    bands = bands_list[0]
    comps = estimate_gap_components(bands)
    sym_json, sym_rows_long, sym_rows_wide = dump_all_bands_at_symmetry_points(
        bands, labels_tex, sym_indices
    )
    bw = bandwidths(bands)

    report: Dict[str, Any] = {
        "material": material,
        "nk_points": nk,
        "nb_bands_analyzed": int(bands.shape[0]),
        "path_symbols": label_symbols,
        "fundamental_gap_eV": comps["fundamental_gap_eV"],
        "gap_is_direct": comps["is_direct"],
        "VBM": comps["VBM"],
        "CBM": comps["CBM"],
        "direct_gap_min_eV": comps["direct_gap_min_eV"],
        "direct_gap_min_k_index": comps["direct_gap_min_k_index"],
        "bandwidths_eV": bw,
        "energies_at_symmetry_points": sym_json,
        "symmetry": {"labels": label_symbols, "indices": sym_indices},
        "input_header": {
            "nb_bands": info0.get("nb_bands"),
            "nonlocal": info0.get("nonlocal"),
            "header_lines": info0.get("header_lines"),
        },
    }

    base = args.report or f"band_metrics_{material}"
    json_path = out_dir / f"{base}.json"
    csv_path = out_dir / f"{base}.csv"
    csv_sym_long = out_dir / f"{base}_symmetry_bands_long.csv"
    csv_sym_wide = out_dir / f"{base}_symmetry_bands_wide.csv"

    print(json.dumps(report, indent=2))
    with open(json_path, "w", encoding="utf-8") as jf:
        json.dump(report, jf, indent=2)
    with open(csv_path, "w", newline="", encoding="utf-8") as cf:
        writer = csv.writer(cf)
        writer.writerow(["k_index", "valence_max_eV", "conduction_min_eV", "direct_gap_Eg_k_eV"])
        valence = np.array(comps["valence_max_per_k_eV"])
        conduction = np.array(comps["conduction_min_per_k_eV"])
        Egk = conduction - valence
        for k in range(nk):
            writer.writerow([k, valence[k], conduction[k], Egk[k]])
    with open(csv_sym_long, "w", newline="", encoding="utf-8") as cf:
        writer = csv.writer(cf)
        writer.writerow(["symmetry_label", "k_index", "band_index", "energy_eV"])
        for lab, kidx, bi, e in sym_rows_long:
            writer.writerow([lab, kidx, bi, e])
    with open(csv_sym_wide, "w", newline="", encoding="utf-8") as cf:
        writer = csv.writer(cf)
        header = ["symmetry_label", "k_index"] + [f"band_{i}" for i in range(bands.shape[0])]
        writer = csv.writer(cf)
        writer.writerow(header)
        for lab, kidx, energies in sym_rows_wide:
            writer.writerow([lab, kidx] + energies)

    print(f"Saved: {json_path}")
    print(f"Saved: {csv_path}")
    print(f"Saved: {csv_sym_long}")
    print(f"Saved: {csv_sym_wide}")


if __name__ == "__main__":
    main()
