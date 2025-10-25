#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Band structure plotting + metrics (enhanced) + full symmetry-point band dump
----------------------------------------------------------------------------
Adds extraction of standard band-structure metrics (JSON/CSV), *plus* a complete
dump of ALL band energies at each symmetry point (both JSON and a tidy CSV).

New outputs:
- JSON: report["bands_at_symmetry_points"][<label>] contains:
    - "k_index": nearest sample index to the symmetry point
    - "bands_eV": list of energies for all analyzed bands at that k
    - "valence_max_eV", "conduction_min_eV", "direct_gap_eV" (redundant but handy)
- CSV (tidy/long): "<basename>_symmetry_bands_long.csv" with rows:
    symmetry_label, k_index, band_index, energy_eV
- CSV (wide): "<basename>_symmetry_bands_wide.csv" with rows per symmetry:
    symmetry_label, k_index, band_0, band_1, ..., band_{NB-1}

Notes:
- Energies are assumed referenced to 0 eV = VBM/Fermi (use --fermi to shift if needed).
- Effective masses still omitted (need lattice in Å).
"""

from __future__ import annotations
import argparse
import json
import csv
from pathlib import Path
from typing import Dict, List, Sequence, Tuple, Any

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

try:
    plt.style.use(['science', 'high-vis', 'no-latex'])
except Exception:
    pass

mpl.rcParams['figure.figsize'] = [3.8, 3.0]
mpl.rcParams['axes.grid'] = True
mpl.rcParams['grid.alpha'] = 0.35
mpl.rcParams['grid.linewidth'] = 0.4

BZ_POINTS = {
    "G":  np.array([0.0, 0.0, 0.0]),
    "X":  np.array([0.0, 0.5, 0.0]),
    "L":  np.array([0.25, 0.25, 0.25]),
    "W":  np.array([0.25, 0.5, 0.0]),
    "U":  np.array([0.125, 0.5, 0.125]),
    "K":  np.array([0.375, 0.375, 0.0]),
}
GAMMA_TEX = r"$\Gamma$"


def parse_path_from_filename(filename: Path) -> Tuple[List[str], List[float], str, List[str]]:
    stem = filename.stem
    parts = stem.split("_")
    try:
        material = parts[parts.index("EPM") + 1]
    except ValueError:
        material = "material"
    try:
        path_str = parts[parts.index("path") + 1]
    except ValueError:
        raise ValueError(f"Could not find 'path' in '{stem}'. Expected ..._path_GXUKL_...")
    labels = list(path_str)
    labels_tex = [GAMMA_TEX if s == "G" else s for s in labels]
    pts = [BZ_POINTS[sym] for sym in labels]
    dists = []
    for i in range(len(pts) - 1):
        if labels[i] == "U" and labels[i + 1] == "K":
            dists.append(0.0)
        else:
            dists.append(float(np.linalg.norm(pts[i + 1] - pts[i])))
    return labels_tex, dists, material, labels


def cumulative_from_dists(dists: Sequence[float]) -> np.ndarray:
    cum = np.zeros(len(dists) + 1, dtype=float)
    for i, d in enumerate(dists, start=1):
        cum[i] = cum[i - 1] + d
    return cum


def load_band_data(filename: Path) -> np.ndarray:
    # Keep user's chosen skiprows=5 convention
    data = np.loadtxt(filename, delimiter=",", comments="#", skiprows=5)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data.T  # (nbands, nk)


def estimate_gap_components(bands: np.ndarray, eps: float = 1e-1):
    nb, nk = bands.shape
    valence_max = np.full(nk, -np.inf, dtype=float)
    valence_band_idx = np.full(nk, -1, dtype=int)
    conduction_min = np.full(nk, np.inf, dtype=float)
    conduction_band_idx = np.full(nk, -1, dtype=int)

    for b in range(nb):
        E = bands[b]
        mask_v = E <= eps
        if np.any(mask_v):
            better = (E > valence_max) & mask_v
            valence_max[better] = E[better]
            valence_band_idx[better] = b
        mask_c = E >= eps
        if np.any(mask_c):
            better = (E < conduction_min) & mask_c
            conduction_min[better] = E[better]
            conduction_band_idx[better] = b

    vbm_k = int(np.nanargmax(valence_max))
    cbm_k = int(np.nanargmin(conduction_min))
    vbm_E = float(valence_max[vbm_k])
    cbm_E = float(conduction_min[cbm_k])
    vbm_band = int(valence_band_idx[vbm_k])
    cbm_band = int(conduction_band_idx[cbm_k])
    fundamental_gap = cbm_E - vbm_E
    is_direct = (vbm_k == cbm_k)

    Eg_k = conduction_min - valence_max
    min_direct_k = int(np.nanargmin(Eg_k))
    min_direct_gap = float(Eg_k[min_direct_k])

    return {
        "VBM": {"E_eV": vbm_E, "band": vbm_band, "k_index": vbm_k},
        "CBM": {"E_eV": cbm_E, "band": cbm_band, "k_index": cbm_k},
        "fundamental_gap_eV": float(fundamental_gap),
        "is_direct": bool(is_direct),
        "direct_gap_min_eV": float(min_direct_gap),
        "direct_gap_min_k_index": min_direct_k,
        "valence_max_per_k_eV": valence_max.tolist(),
        "conduction_min_per_k_eV": conduction_min.tolist(),
    }


def nearest_indices_for_symmetry_ticks(nk: int, dists: Sequence[float]) -> np.ndarray:
    cum = cumulative_from_dists(dists)
    if cum[-1] == 0.0:
        ticks = np.linspace(0, nk - 1, len(cum))
    else:
        ticks = (cum / cum[-1]) * (nk - 1)
    return np.rint(ticks).astype(int)


def dump_all_bands_at_symmetry_points(
    bands: np.ndarray,
    labels_tex: Sequence[str],
    dists: Sequence[float],
) -> Tuple[Dict[str, Any], List[Tuple[str, int, int, float]], List[Tuple[str, int, List[float]]]]:
    nk = bands.shape[1]
    tick_idx = nearest_indices_for_symmetry_ticks(nk, dists)
    dict_json: Dict[str, Any] = {}
    rows_long: List[Tuple[str, int, int, float]] = []
    rows_wide: List[Tuple[str, int, List[float]]] = []

    for lab, kidx in zip(labels_tex, tick_idx):
        col = bands[:, kidx]  # (nb,)
        vmax = float(col[col <= 0.0].max()) if np.any(col <= 0.0) else float('nan')
        cmin = float(col[col >= 0.0].min()) if np.any(col >= 0.0) else float('nan')
        dgap = (cmin - vmax) if np.isfinite(vmax) and np.isfinite(cmin) else float('nan')
        dict_json[lab] = {
            "k_index": int(kidx),
            "bands_eV": col.astype(float).tolist(),
            "valence_max_eV": vmax,
            "conduction_min_eV": cmin,
            "direct_gap_eV": float(dgap),
        }
        for bi, e in enumerate(col):
            rows_long.append((lab, int(kidx), int(bi), float(e)))
        rows_wide.append((lab, int(kidx), col.astype(float).tolist()))
    return dict_json, rows_long, rows_wide


def plot_and_save_figs(
    bands_list: List[np.ndarray],
    labels_tex,
    dists,
    material,
    out_dir: Path,
    nbbands: int,
    show_plot: bool,
    ylim: Tuple[float, float] | None,
    zoom_pad: float,
) -> Tuple[int, float]:
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

    cum = cumulative_from_dists(dists)
    ticks = (cum / cum[-1]) * (nk_ref - 1) if cum[-1] != 0 else np.linspace(0, nk_ref - 1, len(cum))
    ax.set_xticks(ticks, labels=labels_tex)
    ax.set_xticks([], minor=True)

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


def bandwidths(bands: np.ndarray) -> List[float]:
    return (bands.max(axis=1) - bands.min(axis=1)).astype(float).tolist()


def main(argv: Sequence[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Plot bands and extract metrics + full symmetry-point bands.")
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
    ap.add_argument("--no-plot", dest="show_plot", action="store_false", help="Do not display the plot window.")
    ap.set_defaults(show_plot=True)
    ap.add_argument("--report", default=None,
                    help="Basename for metrics reports (without extension). Defaults to 'band_metrics_<material>'.")
    args = ap.parse_args(argv)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    files = [Path(f) for f in args.band_files]
    labels_tex, dists, material, label_symbols = parse_path_from_filename(files[0])

    bands_list = []
    for f in files:
        arr = load_band_data(f)
        nb = min(args.nb_bands, arr.shape[0])
        arr = arr[:nb, :]
        if args.fermi != 0.0:
            arr = arr - args.fermi
        bands_list.append(arr)

    nk, _ = plot_and_save_figs(bands_list, labels_tex, dists, material, out_dir, args.nb_bands, args.show_plot, tuple(args.ylim) if args.ylim else None, args.zoom_pad)

    # Metrics from the FIRST file
    bands = bands_list[0]
    comps = estimate_gap_components(bands)
    sym_json, sym_rows_long, sym_rows_wide = dump_all_bands_at_symmetry_points(bands, labels_tex, dists)
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
    }

    base = args.report or f"band_metrics_{material}"
    json_path = out_dir / f"{base}.json"
    csv_path = out_dir / f"{base}.csv"
    csv_sym_long = out_dir / f"{base}_symmetry_bands_long.csv"
    csv_sym_wide = out_dir / f"{base}_symmetry_bands_wide.csv"

    print("Band structure metrics report:")
    print(json.dumps(report, indent=2))
    with open(json_path, "w", encoding="utf-8") as jf:
        json.dump(report, jf, indent=2)
    print(f"Saved metrics JSON: {json_path}")

    valence = np.array(comps["valence_max_per_k_eV"])
    conduction = np.array(comps["conduction_min_per_k_eV"])
    Egk = conduction - valence
    with open(csv_path, "w", newline="", encoding="utf-8") as cf:
        writer = csv.writer(cf)
        writer.writerow(["k_index", "valence_max_eV", "conduction_min_eV", "direct_gap_Eg_k_eV"])
        for k in range(nk):
            writer.writerow([k, valence[k], conduction[k], Egk[k]])
    print(f"Saved metrics CSV: {csv_path}")

    with open(csv_sym_long, "w", newline="", encoding="utf-8") as cf:
        writer = csv.writer(cf)
        writer.writerow(["symmetry_label", "k_index", "band_index", "energy_eV"])
        for lab, kidx, bi, e in sym_rows_long:
            writer.writerow([lab, kidx, bi, e])
    print(f"Saved symmetry ALL-bands (long) CSV: {csv_sym_long}")

    with open(csv_sym_wide, "w", newline="", encoding="utf-8") as cf:
        writer = csv.writer(cf)
        header = ["symmetry_label", "k_index"] + [f"band_{i}" for i in range(bands.shape[0])]
        writer.writerow(header)
        for lab, kidx, energies in sym_rows_wide:
            writer.writerow([lab, kidx] + energies)
    print(f"Saved symmetry ALL-bands (wide) CSV: {csv_sym_wide}")


if __name__ == "__main__":
    print("Enhanced Band Structure Plotter + Metrics Extractor (with full symmetry-point bands)")
    main()