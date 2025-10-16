#!/usr/bin/env python3
# plot_mobility_folder.py
import argparse, re, csv
from pathlib import Path
from typing import Optional, Tuple, List
import numpy as np
import matplotlib.pyplot as plt

import scienceplots
plt.style.use(['science', 'notebook', 'grid'])

def parse_mu_iso(path: Path) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """
    Returns: (mu_iso_cm2Vs, temperature_K, Ef_eV) if found; otherwise None for missing.
    """
    try:
        lines = path.read_text(errors="ignore").splitlines()
    except Exception:
        return None, None, None

    T = None
    Ef = None
    for line in lines:
        if line.startswith("# Temperature in Kelvin"):
            m = re.search(r":\s*([+-]?\d+(?:\.\d+)?)", line)
            if m: T = float(m.group(1))
        elif line.startswith("# Fermi level in eV"):
            m = re.search(r":\s*([+-]?\d+(?:\.\d+)?)", line)
            if m: Ef = float(m.group(1))

    # Find the isotropic header, then take the next numeric line
    mu_iso = None
    iso_idx = None
    for i, line in enumerate(lines):
        if line.strip().lower().startswith("# isotropic mobility"):
            iso_idx = i
            break
    if iso_idx is not None:
        for j in range(iso_idx + 1, len(lines)):
            s = lines[j].strip()
            if s and not s.startswith("#"):
                m = re.search(r"([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", s)
                if m:
                    mu_iso = float(m.group(1))  # already cm^2/(V·s) per your file
                break
    return mu_iso, T, Ef

def guess_numeric_x_from_stem(stem: str) -> Optional[float]:
    """
    Try to pull a number from the filename stem.
    Accepts formats like: 0.02, 0p02, h_0p02, mesh0p008, etc.
    Returns float or None.
    """
    # Replace 'p' with '.' in number-like chunks (e.g., 0p02 -> 0.02)
    s = re.sub(r'(?<=\d)p(?=\d)', '.', stem)
    # Grab the first float-looking token
    m = re.search(r"([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", s)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None

def main():
    ap = argparse.ArgumentParser(description="Plot isotropic mobility from *_mobility.txt files in a folder.")
    ap.add_argument("folder", help="Path to folder with *_mobility.txt files")
    ap.add_argument("--out", default="mobility_plot.png", help="Output plot filename (PNG)")
    ap.add_argument("--csv", default="summary.csv", help="Output CSV filename")
    args = ap.parse_args()

    folder = Path(args.folder)
    files = sorted(folder.glob("*_mobility.txt"))
    if not files:
        print(f"No *_mobility.txt files found in {folder}")
        return

    labels: List[str] = []
    xs: List[Optional[float]] = []
    mus: List[Optional[float]] = []
    Ts: List[Optional[float]] = []
    Efs: List[Optional[float]] = []

    for f in files:
        stem = f.stem
        # remove trailing "_mobility" if present for cleaner label
        label = stem[:-9] if stem.endswith("_mobility") else stem
        mu, T, Ef = parse_mu_iso(f)
        x = guess_numeric_x_from_stem(label)
        labels.append(label)
        xs.append(x)
        mus.append(mu)
        Ts.append(T)
        Efs.append(Ef)

    # Write CSV
    csv_path = folder / args.csv
    with csv_path.open("w", newline="") as g:
        w = csv.writer(g)
        w.writerow(["label", "x_guess", "mu_iso_cm2Vs", "T_K", "Ef_eV", "file"])
        for label, x, mu, T, Ef, f in zip(labels, xs, mus, Ts, Efs, files):
            w.writerow([label, "" if x is None else x,
                        "" if mu is None else mu,
                        "" if T is None else T,
                        "" if Ef is None else Ef,
                        str(f)])
    print(f"Wrote {csv_path}")

    # Prepare data for plotting
    # Prefer numeric X if all present; otherwise bar chart by label
    has_all_mu = all(m is not None for m in mus)
    if not has_all_mu:
        print("Warning: some files missing μ_iso; those will be skipped in plot.")

    plot_pairs = [(x, m, lbl) for x, m, lbl in zip(xs, mus, labels) if m is not None]

    if plot_pairs and all(p[0] is not None for p in plot_pairs):
        # Scatter/line vs numeric x
        plot_pairs.sort(key=lambda t: t[0])
        x_vals = [p[0] for p in plot_pairs]
        y_vals = [p[1] for p in plot_pairs]

        plt.figure(figsize=(7,5))
        plt.plot(x_vals, y_vals, marker="o")
        plt.xlabel("x (guessed from filename)")
        plt.ylabel("Isotropic μ (cm²/V·s)")
        plt.title("Isotropic mobility vs filename-derived x")
        plt.tight_layout()
        plt.savefig(folder / args.out, dpi=150)
        plt.show()
        plt.close()
        print(f"Saved plot to {folder / args.out}")
    else:
        # Bar chart by label
        labels_plot = [lbl for (_, m, lbl) in plot_pairs]
        y_vals = [m for (_, m, _) in plot_pairs]
        if not y_vals:
            print("No valid μ_iso values to plot.")
            return
        plt.figure(figsize=(max(6, 0.7*len(labels_plot)), 5))
        pos = np.arange(len(labels_plot))
        plt.bar(pos, y_vals)
        plt.xticks(pos, labels_plot, rotation=45, ha="right")
        plt.ylabel("Isotropic μ (cm²/V·s)")
        plt.title("Isotropic mobility by file")
        plt.tight_layout()
        plt.show()
        plt.savefig(folder / args.out, dpi=150)
        plt.close()
        print(f"Saved plot to {folder / args.out}")

if __name__ == "__main__":
    main()
    plt.show()
