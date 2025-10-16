#!/usr/bin/env python3
# plot_mobility_folder.py
import argparse, re, csv
from pathlib import Path
from typing import Optional, Tuple, List
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science', 'no-latex', 'grid'])

NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"

def parse_mobility_file(path: Path):
    """
    Parse a file written by your C++ snippet.
    Returns a dict with:
      label, file, nb_vtx, nb_cb, nb_vb, material, band_gap_eV, energy_range_eV,
      T_K, Ef_eV, mu_tensor_cm2 (3x3), mu_iso_cm2
    Missing fields are set to None.
    """
    out = {
        "label": path.stem[:-9] if path.stem.endswith("_mobility") else path.stem,
        "file": str(path),
        "nb_vtx": None,
        "nb_cb": None,
        "nb_vb": None,
        "material": None,
        "band_gap_eV": None,
        "energy_range_eV": None,
        "T_K": None,
        "Ef_eV": None,
        "mu_tensor_cm2": None,
        "mu_iso_cm2": None,
    }
    lines = path.read_text(errors="ignore").splitlines()

    # simple helpers
    def grab_after(prefix):
        for s in lines:
            if s.strip().startswith(prefix):
                m = re.search(r":\s*(" + NUM + ")", s)
                if m: return float(m.group(1))
        return None

    out["nb_vtx"] = grab_after("# Number of vertices")
    out["nb_cb"]  = grab_after("# Number of conduction bands")
    out["nb_vb"]  = grab_after("# Number of valence bands")
    out["material"] = next((s.split(":",1)[1].strip()
                            for s in lines if s.strip().startswith("# Material")), None)
    out["band_gap_eV"] = grab_after("# Band gap (after scissor) in eV")
    out["energy_range_eV"] = grab_after("# Energy range in eV")
    out["T_K"] = grab_after("# Temperature in Kelvin")
    out["Ef_eV"] = grab_after("# Fermi level in eV")

    # Tensor block: after '# Mobility tensor in cm^2/(V·s)' take next 3 numeric lines
    mu_tensor = []
    try:
        idx_t = next(i for i,s in enumerate(lines)
                     if s.strip().lower().startswith("# mobility tensor in"))
        i = idx_t + 1
        while i < len(lines) and len(mu_tensor) < 3:
            s = lines[i].strip()
            if s and not s.startswith("#") and re.search(NUM, s):
                row = [float(x) for x in re.findall(NUM, s)]
                if len(row) >= 3:
                    mu_tensor.append(row[:3])
            i += 1
        if len(mu_tensor) == 3:
            out["mu_tensor_cm2"] = np.array(mu_tensor, dtype=float)
    except StopIteration:
        pass

    # Isotropic: after '# Isotropic mobility in cm^2/(V·s)' take next numeric
    try:
        idx_iso = next(i for i,s in enumerate(lines)
                       if s.strip().lower().startswith("# isotropic mobility"))
        for j in range(idx_iso+1, len(lines)):
            s = lines[j].strip()
            if s and not s.startswith("#"):
                m = re.search(NUM, s)
                if m:
                    out["mu_iso_cm2"] = float(m.group(1))
                break
    except StopIteration:
        pass

    return out

def guess_numeric_x_from_label(label: str) -> Optional[float]:
    """
    Try to pull a number from label, accepting 0p008 as 0.008 etc.
    """
    s = re.sub(r'(?<=\d)p(?=\d)', '.', label)
    m = re.search(NUM, s)
    return float(m.group(0)) if m else None

def main():
    ap = argparse.ArgumentParser(description="Read *_mobility.txt files and plot isotropic mobility.")
    ap.add_argument("folder", help="Folder containing *_mobility.txt")
    ap.add_argument("--out", default="mu_plot.png", help="Output plot filename")
    ap.add_argument("--csv", default="summary.csv", help="Summary CSV filename")
    args = ap.parse_args()

    folder = Path(args.folder)
    files = sorted(folder.glob("*_mobility.txt"))
    if not files:
        print(f"No *_mobility.txt found in {folder}")
        return

    rows: List[dict] = [parse_mobility_file(f) for f in files]
    # write CSV
    csv_path = folder / args.csv
    with csv_path.open("w", newline="") as g:
        w = csv.writer(g)
        w.writerow([
            "label","file","x_guess",
            "mu_iso_cm2Vs",
            "T_K","Ef_eV","nb_vtx","nb_cb","nb_vb","material","band_gap_eV","energy_range_eV",
            "mu_tensor_00","mu_tensor_01","mu_tensor_02",
            "mu_tensor_10","mu_tensor_11","mu_tensor_12",
            "mu_tensor_20","mu_tensor_21","mu_tensor_22",
        ])
        for r in rows:
            x_guess = guess_numeric_x_from_label(r["label"])
            mt = r["mu_tensor_cm2"]
            mt_vals = (mt.flatten().tolist() if isinstance(mt, np.ndarray) else [None]*9)
            w.writerow([
                r["label"], r["file"], "" if x_guess is None else x_guess,
                "" if r["mu_iso_cm2"] is None else r["mu_iso_cm2"],
                "" if r["T_K"] is None else r["T_K"],
                "" if r["Ef_eV"] is None else r["Ef_eV"],
                "" if r["nb_vtx"] is None else r["nb_vtx"],
                "" if r["nb_cb"] is None else r["nb_cb"],
                "" if r["nb_vb"] is None else r["nb_vb"],
                "" if r["material"] is None else r["material"],
                "" if r["band_gap_eV"] is None else r["band_gap_eV"],
                "" if r["energy_range_eV"] is None else r["energy_range_eV"],
                *mt_vals
            ])
    print(f"Wrote {csv_path}")

    # Plot μ_iso vs x_guess if numbers exist, else bar by label
    pairs = []
    for r in rows:
        mu = r["mu_iso_cm2"]
        if mu is None: continue
        x = guess_numeric_x_from_label(r["label"])
        pairs.append((x, mu, r["label"]))

    out_path = folder / args.out
    if pairs and all(p[0] is not None for p in pairs):
        pairs.sort(key=lambda t: t[0])
        xs = [p[0] for p in pairs]
        ys = [p[1] for p in pairs]
        plt.figure(figsize=(7,5))
        plt.plot(xs, ys, marker="o")
        plt.xlabel("x (guessed from filename)")
        plt.ylabel("μ_iso (cm²/V·s)")
        plt.title("Isotropic mobility vs x")
        plt.tight_layout()
        plt.savefig(out_path, dpi=150)
        plt.close()
    else:
        # bar chart by label
        labels = [p[2] for p in pairs]
        ys = [p[1] for p in pairs]
        if not ys:
            print("No μ_iso values found to plot.")
            return
        plt.figure(figsize=(max(6, 0.7*len(labels)), 5))
        pos = np.arange(len(labels))
        plt.bar(pos, ys)
        plt.xticks(pos, labels, rotation=45, ha="right")
        plt.ylabel("μ_iso (cm²/V·s)")
        plt.title("Isotropic mobility by file")
        plt.tight_layout()
        plt.savefig(out_path, dpi=150)
        plt.close()

    print(f"Saved plot to {out_path}")

if __name__ == "__main__":
    main()
