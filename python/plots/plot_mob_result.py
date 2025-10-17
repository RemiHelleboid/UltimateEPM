#!/usr/bin/env python3
# plot_mobility_folder.py
import argparse, re, csv
from pathlib import Path
from typing import Optional, Tuple, List
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science', 'grid', 'no-latex'])

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
        # if len(mu_tensor) == 3:
            # out["mu_tensor_cm2"] = np.array(mu_tensor, dtype=float)
    except StopIteration:
        pass

    # Isotropic: after '# Isotropic mobility in cm^2/(V·s)' take next numeric
    try:
        idx_iso = next(i for i,s in enumerate(lines)
                       if s.strip().lower().startswith("# isotropic mobility"))
        for j in range(idx_iso+1, len(lines)):
            s = lines[j].strip()
            if s and not s.startswith("#"):
                print(s)
                m = re.search(NUM, s)
                print(m)
                if m:
                    print(m.group(0))
                    out["mu_iso_cm2"] = float(m.group(0))
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

    keys = rows[0].keys()
    print(f"Parsed keys: {keys}")
    data_keys = [k for k in keys if k not in ("label", "file", "mu_tensor_cm2")]
    print(f"Data keys: {data_keys}")
    data = np.array([[r[k] if r[k] is not None else np.nan for k in data_keys] for r in rows])
    print(f"Data shape: {data.shape}")
    print(f"Data (first 5 rows):\n{data[:5]}")
    data_T = data.T
    print(f"Data transposed shape: {data_T}")
    key_plot = ""
    for idx, k in enumerate(data_keys):
        row = data_T[idx]
        nb_unique = len(np.unique(np.array(row)))
        print(f"Key '{k}' has {nb_unique} unique values")
        if nb_unique > 1 and k != "mu_iso_cm2":
            key_plot = k
            break
    if key_plot:
        print(f"Will plot μ_iso vs {key_plot}")
    else:
        print("No varying key found to plot against μ_iso")


    print(f"Parsed ", rows)
    # write CSV
    csv_path = folder / args.csv
    with csv_path.open("w", newline="") as g:
        w = csv.writer(g)
        w.writerow([
            "label","file","x_guess",
            "mu_iso_cm2Vs",
            "T_K","Ef_eV","nb_vtx","nb_cb","nb_vb","material","band_gap_eV","energy_range_eV"
        ])
        for r in rows:
            x_guess = guess_numeric_x_from_label(r["label"])

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
                "" if r["energy_range_eV"] is None else r["energy_range_eV"]
            ])
    print(f"Wrote {csv_path}")
    if not key_plot:
        print("No plot generated.")
        return
    # plot
    fig, ax = plt.subplots(figsize=(6,4))
    x = []
    y = []
    for r in rows:
        if r["mu_iso_cm2"] is not None and r[key_plot] is not None:
            x_val = r[key_plot]
            if key_plot == "T_K":
                x_val = float(x_val)
            x.append(x_val)
            y.append(r["mu_iso_cm2"])
    x = np.array(x)
    y = np.array(y)
    if (key_plot == "nb_vtx"):
        # Cubic  root scale for number of k-points
        x = 1.0 / np.cbrt(x)
    key_plot = "1 / nb_vtx"
    print(f"Plotting {len(x)} points")
    # sort by x
    idx_sort = np.argsort(x)
    x = x[idx_sort]
    y = y[idx_sort]
    ax.plot(x, y, marker='o', linestyle='-')
    ax.set_xlabel(key_plot.replace("_", " "))
    ax.set_ylabel(r"Isotropic mobility $\mu_{iso}$ (cm$^2$/V·s)")
    ax.set_title(f"Isotropic mobility vs {key_plot}")
    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    print(f"Wrote plot to {args.out}")
    plt.show()

if __name__ == "__main__":
    main()
