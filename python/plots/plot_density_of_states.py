#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from argparse import ArgumentParser

import matplotlib.pyplot as plt


# ---- optional styles (won't crash if missing) ----
try:
    import matplotlib as mpl
    import scienceplots  # noqa: F401
    plt.style.use(['science', 'muted'])
    mpl.rcParams['figure.figsize'] = [3.5, 2.8]
    mpl.rcParams['figure.dpi'] = 300
except Exception:
    pass

# ---- Silicon defaults ----
# Lattice constant (m)
A_SI = 5.431e-10
# Primitive cell volume for diamond/fcc
VCELL_SI = (A_SI**3) / 4.0
# Spin degeneracy
GS_SI = 2


def load_bands_csv(filename: str):
    arr = np.loadtxt(filename, delimiter=',', skiprows=1)
    if arr.ndim == 1:
        raise ValueError("File has one row only; need at least 2 rows.")
    nb_cols = arr.shape[1]
    if nb_cols % 2 != 0:
        nb_cols -= 1
        arr = arr[:, :nb_cols]

    nbands = nb_cols // 2
    energies, doss = [], []

    for b in range(nbands):
        E = arr[:, 2 * b]
        G = arr[:, 2 * b + 1]
        order = np.argsort(E)
        energies.append(E[order])
        doss.append(G[order])
    return energies, doss


def band_filter_mask(E: np.ndarray, band_type: str) -> bool:
    mu = float(np.nanmean(E))
    if band_type == "conduction":
        return mu >= 0.0
    if band_type == "valence":
        return mu <= 0.0
    return True


def plot_dos_per_band(filename: str, ax=None, band_type: str = "all"):
    if ax is None:
        _, ax = plt.subplots()
    energies, doss = load_bands_csv(filename)
    kept = 0
    for idx, (E, G) in enumerate(zip(energies, doss), 1):
        if not band_filter_mask(E, band_type):
            continue
        kept += 1
        ax.plot(E, G, lw=0.9, label=f"{kept}")
    if kept == 0:
        ax.text(0.5, 0.5, "No bands selected",
                ha="center", va="center", transform=ax.transAxes)

    ax.legend(fontsize='x-small', title_fontsize='x-small',
              title="Band index", fancybox=True, ncol=2)
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("Density of states (a.u.)")
    ax.set_ylim(bottom=0.0)

    if band_type == "conduction":
        ax.set_xlim(left=-1.0)
    elif band_type == "valence":
        ax.set_xlim(right=1.0)
    return ax

def plot_dos_sum_bands(
    filename: str,
    ax: Optional[plt.Axes] = None,
    band_type: str = "all",
    emin: Optional[float] = None,
    emax: Optional[float] = None,
    npts: int = 2000
) -> Tuple[plt.Axes, np.ndarray, np.ndarray]:
    """
    Plot the summed DOS across bands matching `band_type`.

    Returns (ax, eplot, gtot).
    """
    if ax is None:
        _, ax = plt.subplots()

    energies, doss = load_bands_csv(filename)  # expects: list/iter of np.ndarray pairs

    # Collect energies from bands that pass the filter (handle empty gracefully)
    filt_energies = [E for E in energies if band_filter_mask(E, band_type)]
    if filt_energies:
        # nan-safe global mins/maxes
        try:
            allE = np.concatenate(filt_energies)
        except ValueError:
            allE = np.array([], dtype=float)
    else:
        allE = np.array([], dtype=float)

    if emin is None:
        if band_type == "conduction":
            emin = 0.0
        else:
            emin = float(np.nanmin(allE)) if allE.size else -10.0
    if emax is None:
        if band_type == "valence":
            emax = 0.0
        else:
            emax = float(np.nanmax(allE)) if allE.size else 10.0

    # Ensure finite bounds in case allE was NaN-only
    if not np.isfinite(emin):
        emin = -10.0
    if not np.isfinite(emax):
        emax = 10.0
    if emax <= emin:
        # widen degenerate/invalid range a bit
        emax = emin + 1e-6

    eplot = np.linspace(emin, emax, int(npts))
    gtot = np.zeros_like(eplot)

    count = 0
    for E, G in zip(energies, doss):
        if not band_filter_mask(E, band_type):
            continue
        # Sort by energy for safe interpolation if inputs aren’t strictly increasing
        order = np.argsort(E)
        Es, Gs = np.asarray(E)[order], np.asarray(G)[order]
        gtot += np.interp(eplot, Es, Gs, left=0.0, right=0.0)
        count += 1

    ax.plot(eplot, gtot, c="darkblue", label="Sum of {} band{}".format(count, "" if count == 1 else "s"))
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("Density of states (a.u.)")
    ax.set_ylim(bottom=0.0)
    ax.set_xlim(emin, emax)
    ax.legend(frameon=False, fontsize="x-small")
    return ax, eplot, gtot

def quick_sum_rule_check(e_eV: np.ndarray, g_eV_m3: np.ndarray,
                         Vcell_m3: float, g_s: int = 2) -> float:
    integral = np.trapz(g_eV_m3, e_eV)                 # states/m^3
    target = g_s / Vcell_m3                            # states/m^3
    return integral / target if target != 0 else np.nan


def main():
    p = ArgumentParser()
    p.add_argument("-f", "--file", required=True, help="CSV file to parse.")
    p.add_argument("-o", "--outputdir", default="./", help="Directory to save figures.")
    p.add_argument("-t", "--type", dest="band_type", default="all",
                   choices=["all", "conduction", "valence"],
                   help="Which bands to include.")
    p.add_argument("--emin", type=float, default=None,
                   help="Min energy for total DOS plot (eV).")
    p.add_argument("--emax", type=float, default=None,
                   help="Max energy for total DOS plot (eV).")
    p.add_argument("--dpi", type=int, default=600, help="Figure DPI.")
    p.add_argument("--gs", type=int, default=None,
                   help="Spin degeneracy for sum-rule check (default: 2 for Si).")
    p.add_argument("--vcell", type=float, default=None,
                   help="Primitive cell volume (m^3). Default = Si diamond value.")
    args = p.parse_args()

    gs = args.gs if args.gs is not None else GS_SI
    vcell = args.vcell if args.vcell is not None else VCELL_SI

    in_path = Path(args.file)
    out_dir = Path(args.outputdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = in_path.stem

    print(f"PLOTTING DOS FOR FILE {in_path}")
    print(f"[Defaults] Using g_s={gs}, Vcell={vcell:.3e} m^3 (Silicon)")

    # Number of bands
    energies, doss = load_bands_csv(str(in_path))
    nbands = len(energies)

    fig1, ax1 = plt.subplots()
    plot_dos_per_band(str(in_path), ax1, args.band_type)
    fig1.tight_layout()
    fig1.savefig(out_dir / f"DOS_PER_BAND_{stem}.pdf", dpi=args.dpi)
    fig1.savefig(out_dir / f"DOS_PER_BAND_{stem}.png", dpi=args.dpi)

    fig2, ax2 = plt.subplots()
    ax2, eplot, gtot = plot_dos_sum_bands(str(in_path), ax2, args.band_type,
                                          emin=args.emin, emax=args.emax)
    fig2.tight_layout()
    fig2.savefig(out_dir / f"DOS_TOTAL_{stem}.pdf", dpi=args.dpi)
    fig2.savefig(out_dir / f"DOS_TOTAL_{stem}.png", dpi=args.dpi)

    ratio = quick_sum_rule_check(eplot, gtot, vcell, g_s=gs) / nbands
    print(f"[Sum-rule check] ∫g(E)dE = {ratio:.3f} × (nbands *g_s/Vcell)")

    plt.show()


if __name__ == "__main__":
    main()
