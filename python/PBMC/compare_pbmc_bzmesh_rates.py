#!/usr/bin/env python3
"""Compare analytical PBMC phonon rates with BZ_MESH/FBMC phonon rates.

The script can generate PBMC rates by calling build/apps/pbmc_rates, then
optionally merge them with an existing BZ_MESH rates-vs-energy CSV produced by
elph.epm. It writes a compact comparison CSV and, when matplotlib is available,
a log-scale PNG plot.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import subprocess
from pathlib import Path


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pbmc-rates", type=Path, help="Existing PBMC rates CSV.")
    parser.add_argument("--bzmesh-rates", type=Path, help="BZ_MESH rates-vs-energy CSV from elph.epm.")
    parser.add_argument("--out-prefix", type=Path, default=root / "rate_model_comparison")
    parser.add_argument("--build-dir", type=Path, default=root / "build")
    parser.add_argument("--material", default="Si")
    parser.add_argument("--pbmc-set", default="default")
    parser.add_argument("--part-type", choices=("electron", "hole"), default="electron")
    parser.add_argument("--valley", type=int, default=0, help="Electron valley or hole band index for PBMC.")
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--max-energy", type=float, default=1.0)
    parser.add_argument("--energy-step", type=float, default=0.005)
    parser.add_argument("--no-plot", action="store_true")
    return parser.parse_args()


def read_table(filename: Path) -> list[dict[str, float]]:
    with filename.open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows: list[dict[str, float]] = []
        for row in reader:
            parsed = {}
            for key, value in row.items():
                if key is not None and value not in (None, ""):
                    parsed[key.strip()] = float(value)
            rows.append(parsed)
    if not rows:
        raise RuntimeError(f"{filename} has no data rows")
    return rows


def run_pbmc_export(args: argparse.Namespace, output_file: Path) -> None:
    exe = args.build_dir / "apps" / "pbmc_rates"
    if not exe.exists():
        raise FileNotFoundError(f"{exe} does not exist. Build target pbmc_rates first.")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(exe),
        "--material",
        args.material,
        "--pbmc-set",
        args.pbmc_set,
        "--part-type",
        args.part_type,
        "--valley",
        str(args.valley),
        "--temperature",
        str(args.temperature),
        "--max-energy",
        str(args.max_energy),
        "--energy-step",
        str(args.energy_step),
        "--out",
        str(output_file),
    ]
    print("running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def pbmc_summary(row: dict[str, float]) -> dict[str, float]:
    inelastic = (
        row.get("intervalley_absorption_s_1", 0.0)
        + row.get("intervalley_emission_s_1", 0.0)
        + row.get("optical_absorption_s_1", 0.0)
        + row.get("optical_emission_s_1", 0.0)
    )
    absorption = row.get("intervalley_absorption_s_1", 0.0) + row.get("optical_absorption_s_1", 0.0)
    emission = row.get("intervalley_emission_s_1", 0.0) + row.get("optical_emission_s_1", 0.0)
    return {
        "energy_eV": row["energy_eV"],
        "pbmc_acoustic_s_1": row.get("acoustic_s_1", 0.0),
        "pbmc_inelastic_s_1": inelastic,
        "pbmc_absorption_s_1": absorption,
        "pbmc_emission_s_1": emission,
        "pbmc_total_s_1": row.get("total_s_1", row.get("acoustic_s_1", 0.0) + inelastic),
    }


def bzmesh_summary(row: dict[str, float]) -> dict[str, float]:
    rate_items = [(key, value) for key, value in row.items() if key.startswith("rate_")]
    acoustic = sum(value for key, value in rate_items if key.startswith("rate_ac_"))
    optical = sum(value for key, value in rate_items if key.startswith("rate_op_"))
    absorption = sum(value for key, value in rate_items if key.endswith("_ab"))
    emission = sum(value for key, value in rate_items if key.endswith("_em"))
    total = sum(value for _, value in rate_items)
    return {
        "energy_eV": row.get("energy_eV", row.get("energy(eV)")),
        "bzmesh_acoustic_s_1": acoustic,
        "bzmesh_inelastic_s_1": optical,
        "bzmesh_absorption_s_1": absorption,
        "bzmesh_emission_s_1": emission,
        "bzmesh_total_s_1": total,
    }


def interpolate(x: list[float], y: list[float], xq: float) -> float | None:
    if xq < x[0] or xq > x[-1]:
        return None
    i = bisect.bisect_left(x, xq)
    if i < len(x) and x[i] == xq:
        return y[i]
    if i == 0 or i == len(x):
        return None
    x0, x1 = x[i - 1], x[i]
    y0, y1 = y[i - 1], y[i]
    if x1 == x0:
        return y0
    t = (xq - x0) / (x1 - x0)
    return y0 + t * (y1 - y0)


def merge_tables(pbmc_rows: list[dict[str, float]], bz_rows: list[dict[str, float]] | None) -> list[dict[str, float]]:
    pbmc = [pbmc_summary(row) for row in pbmc_rows]
    if bz_rows is None:
        return pbmc

    bz = [bzmesh_summary(row) for row in bz_rows]
    bz = [row for row in bz if row["energy_eV"] is not None]
    bz.sort(key=lambda row: row["energy_eV"])
    x = [row["energy_eV"] for row in bz]

    merged = []
    bz_keys = [key for key in bz[0] if key != "energy_eV"]
    columns = {key: [row[key] for row in bz] for key in bz_keys}
    for row in pbmc:
        out = dict(row)
        for key in bz_keys:
            value = interpolate(x, columns[key], row["energy_eV"])
            out[key] = value if value is not None else float("nan")
        if out["bzmesh_total_s_1"] == out["bzmesh_total_s_1"] and out["bzmesh_total_s_1"] != 0.0:
            out["total_ratio_pbmc_over_bzmesh"] = out["pbmc_total_s_1"] / out["bzmesh_total_s_1"]
        else:
            out["total_ratio_pbmc_over_bzmesh"] = float("nan")
        merged.append(out)
    return merged


def write_csv(rows: list[dict[str, float]], filename: Path) -> None:
    filename.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with filename.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: f"{row.get(key, float('nan')):.12g}" for key in keys})


def write_plot(rows: list[dict[str, float]], filename: Path) -> bool:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    energy = [row["energy_eV"] for row in rows]
    plt.figure(figsize=(8.0, 5.0))
    plt.semilogy(energy, [row["pbmc_total_s_1"] for row in rows], label="PBMC total")
    plt.semilogy(energy, [row["pbmc_acoustic_s_1"] for row in rows], label="PBMC acoustic")
    plt.semilogy(energy, [row["pbmc_inelastic_s_1"] for row in rows], label="PBMC inelastic")

    if "bzmesh_total_s_1" in rows[0]:
        plt.semilogy(energy, [row["bzmesh_total_s_1"] for row in rows], "--", label="BZ_MESH total")
        plt.semilogy(energy, [row["bzmesh_acoustic_s_1"] for row in rows], "--", label="BZ_MESH acoustic")
        plt.semilogy(energy, [row["bzmesh_inelastic_s_1"] for row in rows], "--", label="BZ_MESH optical")

    plt.xlabel("Carrier kinetic energy (eV)")
    plt.ylabel("Scattering rate (s$^{-1}$)")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend()
    plt.tight_layout()
    filename.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(filename, dpi=180)
    plt.close()
    return True


def main() -> None:
    args = parse_args()
    pbmc_file = args.pbmc_rates or args.out_prefix.with_name(args.out_prefix.name + "_pbmc_rates.csv")
    if args.pbmc_rates is None:
        run_pbmc_export(args, pbmc_file)

    pbmc_rows = read_table(pbmc_file)
    bz_rows = read_table(args.bzmesh_rates) if args.bzmesh_rates else None

    merged = merge_tables(pbmc_rows, bz_rows)
    comparison_file = args.out_prefix.with_name(args.out_prefix.name + "_comparison.csv")
    write_csv(merged, comparison_file)
    print(f"wrote {comparison_file}")

    if not args.no_plot:
        plot_file = args.out_prefix.with_name(args.out_prefix.name + "_comparison.png")
        if write_plot(merged, plot_file):
            print(f"wrote {plot_file}")
        else:
            print("matplotlib is not available; skipped plot")


if __name__ == "__main__":
    main()
