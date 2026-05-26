#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import subprocess
import time
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# Style matplotlib with a clean, modern look.
plt.style.use("seaborn-v0_8-whitegrid")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a bulk AMC electric-field sweep and plot velocity/mobility."
    )

    parser.add_argument(
        "--exe",
        required=True,
        type=Path,
        help="Path to the bulk AMC executable.",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Directory where sweep outputs will be written.",
    )

    parser.add_argument(
        "--material",
        default="Si",
        help="Material symbol passed to the simulator.",
    )

    parser.add_argument(
        "--particle",
        default="electron",
        choices=["electron", "hole"],
        help="Carrier type.",
    )

    parser.add_argument(
        "--runner",
        default="self-scattering",
        choices=["self-scattering", "fixed-step"],
        help="Simulation runner.",
    )

    parser.add_argument(
        "--fields",
        type=float,
        nargs="+",
        required=True,
        help="Electric fields in V/cm.",
    )

    parser.add_argument(
        "--npart",
        type=int,
        default=10000,
        help="Number of particles per simulation.",
    )

    parser.add_argument(
        "--time",
        type=float,
        default=5.0e-12,
        help="Simulation final time in seconds.",
    )

    parser.add_argument(
        "--temperature",
        type=float,
        default=300.0,
        help="Lattice temperature in K.",
    )

    parser.add_argument(
        "--max-energy",
        type=float,
        default=10.0,
        help="Maximum energy in eV used for gamma_max computation.",
    )

    parser.add_argument(
        "--gamma-safety",
        type=float,
        default=1.2,
        help="Safety factor for gamma_max.",
    )

    parser.add_argument(
        "--gamma-samples",
        type=int,
        default=1000,
        help="Number of energy samples for gamma_max.",
    )

    parser.add_argument(
        "--warmup",
        type=float,
        default=0.2,
        help="Warmup fraction ignored in steady-state averages.",
    )

    parser.add_argument(
        "--dt",
        type=float,
        default=5.0e-15,
        help="Fixed-step time step in seconds.",
    )

    return parser.parse_args()


def run_one_field(args: argparse.Namespace, field_v_per_cm: float) -> Path:
    run_dir = args.outdir / f"Ex_{field_v_per_cm:.6e}_Vcm"
    run_dir.mkdir(parents=True, exist_ok=True)

    command = [
        str(args.exe),
        "--material",
        args.material,
        "--part-type",
        args.particle,
        "--runner",
        args.runner,
        "--npart",
        str(args.npart),
        "--time",
        str(args.time),
        "--temperature",
        str(args.temperature),
        "--Ex",
        str(field_v_per_cm),
        "--max-energy",
        str(args.max_energy),
        "--gamma-safety",
        str(args.gamma_safety),
        "--gamma-samples",
        str(args.gamma_samples),
        "--warmup",
        str(args.warmup),
        "--dt",
        str(args.dt),
        "--outdir",
        str(run_dir),
    ]

    log_file = run_dir / "stdout.log"

    with log_file.open("w", encoding="utf-8") as stream:
        completed = subprocess.run(
            command,
            stdout=stream,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )

    if completed.returncode != 0:
        raise RuntimeError(
            f"Simulation failed for Ex={field_v_per_cm:.6e} V/cm. "
            f"See log: {log_file}"
        )

    observables_file = run_dir / "observables.csv"

    if not observables_file.exists():
        raise FileNotFoundError(
            f"Simulation completed but did not create {observables_file}"
        )

    return observables_file


def read_last_observable_row(path: Path) -> dict[str, float]:
    with path.open("r", newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    if not rows:
        raise RuntimeError(f"No observable rows found in {path}")

    row = rows[-1]

    return {
        "electric_field_V_per_m": float(row["electric_field_V_per_m"]),
        "mean_velocity_x_m_per_s": float(row["mean_velocity_x_m_per_s"]),
        "mean_kinetic_energy_eV": float(row["mean_kinetic_energy_eV"]),
        "sample_count": float(row["sample_count"]),
    }


def build_sweep_dataframe(args: argparse.Namespace) -> pd.DataFrame:
    records = []

    for field_v_per_cm in args.fields:
        print(f"Running Ex = {field_v_per_cm:.6e} V/cm")

        started = time.perf_counter()
        observables_file = run_one_field(args, field_v_per_cm)
        elapsed = time.perf_counter() - started

        row = read_last_observable_row(observables_file)

        electric_field_v_per_m = row["electric_field_V_per_m"]
        velocity_x_m_per_s = row["mean_velocity_x_m_per_s"]

        mobility_m2_per_v_s = (
            abs(velocity_x_m_per_s) / abs(electric_field_v_per_m)
            if electric_field_v_per_m != 0.0
            else float("nan")
        )

        records.append(
            {
                "field_V_per_cm": field_v_per_cm,
                "field_V_per_m": electric_field_v_per_m,
                "mean_velocity_x_m_per_s": velocity_x_m_per_s,
                "mean_speed_x_abs_m_per_s": abs(velocity_x_m_per_s),
                "mobility_m2_per_V_s": mobility_m2_per_v_s,
                "mobility_cm2_per_V_s": mobility_m2_per_v_s * 1.0e4,
                "mean_kinetic_energy_eV": row["mean_kinetic_energy_eV"],
                "sample_count": row["sample_count"],
                "runtime_s": elapsed,
                "observables_file": str(observables_file),
            }
        )

    return pd.DataFrame.from_records(records).sort_values("field_V_per_cm")


def plot_velocity(df: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots()

    ax.plot(
        df["field_V_per_cm"],
        df["mean_velocity_x_m_per_s"],
        marker="o",
    )

    ax.set_xscale("log")
    ax.set_xlabel("Electric field (V/cm)")
    ax.set_ylabel("Mean drift velocity x (m/s)")
    ax.set_title("Bulk AMC drift velocity versus electric field")
    ax.grid(True, which="both")

    fig.tight_layout()
    fig.savefig(outdir / "velocity_vs_field.png", dpi=200)
    fig.savefig(outdir / "velocity_vs_field.pdf")
    plt.close(fig)


def plot_mobility(df: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots()

    ax.plot(
        df["field_V_per_cm"],
        df["mobility_cm2_per_V_s"],
        marker="o",
    )

    ax.set_xscale("log")
    ax.set_xlabel("Electric field (V/cm)")
    ax.set_ylabel("Mobility (cm²/V/s)")
    ax.set_title("Bulk AMC mobility versus electric field")
    ax.grid(True, which="both")

    fig.tight_layout()
    fig.savefig(outdir / "mobility_vs_field.png", dpi=200)
    fig.savefig(outdir / "mobility_vs_field.pdf")
    plt.close(fig)


def plot_energy(df: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots()

    ax.plot(
        df["field_V_per_cm"],
        df["mean_kinetic_energy_eV"],
        marker="o",
    )

    ax.set_xscale("log")
    ax.set_xlabel("Electric field (V/cm)")
    ax.set_ylabel("Mean kinetic energy (eV)")
    ax.set_title("Bulk AMC mean energy versus electric field")
    ax.grid(True, which="both")

    fig.tight_layout()
    fig.savefig(outdir / "energy_vs_field.png", dpi=200)
    fig.savefig(outdir / "energy_vs_field.pdf")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = build_sweep_dataframe(args)

    results_csv = args.outdir / "field_sweep_results.csv"
    df.to_csv(results_csv, index=False)

    plot_velocity(df, args.outdir)
    plot_mobility(df, args.outdir)
    plot_energy(df, args.outdir)

    print(f"Wrote {results_csv}")
    print(f"Wrote {args.outdir / 'velocity_vs_field.png'}")
    print(f"Wrote {args.outdir / 'mobility_vs_field.png'}")
    print(f"Wrote {args.outdir / 'energy_vs_field.png'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())