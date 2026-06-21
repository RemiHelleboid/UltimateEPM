#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
import time
from pathlib import Path

if "MPLCONFIGDIR" not in os.environ:
    matplotlib_cache = Path(os.environ.get("TMPDIR", "/tmp")) / "uepm_matplotlib"
    matplotlib_cache.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(matplotlib_cache)

if "XDG_CACHE_HOME" not in os.environ:
    xdg_cache = Path(os.environ.get("TMPDIR", "/tmp")) / "uepm_cache"
    xdg_cache.mkdir(parents=True, exist_ok=True)
    os.environ["XDG_CACHE_HOME"] = str(xdg_cache)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a bulk FBMC electric-field sweep and extract low-field mobility."
    )

    parser.add_argument("--exe", required=True, type=Path, help="Path to fbmc.epm.")
    parser.add_argument(
        "--mesh",
        required=True,
        type=Path,
        help="BZ mesh and band-energy .msh file passed to fbmc.epm.",
    )
    parser.add_argument(
        "--phonon-rates",
        default=None,
        type=Path,
        help="CSV file with full-band phonon scattering rates. If omitted, auto-detect in the current directory.",
    )
    parser.add_argument("--outdir", required=True, type=Path, help="Sweep output directory.")
    parser.add_argument("--material", default="Si", help="Material symbol passed to fbmc.epm.")
    parser.add_argument(
        "--bz-domain",
        choices=("full", "octant"),
        default="full",
        help="Stored Brillouin-zone domain used by fbmc.epm.",
    )
    parser.add_argument(
        "--phonon-params",
        default=None,
        help="Electron-phonon parameter set passed to fbmc.epm, e.g. kamakura, michaillat, or fischetti.",
    )
    parser.add_argument(
        "--phonon-params-file",
        type=Path,
        default=None,
        help="Explicit electron-phonon YAML profile passed to fbmc.epm.",
    )
    parser.add_argument(
        "--fields",
        required=True,
        help="Electric fields in V/cm, comma-separated. Example: -1000,-500,500,1000.",
    )
    parser.add_argument("--npart", type=int, default=1000, help="Number of particles per field.")
    parser.add_argument("--time", type=float, default=5.0e-12, help="Simulation final time in seconds.")
    parser.add_argument(
        "--warmup",
        type=float,
        default=0.2,
        help="Warmup fraction ignored in steady-state averages.",
    )
    parser.add_argument("--temperature", type=float, default=300.0, help="Lattice temperature in K.")
    parser.add_argument("--max-energy", type=float, default=10.0, help="Maximum energy in eV.")
    parser.add_argument("--gamma-safety", type=float, default=1.2, help="Self-scattering rate safety factor.")
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Base random seed. If omitted, each FBMC process generates a random seed.",
    )
    parser.add_argument(
        "--enable-impact-ionization",
        action="store_true",
        help="Enable impact-ionization scattering.",
    )
    parser.add_argument(
        "--ncbands",
        type=int,
        default=-1,
        help="Number of conduction bands to read from the mesh. Use -1 for simulator default.",
    )
    parser.add_argument(
        "--nvbands",
        type=int,
        default=-1,
        help="Number of valence bands to read from the mesh. Use -1 for simulator default.",
    )
    parser.add_argument("--nbthreads", type=int, default=1, help="OpenMP threads per FBMC run.")
    parser.add_argument(
        "--mobility-fit-max-field",
        type=float,
        default=2.0e3,
        help="Maximum |Ex| in V/cm used for low-field mobility extraction.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse an existing *_stats.csv in a field directory when available.",
    )
    parser.add_argument("--show", action="store_true", help="Show plots interactively.")

    args = parser.parse_args()
    args.fields = [float(value) for value in args.fields.split(",") if value.strip()]
    return args


def validate_args(args: argparse.Namespace) -> None:
    if not args.exe.is_file():
        raise FileNotFoundError(f"Executable not found: {args.exe}")
    if not os.access(args.exe, os.X_OK):
        raise PermissionError(f"Executable is not executable: {args.exe}")
    if not args.mesh.is_file():
        raise FileNotFoundError(f"Mesh file not found: {args.mesh}")
    if args.phonon_params is not None and args.phonon_params_file is not None:
        raise ValueError("--phonon-params and --phonon-params-file are mutually exclusive.")
    if args.phonon_params_file is not None and not args.phonon_params_file.is_file():
        raise FileNotFoundError(f"Phonon parameter file not found: {args.phonon_params_file}")
    if args.phonon_rates is None:
        args.phonon_rates = detect_phonon_rates_file(Path.cwd())
        print(f"Auto-detected phonon-rate file: {args.phonon_rates}", flush=True)
    elif not args.phonon_rates.is_file():
        raise FileNotFoundError(f"Phonon-rate file not found: {args.phonon_rates}")

    if len(args.fields) < 2:
        raise ValueError("At least two electric-field values are required.")
    if args.npart <= 0:
        raise ValueError("--npart must be positive.")
    if args.time <= 0.0 or not math.isfinite(args.time):
        raise ValueError("--time must be finite and positive.")
    if args.warmup < 0.0 or args.warmup >= 1.0 or not math.isfinite(args.warmup):
        raise ValueError("--warmup must be finite and in [0, 1).")
    if args.temperature <= 0.0 or not math.isfinite(args.temperature):
        raise ValueError("--temperature must be finite and positive.")
    if args.max_energy <= 0.0 or not math.isfinite(args.max_energy):
        raise ValueError("--max-energy must be finite and positive.")
    if args.gamma_safety < 1.0 or not math.isfinite(args.gamma_safety):
        raise ValueError("--gamma-safety must be finite and at least one.")
    if args.seed is not None and args.seed < 0:
        raise ValueError("--seed must be non-negative.")
    if args.nbthreads <= 0:
        raise ValueError("--nbthreads must be positive.")
    if args.mobility_fit_max_field <= 0.0:
        raise ValueError("--mobility-fit-max-field must be positive.")


def field_directory_name(field_v_per_cm: float) -> str:
    return f"Ex_{field_v_per_cm:+.6e}_Vcm".replace("+", "p").replace("-", "m")


def looks_like_phonon_rates_csv(path: Path) -> bool:
    if not path.is_file() or path.suffix != ".csv":
        return False

    try:
        header = path.open("r", encoding="utf-8").readline()
    except OSError:
        return False

    required_tokens = [
        "vertex_index",
        "local_band_index",
        "energy_eV",
        "rate_ac_L_ab",
        "rate_op_T_em",
    ]
    return all(token in header for token in required_tokens)


def detect_phonon_rates_file(directory: Path) -> Path:
    matches = sorted(path for path in directory.iterdir() if looks_like_phonon_rates_csv(path))
    if not matches:
        raise FileNotFoundError(
            f"No phonon-rate CSV was provided and no matching file was found in {directory}. "
            "Pass --phonon-rates or generate phonon_rates.csv with elph.epm --export-rates."
        )
    if len(matches) > 1:
        match_text = "\n  ".join(str(path) for path in matches)
        raise RuntimeError(
            f"Multiple phonon-rate CSV files were found in {directory}. "
            f"Pass --phonon-rates explicitly. Matches:\n  {match_text}"
        )
    return matches[0]


def newest_stats_file(run_dir: Path) -> Path | None:
    observables_file = run_dir / "observables.csv"
    if observables_file.exists():
        return observables_file

    stats_files = sorted(
        run_dir.glob("simulation_results_*_stats.csv"),
        key=lambda path: path.stat().st_mtime,
    )
    return stats_files[-1] if stats_files else None


def run_one_field(args: argparse.Namespace, field_v_per_cm: float) -> Path:
    run_dir = args.outdir / field_directory_name(field_v_per_cm)
    run_dir.mkdir(parents=True, exist_ok=True)

    if args.resume:
        existing_observables = newest_stats_file(run_dir)
        if existing_observables is not None:
            print(f"Reusing Ex = {field_v_per_cm:.6e} V/cm from {existing_observables}", flush=True)
            return existing_observables

    command = [
        str(args.exe),
        "--meshbandfile",
        str(args.mesh),
        "--material",
        args.material,
        "--bz-domain",
        args.bz_domain,
        "--outdir",
        str(run_dir),
        "--npart",
        str(args.npart),
        "--ncbands",
        str(args.ncbands),
        "--nvbands",
        str(args.nvbands),
        "--nthreads",
        str(args.nbthreads),
        "--maxenergy",
        str(args.max_energy),
        "--gamma-safety",
        str(args.gamma_safety),
        "--time",
        str(args.time),
        "--warmup",
        str(args.warmup),
        "--temperature",
        str(args.temperature),
        "--Ex",
        str(field_v_per_cm),
    ]
    if args.phonon_params_file is not None:
        command.extend(["--phonon-params-file", str(args.phonon_params_file.resolve())])
    else:
        command.extend(["--phonon-params", args.phonon_params or "remi-2026"])
    if args.phonon_rates is not None:
        command.extend(["--phononfile", str(args.phonon_rates)])
    if args.seed is not None:
        command.extend(["--seed", str(args.seed)])
    if args.enable_impact_ionization:
        command.append("--enable-impact-ionization")

    log_file = run_dir / "stdout.log"
    print(f"Running Ex = {field_v_per_cm:.6e} V/cm", flush=True)
    print("  " + " ".join(command), flush=True)

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
            f"FBMC failed for Ex={field_v_per_cm:.6e} V/cm. See log: {log_file}"
        )

    stats_file = newest_stats_file(run_dir)
    if stats_file is None:
        raise FileNotFoundError(f"FBMC completed but no observables.csv or *_stats.csv was found in {run_dir}")

    return stats_file


def read_stats_row(path: Path) -> dict[str, float]:
    with path.open("r", newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    if not rows:
        raise RuntimeError(f"No stats rows found in {path}")

    row = rows[-1]

    velocity_column = None
    for candidate in ("mean_velocity_x_m_per_s", "mean_x_velocity_m_per_s", "mean_velocity_norm_m_per_s"):
        if candidate in row:
            velocity_column = candidate
            break

    if velocity_column is None:
        raise KeyError(
            f"Missing mean velocity column in {path}; expected mean_x_velocity_m_per_s "
            "or legacy mean_velocity_norm_m_per_s"
        )

    ionization_column = None
    for candidate in ("impact_ionization_coefficient_cm_1", "mean_ionization_coeff_1_per_m", "mean_ionization_coeff_1_per_s"):
        if candidate in row:
            ionization_column = candidate
            break

    if ionization_column is None:
        raise KeyError(
            f"Missing ionization coefficient column in {path}; expected "
            "mean_ionization_coeff_1_per_m or legacy mean_ionization_coeff_1_per_s"
        )

    energy_column = "mean_kinetic_energy_eV" if "mean_kinetic_energy_eV" in row else "mean_energy_eV"
    if energy_column not in row:
        raise KeyError(f"Missing mean kinetic energy column in {path}")

    ionization_value = float(row[ionization_column])
    if ionization_column == "impact_ionization_coefficient_cm_1":
        ionization_value *= 100.0

    discarded_carriers = int(float(row.get("discarded_carriers_over_max_energy", 0)))
    run_complete = bool(int(float(row.get("run_complete", 1))))
    if discarded_carriers > 0:
        print(
            f"Warning: {path} discarded {discarded_carriers} carrier(s) above --max-energy; "
            "this field point is incomplete.",
            flush=True,
        )

    return {
        "mean_energy_eV": float(row[energy_column]),
        "mean_velocity_x_m_per_s": float(row[velocity_column]),
        "mean_ionization_coeff_1_per_m": ionization_value,
        "discarded_carriers_over_max_energy": discarded_carriers,
        "run_complete": run_complete,
    }


def build_sweep_dataframe(args: argparse.Namespace) -> pd.DataFrame:
    records: list[dict[str, float | str]] = []

    for field_v_per_cm in args.fields:
        started = time.perf_counter()
        stats_file = run_one_field(args, field_v_per_cm)
        elapsed = time.perf_counter() - started
        stats = read_stats_row(stats_file)

        field_v_per_m = field_v_per_cm * 100.0
        velocity_x_m_per_s = stats["mean_velocity_x_m_per_s"]

        if field_v_per_m != 0.0:
            signed_mobility_m2_per_v_s = velocity_x_m_per_s / field_v_per_m
            mobility_cm2_per_v_s = abs(signed_mobility_m2_per_v_s) * 1.0e4
        else:
            signed_mobility_m2_per_v_s = float("nan")
            mobility_cm2_per_v_s = float("nan")

        records.append(
            {
                "field_V_per_cm": field_v_per_cm,
                "field_V_per_m": field_v_per_m,
                "field_abs_V_per_cm": abs(field_v_per_cm),
                "field_abs_V_per_m": abs(field_v_per_m),
                "mean_velocity_x_m_per_s": velocity_x_m_per_s,
                "mean_velocity_abs_m_per_s": abs(velocity_x_m_per_s),
                "signed_mobility_m2_per_V_s": signed_mobility_m2_per_v_s,
                "signed_mobility_cm2_per_V_s": signed_mobility_m2_per_v_s * 1.0e4,
                "mobility_m2_per_V_s": abs(signed_mobility_m2_per_v_s),
                "mobility_cm2_per_V_s": mobility_cm2_per_v_s,
                "mean_energy_eV": stats["mean_energy_eV"],
                "mean_ionization_coeff_1_per_m": stats["mean_ionization_coeff_1_per_m"],
                "mean_ionization_coeff_cm_1": stats["mean_ionization_coeff_1_per_m"] / 100.0,
                "discarded_carriers_over_max_energy": stats["discarded_carriers_over_max_energy"],
                "run_complete": stats["run_complete"],
                "runtime_s": elapsed,
                "stats_file": str(stats_file),
            }
        )

    return pd.DataFrame.from_records(records).sort_values("field_V_per_cm")


def select_fit_data(df: pd.DataFrame, max_field_v_per_cm: float) -> pd.DataFrame:
    data = df[df["field_V_per_m"].abs() > 0.0].copy()
    data = data[data["field_V_per_cm"].abs() <= max_field_v_per_cm]
    data = data.sort_values("field_V_per_cm")

    if len(data) < 2:
        raise RuntimeError(
            "Cannot extract mobility: at least two non-zero fit points are required. "
            "Increase --mobility-fit-max-field or provide more low-field points."
        )

    return data


def extract_low_field_mobility(
    df: pd.DataFrame,
    max_field_v_per_cm: float,
) -> tuple[float, float, float, pd.DataFrame]:
    fit_data = select_fit_data(df, max_field_v_per_cm)

    field = fit_data["field_V_per_m"].to_numpy(dtype=float)
    velocity = fit_data["mean_velocity_x_m_per_s"].to_numpy(dtype=float)

    denominator = float(np.sum(field * field))
    if denominator <= 0.0:
        raise RuntimeError("Cannot extract mobility: zero field denominator.")

    mobility_m2_per_v_s = float(np.sum(field * velocity) / denominator)
    mobility_cm2_per_v_s = abs(mobility_m2_per_v_s) * 1.0e4
    intercept_m_per_s = 0.0

    free_slope_m2_per_v_s, free_intercept_m_per_s = np.polyfit(field, velocity, deg=1)

    fit_data = fit_data.copy()
    fit_data["fitted_velocity_x_m_per_s"] = mobility_m2_per_v_s * fit_data["field_V_per_m"]
    fit_data["fitted_mobility_cm2_per_V_s"] = mobility_cm2_per_v_s
    fit_data["fit_intercept_m_per_s"] = intercept_m_per_s
    fit_data["free_slope_diagnostic_cm2_per_V_s"] = abs(float(free_slope_m2_per_v_s)) * 1.0e4
    fit_data["free_intercept_diagnostic_m_per_s"] = float(free_intercept_m_per_s)

    return mobility_m2_per_v_s, mobility_cm2_per_v_s, intercept_m_per_s, fit_data


def build_fit_curve(
    fit_data: pd.DataFrame,
    mobility_m2_per_v_s: float,
    intercept_m_per_s: float,
) -> pd.DataFrame:
    min_field = float(fit_data["field_V_per_cm"].min())
    max_field = float(fit_data["field_V_per_cm"].max())
    field_v_per_cm = np.linspace(min_field, max_field, 200)
    field_v_per_m = field_v_per_cm * 100.0

    return pd.DataFrame(
        {
            "field_V_per_cm": field_v_per_cm,
            "field_V_per_m": field_v_per_m,
            "fitted_velocity_x_m_per_s": mobility_m2_per_v_s * field_v_per_m + intercept_m_per_s,
        }
    )


def plot_velocity(
    df: pd.DataFrame,
    fit_data: pd.DataFrame,
    fit_curve: pd.DataFrame,
    mobility_cm2_per_v_s: float,
    outdir: Path,
    show: bool,
) -> None:
    fig, ax = plt.subplots()
    ax.plot(df["field_V_per_cm"], df["mean_velocity_x_m_per_s"], marker="o", label="FBMC signed vx")
    ax.plot(
        fit_curve["field_V_per_cm"],
        fit_curve["fitted_velocity_x_m_per_s"],
        linestyle="--",
        label=f"zero-intercept fit: mu = {mobility_cm2_per_v_s:.1f} cm^2/V/s",
    )
    ax.scatter(fit_data["field_V_per_cm"], fit_data["mean_velocity_x_m_per_s"], marker="s", label="Fit points")
    ax.axhline(0.0, linewidth=0.8)
    ax.axvline(0.0, linewidth=0.8)
    ax.set_xlabel("Electric field Ex (V/cm)")
    ax.set_ylabel("Mean drift velocity vx (m/s)")
    ax.set_title("Bulk FBMC signed drift velocity versus electric field")
    ax.grid(True, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "velocity_vs_field.png", dpi=200)
    fig.savefig(outdir / "velocity_vs_field.pdf")
    if show:
        plt.show()
    plt.close(fig)


def plot_mobility(
    df: pd.DataFrame,
    fit_data: pd.DataFrame,
    mobility_cm2_per_v_s: float,
    outdir: Path,
    show: bool,
) -> None:
    fig, ax = plt.subplots()
    ax.plot(df["field_V_per_cm"], df["mobility_cm2_per_V_s"], marker="o", label="Pointwise mobility")
    ax.axhline(
        mobility_cm2_per_v_s,
        linestyle="--",
        label=f"zero-intercept fit mu = {mobility_cm2_per_v_s:.1f} cm^2/V/s",
    )
    ax.scatter(fit_data["field_V_per_cm"], fit_data["mobility_cm2_per_V_s"], marker="s", label="Fit points")
    ax.set_xlabel("Electric field Ex (V/cm)")
    ax.set_ylabel("Mobility (cm^2/V/s)")
    ax.set_title("Bulk FBMC mobility versus electric field")
    ax.grid(True, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "mobility_vs_field.png", dpi=200)
    fig.savefig(outdir / "mobility_vs_field.pdf")
    if show:
        plt.show()
    plt.close(fig)


def plot_energy(df: pd.DataFrame, outdir: Path, show: bool) -> None:
    fig, ax = plt.subplots()
    ax.plot(df["field_V_per_cm"], df["mean_energy_eV"], marker="o")
    ax.axvline(0.0, linewidth=0.8)
    ax.set_xlabel("Electric field Ex (V/cm)")
    ax.set_ylabel("Mean energy (eV)")
    ax.set_title("Bulk FBMC mean energy versus electric field")
    ax.grid(True, which="both")
    fig.tight_layout()
    fig.savefig(outdir / "energy_vs_field.png", dpi=200)
    fig.savefig(outdir / "energy_vs_field.pdf")
    if show:
        plt.show()
    plt.close(fig)


def write_summary(
    outdir: Path,
    mobility_m2_per_v_s: float,
    mobility_cm2_per_v_s: float,
    intercept_m_per_s: float,
    fit_data: pd.DataFrame,
) -> None:
    summary_file = outdir / "mobility_summary.txt"
    with summary_file.open("w", encoding="utf-8") as stream:
        stream.write(f"low_field_mobility_m2_per_V_s = {mobility_m2_per_v_s:.8e}\n")
        stream.write(f"low_field_mobility_cm2_per_V_s = {mobility_cm2_per_v_s:.8e}\n")
        stream.write("mobility_extraction_method = signed_zero_intercept_vx_vs_requested_Ex\n")
        stream.write(f"linear_fit_intercept_m_per_s = {intercept_m_per_s:.8e}\n")
        stream.write(f"fit_field_min_V_per_cm = {float(fit_data['field_V_per_cm'].min()):.8e}\n")
        stream.write(f"fit_field_max_V_per_cm = {float(fit_data['field_V_per_cm'].max()):.8e}\n")
        stream.write(f"fit_points = {len(fit_data)}\n")


def main() -> int:
    args = parse_args()
    validate_args(args)
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = build_sweep_dataframe(args)
    mobility_m2_per_v_s, mobility_cm2_per_v_s, intercept_m_per_s, fit_data = extract_low_field_mobility(
        df,
        args.mobility_fit_max_field,
    )
    fit_curve = build_fit_curve(fit_data, mobility_m2_per_v_s, intercept_m_per_s)

    df["low_field_mobility_m2_per_V_s"] = mobility_m2_per_v_s
    df["low_field_mobility_cm2_per_V_s"] = mobility_cm2_per_v_s
    df["low_field_fit_intercept_m_per_s"] = intercept_m_per_s

    df.to_csv(args.outdir / "field_sweep_results.csv", index=False)
    fit_data.to_csv(args.outdir / "mobility_fit_points.csv", index=False)
    fit_curve.to_csv(args.outdir / "mobility_fit_curve.csv", index=False)

    write_summary(args.outdir, mobility_m2_per_v_s, mobility_cm2_per_v_s, intercept_m_per_s, fit_data)
    plot_velocity(df, fit_data, fit_curve, mobility_cm2_per_v_s, args.outdir, args.show)
    plot_mobility(df, fit_data, mobility_cm2_per_v_s, args.outdir, args.show)
    plot_energy(df, args.outdir, args.show)

    print(f"Extracted FBMC low-field mobility: {mobility_cm2_per_v_s:.3f} cm^2/V/s")
    print(f"Wrote {args.outdir / 'field_sweep_results.csv'}")
    print(f"Wrote {args.outdir / 'mobility_summary.txt'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
