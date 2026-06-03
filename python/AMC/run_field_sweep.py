#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import subprocess
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd



def van_overstraeten_de_man_alpha_n(electric_field_V_per_m, temperature_K=300.0):
    electric_field_V_per_cm = np.asarray(electric_field_V_per_m) / 100.0

    a_inf = 7.03e5  # cm^-1
    b = 1.231e6     # V/cm

    hbar_omega_op_eV = 0.063
    k_B_eV_per_K = 8.617333262145e-5
    T0 = 300.0

    gamma = np.tanh(hbar_omega_op_eV / (2.0 * k_B_eV_per_K * T0)) / \
            np.tanh(hbar_omega_op_eV / (2.0 * k_B_eV_per_K * temperature_K))

    return gamma * a_inf * np.exp(-gamma * b / electric_field_V_per_cm)



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a bulk AMC electric-field sweep and extract low-field mobility."
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
        "--enable-impurity-scattering",
        action="store_true",
        help="Enable impurity scattering in the bulk AMC executable.",
    )

    parser.add_argument(
        "--impurity-density",
        type=float,
        default=0.0,
        help="Uniform ionized impurity density in cm^-3.",
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

    parser.add_argument(
        "--mobility-fit-max-field",
        type=float,
        default=2e3,
        help=(
            "Maximum electric field in V/cm used for low-field mobility extraction. "
            "If omitted, the lowest third of non-zero field points is used."
        ),
    )
    
    parser.add_argument(
        "--nbthreads",
        type=int,
        default=1,
        help="Number of threads for parallel execution.",
    )

    parser.add_argument(
        "--show",
        action="store_true",
        help="Show plots interactively after saving them.",
    )

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not args.exe.exists():
        raise FileNotFoundError(f"Executable not found: {args.exe}")
    
    if args.impurity_density < 0.0:
        raise ValueError("--impurity-density must be non-negative.")

    if args.impurity_density > 0.0 and not args.enable_impurity_scattering:
        print(
            "Warning: --impurity-density is positive but "
            "--enable-impurity-scattering was not provided."
        )

    if args.npart <= 0:
        raise ValueError("--npart must be positive.")

    if args.time <= 0.0:
        raise ValueError("--time must be positive.")

    if args.dt <= 0.0:
        raise ValueError("--dt must be positive.")

    if args.temperature < 0.0:
        raise ValueError("--temperature must be non-negative.")

    if args.max_energy <= 0.0:
        raise ValueError("--max-energy must be positive.")

    if args.gamma_safety <= 0.0:
        raise ValueError("--gamma-safety must be positive.")

    if args.gamma_samples < 2:
        raise ValueError("--gamma-samples must be at least 2.")

    if args.warmup < 0.0 or args.warmup >= 1.0:
        raise ValueError("--warmup must be in [0, 1).")

    if len(args.fields) < 2:
        raise ValueError("At least two electric-field values are required.")

    if args.mobility_fit_max_field is not None and args.mobility_fit_max_field <= 0.0:
        raise ValueError("--mobility-fit-max-field must be positive.")


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
        "--impurity-density",
        str(args.impurity_density),
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
        "-j",
        str(args.nbthreads)
    ]
    if args.enable_impurity_scattering:
        command.append("--enable-impurity-scattering")

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

    required_columns = [
        "electric_field_V_per_m",
        "mean_velocity_x_m_per_s",
        "mean_kinetic_energy_eV",
        "sample_count",
        "impact_ionization_events",
        "impact_ionization_rate_per_carrier_s_1",
        "impact_ionization_drift_velocity_m_per_s",
        "impact_ionization_coefficient_cm_1",
    ]

    for column in required_columns:
        if column not in row:
            raise KeyError(f"Missing column '{column}' in {path}")

    return {
        "electric_field_V_per_m": float(row["electric_field_V_per_m"]),
        "mean_velocity_x_m_per_s": float(row["mean_velocity_x_m_per_s"]),
        "mean_kinetic_energy_eV": float(row["mean_kinetic_energy_eV"]),
        "sample_count": float(row["sample_count"]),
        "impact_ionization_events": float(row["impact_ionization_events"]),
        "impact_ionization_rate_per_carrier_s_1": float(row["impact_ionization_rate_per_carrier_s_1"]),
        "impact_ionization_drift_velocity_m_per_s": float(row["impact_ionization_drift_velocity_m_per_s"]),
        "impact_ionization_coefficient_cm_1": float(row["impact_ionization_coefficient_cm_1"]),
    }


def build_sweep_dataframe(args: argparse.Namespace) -> pd.DataFrame:
    records: list[dict[str, float | str]] = []

    for field_v_per_cm in args.fields:
        print(f"Running Ex = {field_v_per_cm:.6e} V/cm")

        started = time.perf_counter()
        observables_file = run_one_field(args, field_v_per_cm)
        elapsed = time.perf_counter() - started

        row = read_last_observable_row(observables_file)

        field_v_per_m = row["electric_field_V_per_m"]
        velocity_x_m_per_s = row["mean_velocity_x_m_per_s"]

        if field_v_per_m != 0.0:
            mobility_m2_per_v_s = abs(velocity_x_m_per_s) / abs(field_v_per_m)
        else:
            mobility_m2_per_v_s = float("nan")

        records.append(
            {
                "field_V_per_cm": field_v_per_cm,
                "field_V_per_m": field_v_per_m,
                "enable_impurity_scattering": args.enable_impurity_scattering,
                "impurity_density_cm_3": args.impurity_density,
                "mean_velocity_x_m_per_s": velocity_x_m_per_s,
                "mean_velocity_abs_m_per_s": abs(velocity_x_m_per_s),
                "mobility_m2_per_V_s": mobility_m2_per_v_s,
                "mobility_cm2_per_V_s": mobility_m2_per_v_s * 1.0e4,
                "mean_kinetic_energy_eV": row["mean_kinetic_energy_eV"],
                "sample_count": row["sample_count"],
                "impact_ionization_events": row["impact_ionization_events"],
                "impact_ionization_rate_per_carrier_s_1": row["impact_ionization_rate_per_carrier_s_1"],
                "impact_ionization_drift_velocity_m_per_s": row["impact_ionization_drift_velocity_m_per_s"],
                "impact_ionization_coefficient_cm_1": row["impact_ionization_coefficient_cm_1"],
                "inverse_field_cm_per_V": 1.0 / abs(field_v_per_cm) if field_v_per_cm != 0.0 else float("nan"),
                "runtime_s": elapsed,
                "observables_file": str(observables_file),
            }
        )

    return pd.DataFrame.from_records(records).sort_values("field_V_per_cm")


def select_fit_data(
    df: pd.DataFrame,
    max_field_v_per_cm: float | None,
) -> pd.DataFrame:
    data = df.copy()
    data = data[data["field_V_per_m"].abs() > 0.0]
    data = data.sort_values("field_V_per_cm")

    if data.empty:
        raise RuntimeError("Cannot extract mobility: all electric fields are zero.")

    if max_field_v_per_cm is not None:
        fit_data = data[data["field_V_per_cm"].abs() <= max_field_v_per_cm]
    else:
        n_fit = max(2, len(data) // 3)
        fit_data = data.head(n_fit)

    if len(fit_data) < 2:
        raise RuntimeError(
            "Cannot extract mobility: at least two fit points are required. "
            "Increase --mobility-fit-max-field or provide more low-field points."
        )

    return fit_data.copy()


def extract_low_field_mobility(
    df: pd.DataFrame,
    max_field_v_per_cm: float | None,
) -> tuple[float, float, float, pd.DataFrame]:
    fit_data = select_fit_data(df, max_field_v_per_cm)

    field = fit_data["field_V_per_m"].abs().to_numpy(dtype=float)
    velocity = fit_data["mean_velocity_abs_m_per_s"].to_numpy(dtype=float)

    mobility_m2_per_v_s, intercept_m_per_s = np.polyfit(field, velocity, deg=1)

    mobility_m2_per_v_s = float(mobility_m2_per_v_s)
    intercept_m_per_s = float(intercept_m_per_s)
    mobility_cm2_per_v_s = mobility_m2_per_v_s * 1.0e4

    fit_data["fitted_velocity_abs_m_per_s"] = (
        mobility_m2_per_v_s * fit_data["field_V_per_m"].abs()
        + intercept_m_per_s
    )
    fit_data["fitted_mobility_cm2_per_V_s"] = mobility_cm2_per_v_s
    fit_data["fit_intercept_m_per_s"] = intercept_m_per_s

    return mobility_m2_per_v_s, mobility_cm2_per_v_s, intercept_m_per_s, fit_data


def build_fit_curve(
    fit_data: pd.DataFrame,
    mobility_m2_per_v_s: float,
    intercept_m_per_s: float,
) -> pd.DataFrame:
    min_field = float(fit_data["field_V_per_cm"].min())
    max_field = float(fit_data["field_V_per_cm"].max())

    field_v_per_cm = np.logspace(
        np.log10(min_field),
        np.log10(max_field),
        200,
    )

    field_v_per_m = field_v_per_cm * 100.0
    velocity = mobility_m2_per_v_s * field_v_per_m + intercept_m_per_s

    return pd.DataFrame(
        {
            "field_V_per_cm": field_v_per_cm,
            "field_V_per_m": field_v_per_m,
            "fitted_velocity_abs_m_per_s": velocity,
        }
    )


def plot_velocity(
    df: pd.DataFrame,
    fit_data: pd.DataFrame,
    fit_curve: pd.DataFrame,
    mobility_cm2_per_v_s: float,
    intercept_m_per_s: float,
    outdir: Path,
    show: bool,
) -> None:
    fig, ax = plt.subplots()

    ax.plot(
        df["field_V_per_cm"],
        df["mean_velocity_abs_m_per_s"],
        marker="o",
        label="AMC data",
    )

    ax.plot(
        fit_curve["field_V_per_cm"],
        fit_curve["fitted_velocity_abs_m_per_s"],
        linestyle="--",
        label=(
            f"Linear fit: μ = {mobility_cm2_per_v_s:.1f} cm²/V/s, "
            f"b = {intercept_m_per_s:.2e} m/s"
        ),
    )

    ax.scatter(
        fit_data["field_V_per_cm"],
        fit_data["mean_velocity_abs_m_per_s"],
        marker="s",
        label="Fit points",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Electric field (V/cm)")
    ax.set_ylabel("|Mean drift velocity x| (m/s)")
    ax.set_title("Bulk AMC drift velocity versus electric field")
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

    ax.plot(
        df["field_V_per_cm"],
        df["mobility_cm2_per_V_s"],
        marker="o",
        label="Pointwise mobility",
    )

    ax.axhline(
        mobility_cm2_per_v_s,
        linestyle="--",
        label=f"Linear-fit μ = {mobility_cm2_per_v_s:.1f} cm²/V/s",
    )

    ax.scatter(
        fit_data["field_V_per_cm"],
        fit_data["mobility_cm2_per_V_s"],
        marker="s",
        label="Fit points",
    )

    ax.set_xscale("log")
    ax.set_xlabel("Electric field (V/cm)")
    ax.set_ylabel("Mobility (cm²/V/s)")
    ax.set_title("Bulk AMC mobility versus electric field")
    ax.grid(True, which="both")
    ax.legend()
    
    print(f"Mobility : {mobility_cm2_per_v_s:.1f} cm²/V/s")

    fig.tight_layout()
    fig.savefig(outdir / "mobility_vs_field.png", dpi=200)
    fig.savefig(outdir / "mobility_vs_field.pdf")

    if show:
        plt.show()

    plt.close(fig)


def plot_energy(
    df: pd.DataFrame,
    outdir: Path,
    show: bool,
) -> None:
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

    if show:
        plt.show()

    plt.close(fig)


def plot_impact_ionization_coefficient(
    df: pd.DataFrame,
    particle: str,
    temperature_K: float,
    outdir: Path,
    show: bool,
) -> None:
    data = df.copy()
    data = data.replace([np.inf, -np.inf], np.nan)
    data = data.dropna(
        subset=[
            "inverse_field_cm_per_V",
            "impact_ionization_coefficient_cm_1",
        ]
    )

    positive_data = data[
        (data["inverse_field_cm_per_V"] > 0.0)
        & (data["impact_ionization_coefficient_cm_1"] > 0.0)
    ].copy()

    fig, ax = plt.subplots()

    if not positive_data.empty:
        positive_data = positive_data.sort_values("inverse_field_cm_per_V")

        ax.plot(
            positive_data["inverse_field_cm_per_V"],
            positive_data["impact_ionization_coefficient_cm_1"],
            marker="o",
            label="AMC data",
        )

    if particle == "electron":
        field_data = data[data["field_V_per_m"].abs() > 0.0].copy()

        if not field_data.empty:
            min_field = float(field_data["field_V_per_m"].abs().min())
            max_field = float(field_data["field_V_per_m"].abs().max())

            if min_field > 0.0 and max_field > min_field:
                field_reference_V_per_m = np.logspace(
                    np.log10(min_field),
                    np.log10(max_field),
                    300,
                )
                inverse_field_reference_cm_per_V = 1.0 / (field_reference_V_per_m / 100.0)
                alpha_reference_cm_1 = van_overstraeten_de_man_alpha_n(
                    field_reference_V_per_m,
                    temperature_K,
                )

                positive_reference = alpha_reference_cm_1 > 0.0

                ax.plot(
                    inverse_field_reference_cm_per_V[positive_reference],
                    alpha_reference_cm_1[positive_reference],
                    linestyle="--",
                    label="Van Overstraeten-de Man electron",
                )

    ax.set_yscale("log")
    # ax.set_xlim(xmax=
    # ax.set_xlim(1.0e-4, 7.0e-4)
    ax.set_xlabel("1 / electric field (cm/V)")
    ax.set_ylabel("Impact ionization coefficient (cm$^{-1}$)")
    ax.set_title("Bulk AMC impact ionization coefficient")
    ax.grid(True, which="both")
    ax.legend()

    fig.tight_layout()
    fig.savefig(outdir / "impact_ionization_coefficient_vs_inverse_field.png", dpi=200)
    fig.savefig(outdir / "impact_ionization_coefficient_vs_inverse_field.pdf")

    if show:
        plt.show()

    plt.close(fig)

    if positive_data.empty:
        print(
            "Warning: no positive impact-ionization coefficients were available; "
            "the log-scale ionization plot contains no AMC points."
        )


def write_summary(
    outdir: Path,
    mobility_m2_per_v_s: float,
    mobility_cm2_per_v_s: float,
    intercept_m_per_s: float,
    fit_data: pd.DataFrame,
) -> None:
    summary_file = outdir / "mobility_summary.txt"

    min_fit_field = float(fit_data["field_V_per_cm"].min())
    max_fit_field = float(fit_data["field_V_per_cm"].max())
    n_fit_points = len(fit_data)

    with summary_file.open("w", encoding="utf-8") as stream:
        stream.write(f"low_field_mobility_m2_per_V_s = {mobility_m2_per_v_s:.8e}\n")
        stream.write(f"low_field_mobility_cm2_per_V_s = {mobility_cm2_per_v_s:.8e}\n")
        stream.write(f"linear_fit_intercept_m_per_s = {intercept_m_per_s:.8e}\n")
        stream.write(f"fit_field_min_V_per_cm = {min_fit_field:.8e}\n")
        stream.write(f"fit_field_max_V_per_cm = {max_fit_field:.8e}\n")
        stream.write(f"fit_points = {n_fit_points}\n")


def main() -> int:
    args = parse_args()
    validate_args(args)
    args.outdir.mkdir(parents=True, exist_ok=True)
    df = build_sweep_dataframe(args)
    mobility_m2_per_v_s, mobility_cm2_per_v_s, intercept_m_per_s, fit_data = (
        extract_low_field_mobility(
            df,
            args.mobility_fit_max_field,
        )
    )

    fit_curve = build_fit_curve(
        fit_data,
        mobility_m2_per_v_s,
        intercept_m_per_s,
    )

    df["low_field_mobility_m2_per_V_s"] = mobility_m2_per_v_s
    df["low_field_mobility_cm2_per_V_s"] = mobility_cm2_per_v_s
    df["low_field_fit_intercept_m_per_s"] = intercept_m_per_s

    results_csv = args.outdir / "field_sweep_results.csv"
    fit_csv = args.outdir / "mobility_fit_points.csv"
    fit_curve_csv = args.outdir / "mobility_fit_curve.csv"

    df.to_csv(results_csv, index=False)
    fit_data.to_csv(fit_csv, index=False)
    fit_curve.to_csv(fit_curve_csv, index=False)

    write_summary(
        args.outdir,
        mobility_m2_per_v_s,
        mobility_cm2_per_v_s,
        intercept_m_per_s,
        fit_data,
    )

    plot_velocity(
        df,
        fit_data,
        fit_curve,
        mobility_cm2_per_v_s,
        intercept_m_per_s,
        args.outdir,
        args.show,
    )

    plot_mobility(
        df,
        fit_data,
        mobility_cm2_per_v_s,
        args.outdir,
        args.show,
    )

    plot_energy(
        df,
        args.outdir,
        args.show,
    )

    plot_impact_ionization_coefficient(
        df,
        args.particle,
        args.temperature,
        args.outdir,
        args.show,
    )

    # print(f"Extracted low-field mobility: {mobility_cm2_per_v_s:.3f} cm²/V/s")
    # print(f"Linear-fit intercept: {intercept_m_per_s:.6e} m/s")
    # print(f"Wrote {results_csv}")
    # print(f"Wrote {fit_csv}")
    # print(f"Wrote {fit_curve_csv}")
    # print(f"Wrote {args.outdir / 'mobility_summary.txt'}")
    # print(f"Wrote {args.outdir / 'velocity_vs_field.png'}")
    # print(f"Wrote {args.outdir / 'mobility_vs_field.png'}")
    # print(f"Wrote {args.outdir / 'energy_vs_field.png'}")
    # print(f"Wrote {args.outdir / 'impact_ionization_coefficient_vs_inverse_field.png'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())