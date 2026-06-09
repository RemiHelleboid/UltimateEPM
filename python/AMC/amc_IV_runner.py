#!/usr/bin/env python3

from __future__ import annotations

import argparse
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a self-consistent device AMC voltage sweep and extract an I/V curve."
    )

    parser.add_argument(
        "--exe",
        required=True,
        type=Path,
        help="Path to the device AMC executable, for example ./apps/device_amc.epm.",
    )

    parser.add_argument(
        "--device-mesh",
        required=True,
        type=Path,
        help="Path to the device mesh .msh file.",
    )

    parser.add_argument(
        "--material-file",
        type=Path,
        default=None,
        help=(
            "Path to the material YAML file used by the Poisson solver. "
            "If omitted, the executable default is used."
        ),
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for the full I/V sweep.",
    )

    parser.add_argument(
        "--vmin",
        required=True,
        type=float,
        help="Minimum anode voltage in V.",
    )

    parser.add_argument(
        "--vmax",
        required=True,
        type=float,
        help="Maximum anode voltage in V.",
    )

    parser.add_argument(
        "--vstep",
        required=True,
        type=float,
        help="Voltage step in V.",
    )

    parser.add_argument(
        "--cathode-voltage",
        type=float,
        default=0.0,
        help="Cathode voltage in V.",
    )

    parser.add_argument(
        "--time",
        type=float,
        default=50.0e-12,
        help="Final simulation time in seconds.",
    )

    parser.add_argument(
        "--dt",
        type=float,
        default=1.0e-15,
        help="Device AMC timestep in seconds.",
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
        help="Maximum carrier energy in eV for gamma_max precomputation.",
    )

    parser.add_argument(
        "--gamma-safety",
        type=float,
        default=1.2,
        help="Self-scattering gamma safety factor.",
    )

    parser.add_argument(
        "--gamma-samples",
        type=int,
        default=1000,
        help="Number of gamma_max energy samples.",
    )

    parser.add_argument(
        "--poisson-frequency",
        type=int,
        default=10,
        help="Number of transport steps between two Poisson updates.",
    )

    parser.add_argument(
        "--x0",
        type=float,
        required=True,
        help="Initial particle x position in microns. Required by the C++ app.",
    )

    parser.add_argument(
        "--y0",
        type=float,
        default=0.0,
        help="Initial particle y position in microns.",
    )

    parser.add_argument(
        "--z0",
        type=float,
        default=0.0,
        help="Initial particle z position in microns.",
    )

    parser.add_argument(
        "--nelectrons",
        type=int,
        default=1,
        help="Initial number of electrons when not using only doping initialization.",
    )

    parser.add_argument(
        "--nholes",
        type=int,
        default=0,
        help="Initial number of holes when not using only doping initialization.",
    )

    parser.add_argument(
        "--max-particles",
        type=int,
        default=1_000_000,
        help="Hard maximum number of active particles.",
    )

    parser.add_argument(
        "--avalanche-voltage-drop",
        type=float,
        default=1.0,
        help="Absolute quench-circuit voltage drop used to detect avalanche, in volts.",
    )

    parser.add_argument(
        "--quench-high-field",
        type=float,
        default=1.0e5,
        help="Particle electric-field threshold that resets the quench quiet window, in V/cm.",
    )

    parser.add_argument(
        "--quench-quiet-time",
        type=float,
        default=1.0e-11,
        help="Required time without high-field particles or impact ionization after avalanche, in seconds.",
    )

    parser.add_argument(
        "--effective-depth",
        type=float,
        default=1.0,
        help="Effective physical depth for 2D simulation, in microns.",
    )

    parser.add_argument(
        "--particle-z-period",
        type=float,
        default=1.0,
        help="Numerical periodic z length for 2D particles, in microns.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Base random seed.",
    )

    parser.add_argument(
        "--transient-fraction",
        type=float,
        default=0.5,
        help="Fraction of each run discarded before current averaging.",
    )

    parser.add_argument(
        "--history-filename",
        default="device_history.csv",
        help="History CSV filename created in each voltage output directory.",
    )

    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of independent voltage simulations to run in parallel.",
    )

    parser.add_argument(
        "--threads-per-run",
        type=int,
        default=1,
        help="Number of C++ threads passed to each individual simulation.",
    )

    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse existing device_history.csv files when available.",
    )

    parser.add_argument(
        "--disable-impact-ionization",
        action="store_true",
        help="Pass --disable-impact-ionization to the simulator.",
    )

    parser.add_argument(
        "--disable-particle-creation",
        action="store_true",
        help="Pass --disable-particle-creation to the simulator.",
    )

    parser.add_argument(
        "--disable-doping-init-particles",
        action="store_true",
        help="Pass --disable-doping-init-particles to the simulator.",
    )

    parser.add_argument(
        "--keep-going-without-electrons",
        action="store_true",
        help="Pass --keep-going-without-electrons to the simulator.",
    )

    parser.add_argument(
        "--export-time-steps",
        action="store_true",
        help="Export VTK/VTP timesteps during each run.",
    )

    parser.add_argument(
        "--export-frequency",
        type=int,
        default=100,
        help="Export one timestep every N iterations.",
    )

    parser.add_argument(
        "--show",
        action="store_true",
        help="Show plots interactively.",
    )

    parser.add_argument(
        "--extra-args",
        nargs=argparse.REMAINDER,
        default=[],
        help="Additional arguments passed directly to the simulator.",
    )

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not args.exe.exists():
        raise FileNotFoundError(f"Executable not found: {args.exe}")

    if not args.device_mesh.exists():
        raise FileNotFoundError(f"Device mesh not found: {args.device_mesh}")

    if args.material_file is not None and not args.material_file.exists():
        raise FileNotFoundError(f"Material file not found: {args.material_file}")

    if args.vstep == 0.0:
        raise ValueError("--vstep must be non-zero.")

    if (args.vmax - args.vmin) * args.vstep < 0.0:
        raise ValueError("--vstep sign is inconsistent with --vmin and --vmax.")

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

    if args.poisson_frequency <= 0:
        raise ValueError("--poisson-frequency must be positive.")

    if args.nelectrons < 0:
        raise ValueError("--nelectrons must be non-negative.")

    if args.nholes < 0:
        raise ValueError("--nholes must be non-negative.")

    if args.max_particles <= 0:
        raise ValueError("--max-particles must be positive.")

    if args.avalanche_voltage_drop <= 0.0:
        raise ValueError("--avalanche-voltage-drop must be positive.")

    if args.quench_high_field <= 0.0:
        raise ValueError("--quench-high-field must be positive.")

    if args.quench_quiet_time <= 0.0:
        raise ValueError("--quench-quiet-time must be positive.")

    if args.effective_depth <= 0.0:
        raise ValueError("--effective-depth must be positive.")

    if args.particle_z_period <= 0.0:
        raise ValueError("--particle-z-period must be positive.")

    if not (0.0 <= args.transient_fraction < 1.0):
        raise ValueError("--transient-fraction must be in [0, 1).")

    if args.jobs <= 0:
        raise ValueError("--jobs must be positive.")

    if args.threads_per_run <= 0:
        raise ValueError("--threads-per-run must be positive.")

    if args.export_frequency <= 0:
        raise ValueError("--export-frequency must be positive.")


def build_voltage_list(vmin: float, vmax: float, vstep: float) -> list[float]:
    voltages: list[float] = []
    value = vmin

    if vstep > 0.0:
        while value <= vmax + 0.5 * abs(vstep):
            voltages.append(round(value, 12))
            value += vstep
    else:
        while value >= vmax - 0.5 * abs(vstep):
            voltages.append(round(value, 12))
            value += vstep

    return voltages


def voltage_directory_name(voltage: float) -> str:
    return f"Va_{voltage:+.6e}_V".replace("+", "p").replace("-", "m")


def build_command(args: argparse.Namespace, voltage: float, run_dir: Path) -> list[str]:
    command = [
        str(args.exe),
        "--device-mesh",
        str(args.device_mesh),
        "--material",
        "Si",
        "--outdir",
        str(run_dir),
        "--time",
        str(args.time),
        "--dt",
        str(args.dt),
        "--temperature",
        str(args.temperature),
        "--max-energy",
        str(args.max_energy),
        "--gamma-safety",
        str(args.gamma_safety),
        "--gamma-samples",
        str(args.gamma_samples),
        "--poisson-frequency",
        str(args.poisson_frequency),
        "--anode-voltage",
        str(voltage),
        "--cathode-voltage",
        str(args.cathode_voltage),
        "--x0",
        str(args.x0),
        "--y0",
        str(args.y0),
        "--z0",
        str(args.z0),
        "--nelectrons",
        str(args.nelectrons),
        "--nholes",
        str(args.nholes),
        "--max-particles",
        str(args.max_particles),
        "--avalanche-voltage-drop",
        str(args.avalanche_voltage_drop),
        "--quench-high-field",
        str(args.quench_high_field),
        "--quench-quiet-time",
        str(args.quench_quiet_time),
        "--effective-depth",
        str(args.effective_depth),
        "--particle-z-period",
        str(args.particle_z_period),
        "--seed",
        str(args.seed),
        "--export-frequency",
        str(args.export_frequency),
        "-j",
        str(args.threads_per_run),
    ]

    if args.material_file is not None:
        command.extend(["--material-file", str(args.material_file)])

    if args.disable_impact_ionization:
        command.append("--disable-impact-ionization")

    if args.disable_particle_creation:
        command.append("--disable-particle-creation")

    if args.disable_doping_init_particles:
        command.append("--disable-doping-init-particles")

    if args.keep_going_without_electrons:
        command.append("--keep-going-without-electrons")

    if args.export_time_steps:
        command.append("--export-time-steps")

    command.extend(args.extra_args)

    return command


def run_one_voltage(args: argparse.Namespace, voltage: float) -> Path:
    run_dir = args.outdir / voltage_directory_name(voltage)
    run_dir.mkdir(parents=True, exist_ok=True)

    history_file = run_dir / args.history_filename

    if args.resume and history_file.exists():
        print(f"Reusing Va = {voltage:.6e} V", flush=True)
        return history_file

    command = build_command(args, voltage, run_dir)

    log_file = run_dir / "stdout.log"
    command_file = run_dir / "command.txt"

    command_file.write_text(" ".join(command) + "\n", encoding="utf-8")

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
            f"Simulation failed for Va={voltage:.6e} V. "
            f"See log: {log_file}"
        )

    if not history_file.exists():
        raise FileNotFoundError(
            f"Simulation completed but did not create {history_file}"
        )

    return history_file


def require_columns(df: pd.DataFrame, columns: list[str], filename: Path) -> None:
    for column in columns:
        if column not in df.columns:
            raise KeyError(f"Missing column '{column}' in {filename}")


def optional_mean_final(
    result: dict[str, float | str],
    df: pd.DataFrame,
    steady: pd.DataFrame,
    column: str,
) -> None:
    if column not in df.columns:
        return

    result[f"mean_{column}"] = float(steady[column].mean())
    result[f"final_{column}"] = float(df[column].iloc[-1])


def optional_mean_max(
    result: dict[str, float | str],
    df: pd.DataFrame,
    steady: pd.DataFrame,
    column: str,
) -> None:
    if column not in df.columns:
        return

    result[f"mean_{column}"] = float(steady[column].mean())
    result[f"max_{column}"] = float(df[column].max())
    result[f"final_{column}"] = float(df[column].iloc[-1])


def extract_iv_point(
    history_file: Path,
    voltage: float,
    transient_fraction: float,
) -> dict[str, float | str]:
    df = pd.read_csv(history_file)

    if df.empty:
        raise RuntimeError(f"History file is empty: {history_file}")

    required_columns = [
        "time",
        "nb_electrons",
        "nb_holes",
        "ramo_current_electron",
        "ramo_current_hole",
        "ramo_current",
    ]

    require_columns(df, required_columns, history_file)

    t_min = float(df["time"].min())
    t_max = float(df["time"].max())
    t_cut = t_min + transient_fraction * (t_max - t_min)

    steady = df[df["time"] >= t_cut].copy()

    if len(steady) < 2:
        raise RuntimeError(
            f"Not enough steady samples in {history_file}. "
            f"Decrease --transient-fraction or run longer."
        )

    current = steady["ramo_current"].to_numpy(dtype=float)
    current_e = steady["ramo_current_electron"].to_numpy(dtype=float)
    current_h = steady["ramo_current_hole"].to_numpy(dtype=float)

    mean_current = float(np.mean(current))
    std_current = float(np.std(current, ddof=1))
    stderr_current = std_current / float(np.sqrt(len(current)))

    result: dict[str, float | str] = {
        "voltage_V": float(voltage),

        "time_min_s": t_min,
        "time_max_s": t_max,
        "transient_cut_s": t_cut,
        "n_samples_total": int(len(df)),
        "n_samples_steady": int(len(steady)),

        "mean_current_A_per_um": mean_current,
        "std_current_A_per_um": std_current,
        "stderr_current_A_per_um": float(stderr_current),
        "mean_abs_current_A_per_um": abs(mean_current),

        "mean_current_electron_A_per_um": float(np.mean(current_e)),
        "mean_current_hole_A_per_um": float(np.mean(current_h)),
        "std_current_electron_A_per_um": float(np.std(current_e, ddof=1)),
        "std_current_hole_A_per_um": float(np.std(current_h, ddof=1)),

        "mean_nb_electrons": float(steady["nb_electrons"].mean()),
        "mean_nb_holes": float(steady["nb_holes"].mean()),
        "final_nb_electrons": float(df["nb_electrons"].iloc[-1]),
        "final_nb_holes": float(df["nb_holes"].iloc[-1]),
        "max_nb_electrons": float(df["nb_electrons"].max()),
        "max_nb_holes": float(df["nb_holes"].max()),

        "history_file": str(history_file),
    }

    optional_mean_final(result, df, steady, "total_electron_weight")
    optional_mean_final(result, df, steady, "total_hole_weight")

    optional_mean_max(result, df, steady, "nb_impact_ionization")
    optional_mean_max(result, df, steady, "max_electric_field")

    return result


def run_and_extract_one_voltage(
    args: argparse.Namespace,
    voltage: float,
) -> dict[str, float | str]:
    print(f"Running Va = {voltage:.6e} V", flush=True)

    started = time.perf_counter()

    history_file = run_one_voltage(args, voltage)

    record = extract_iv_point(
        history_file,
        voltage,
        args.transient_fraction,
    )

    record["runtime_s"] = float(time.perf_counter() - started)

    print(f"Finished Va = {voltage:.6e} V", flush=True)

    return record


def plot_iv_signed(df: pd.DataFrame, outdir: Path, show: bool) -> None:
    fig, ax = plt.subplots()

    ax.errorbar(
        df["voltage_V"],
        df["mean_current_A_per_um"],
        yerr=df["stderr_current_A_per_um"],
        marker="o",
        linestyle="-",
        capsize=3,
        label="total",
    )

    ax.plot(
        df["voltage_V"],
        df["mean_current_electron_A_per_um"],
        marker="s",
        linestyle="--",
        label="electron",
    )

    ax.plot(
        df["voltage_V"],
        df["mean_current_hole_A_per_um"],
        marker="^",
        linestyle="--",
        label="hole",
    )

    ax.set_xlabel("Anode voltage (V)")
    ax.set_ylabel("Mean Ramo current (A/um)")
    ax.set_title("Device AMC I/V curve")
    ax.grid(True)
    ax.legend()

    fig.tight_layout()
    fig.savefig(outdir / "iv_curve_signed.png", dpi=200)
    fig.savefig(outdir / "iv_curve_signed.pdf")

    if show:
        plt.show()

    plt.close(fig)


def plot_iv_abs_log(df: pd.DataFrame, outdir: Path, show: bool) -> None:
    data = df[df["mean_abs_current_A_per_um"] > 0.0].copy()

    if data.empty:
        print("No positive absolute current values available for log plot.")
        return

    fig, ax = plt.subplots()

    ax.plot(
        data["voltage_V"],
        data["mean_abs_current_A_per_um"],
        marker="o",
        linestyle="-",
    )

    ax.set_yscale("log")
    ax.set_xlabel("Anode voltage (V)")
    ax.set_ylabel("|Mean Ramo current| (A/um)")
    ax.set_title("Device AMC I/V curve")
    ax.grid(True, which="both")

    fig.tight_layout()
    fig.savefig(outdir / "iv_curve_abs_log.png", dpi=200)
    fig.savefig(outdir / "iv_curve_abs_log.pdf")

    if show:
        plt.show()

    plt.close(fig)


def plot_particle_counts(df: pd.DataFrame, outdir: Path, show: bool) -> None:
    fig, ax = plt.subplots()

    ax.plot(
        df["voltage_V"],
        df["mean_nb_electrons"],
        marker="o",
        label="mean electrons",
    )

    ax.plot(
        df["voltage_V"],
        df["mean_nb_holes"],
        marker="s",
        label="mean holes",
    )

    ax.set_xlabel("Anode voltage (V)")
    ax.set_ylabel("Mean numerical particle count")
    ax.set_title("Mean particle population versus voltage")
    ax.grid(True)
    ax.legend()

    fig.tight_layout()
    fig.savefig(outdir / "particle_count_vs_voltage.png", dpi=200)
    fig.savefig(outdir / "particle_count_vs_voltage.pdf")

    if show:
        plt.show()

    plt.close(fig)


def plot_represented_carriers_if_available(
    df: pd.DataFrame,
    outdir: Path,
    show: bool,
) -> None:
    columns = [
        "mean_total_electron_weight",
        "mean_total_hole_weight",
    ]

    if any(column not in df.columns for column in columns):
        return

    fig, ax = plt.subplots()

    ax.plot(
        df["voltage_V"],
        df["mean_total_electron_weight"],
        marker="o",
        label="electrons",
    )

    ax.plot(
        df["voltage_V"],
        df["mean_total_hole_weight"],
        marker="s",
        label="holes",
    )

    ax.set_xlabel("Anode voltage (V)")
    ax.set_ylabel("Mean represented carrier count")
    ax.set_title("Mean represented carriers versus voltage")
    ax.grid(True)
    ax.legend()

    fig.tight_layout()
    fig.savefig(outdir / "represented_carriers_vs_voltage.png", dpi=200)
    fig.savefig(outdir / "represented_carriers_vs_voltage.pdf")

    if show:
        plt.show()

    plt.close(fig)


def plot_max_field_if_available(
    df: pd.DataFrame,
    outdir: Path,
    show: bool,
) -> None:
    if "mean_max_electric_field" not in df.columns:
        return

    fig, ax = plt.subplots()

    ax.plot(
        df["voltage_V"],
        df["mean_max_electric_field"],
        marker="o",
    )

    ax.set_xlabel("Anode voltage (V)")
    ax.set_ylabel("Mean max electric field")
    ax.set_title("Maximum electric field versus voltage")
    ax.grid(True)

    fig.tight_layout()
    fig.savefig(outdir / "max_electric_field_vs_voltage.png", dpi=200)
    fig.savefig(outdir / "max_electric_field_vs_voltage.pdf")

    if show:
        plt.show()

    plt.close(fig)


def write_manifest(args: argparse.Namespace, voltages: list[float]) -> None:
    manifest_file = args.outdir / "iv_sweep_manifest.txt"

    with manifest_file.open("w", encoding="utf-8") as stream:
        stream.write("Device AMC I/V sweep\n")
        stream.write("====================\n\n")

        stream.write(f"exe = {args.exe}\n")
        stream.write(f"device_mesh = {args.device_mesh}\n")
        stream.write(f"material_file = {args.material_file}\n")
        stream.write(f"outdir = {args.outdir}\n\n")

        stream.write(f"vmin = {args.vmin:.8e}\n")
        stream.write(f"vmax = {args.vmax:.8e}\n")
        stream.write(f"vstep = {args.vstep:.8e}\n")
        stream.write("voltages = " + " ".join(f"{v:.8e}" for v in voltages) + "\n\n")

        stream.write(f"cathode_voltage = {args.cathode_voltage:.8e}\n")
        stream.write(f"time = {args.time:.8e}\n")
        stream.write(f"dt = {args.dt:.8e}\n")
        stream.write(f"temperature = {args.temperature:.8e}\n")
        stream.write(f"max_energy = {args.max_energy:.8e}\n")
        stream.write(f"gamma_safety = {args.gamma_safety:.8e}\n")
        stream.write(f"gamma_samples = {args.gamma_samples}\n")
        stream.write(f"poisson_frequency = {args.poisson_frequency}\n")
        stream.write(f"transient_fraction = {args.transient_fraction:.8e}\n\n")

        stream.write(f"x0 = {args.x0:.8e}\n")
        stream.write(f"y0 = {args.y0:.8e}\n")
        stream.write(f"z0 = {args.z0:.8e}\n")
        stream.write(f"nelectrons = {args.nelectrons}\n")
        stream.write(f"nholes = {args.nholes}\n")
        stream.write(f"max_particles = {args.max_particles}\n")
        stream.write(f"avalanche_threshold = {args.avalanche_threshold}\n")
        stream.write(f"effective_depth = {args.effective_depth:.8e}\n")
        stream.write(f"particle_z_period = {args.particle_z_period:.8e}\n\n")

        stream.write(f"jobs = {args.jobs}\n")
        stream.write(f"threads_per_run = {args.threads_per_run}\n")
        stream.write(f"resume = {args.resume}\n")
        stream.write(f"seed = {args.seed}\n\n")

        stream.write(f"disable_impact_ionization = {args.disable_impact_ionization}\n")
        stream.write(f"disable_particle_creation = {args.disable_particle_creation}\n")
        stream.write(f"disable_doping_init_particles = {args.disable_doping_init_particles}\n")
        stream.write(f"keep_going_without_electrons = {args.keep_going_without_electrons}\n")
        stream.write(f"export_time_steps = {args.export_time_steps}\n")
        stream.write(f"export_frequency = {args.export_frequency}\n")
        stream.write("extra_args = " + " ".join(args.extra_args) + "\n")


def main() -> int:
    args = parse_args()
    validate_args(args)

    args.outdir.mkdir(parents=True, exist_ok=True)

    voltages = build_voltage_list(
        args.vmin,
        args.vmax,
        args.vstep,
    )

    if not voltages:
        raise RuntimeError("Voltage list is empty.")

    write_manifest(args, voltages)

    print("Voltage sweep:")
    for voltage in voltages:
        print(f"  {voltage:.6e} V")

    print(
        f"Running with jobs={args.jobs}, "
        f"threads_per_run={args.threads_per_run}",
        flush=True,
    )

    records: list[dict[str, float | str]] = []

    if args.jobs == 1:
        for voltage in voltages:
            records.append(run_and_extract_one_voltage(args, voltage))
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            futures = {
                executor.submit(run_and_extract_one_voltage, args, voltage): voltage
                for voltage in voltages
            }

            for future in as_completed(futures):
                voltage = futures[future]

                try:
                    records.append(future.result())
                except Exception as exc:
                    raise RuntimeError(
                        f"Voltage run failed for Va = {voltage:.6e} V"
                    ) from exc

    df = pd.DataFrame.from_records(records)
    df = df.sort_values("voltage_V")

    results_file = args.outdir / "iv_curve_results.csv"
    df.to_csv(results_file, index=False)

    plot_iv_signed(df, args.outdir, args.show)
    plot_iv_abs_log(df, args.outdir, args.show)
    plot_particle_counts(df, args.outdir, args.show)
    plot_represented_carriers_if_available(df, args.outdir, args.show)
    plot_max_field_if_available(df, args.outdir, args.show)

    print(f"Wrote {results_file}")
    print(f"Wrote {args.outdir / 'iv_curve_signed.png'}")
    print(f"Wrote {args.outdir / 'iv_curve_abs_log.png'}")
    print(f"Wrote {args.outdir / 'particle_count_vs_voltage.png'}")

    if "mean_total_electron_weight" in df.columns:
        print(f"Wrote {args.outdir / 'represented_carriers_vs_voltage.png'}")

    if "mean_max_electric_field" in df.columns:
        print(f"Wrote {args.outdir / 'max_electric_field_vs_voltage.png'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
