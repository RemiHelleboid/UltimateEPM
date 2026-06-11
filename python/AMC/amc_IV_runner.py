#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import os
import shlex
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
        "--config",
        required=True,
        type=Path,
        help="Base device AMC YAML configuration file.",
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
        help="Minimum voltage applied to the swept contact, in V.",
    )

    parser.add_argument(
        "--vmax",
        required=True,
        type=float,
        help="Maximum voltage applied to the swept contact, in V.",
    )

    parser.add_argument(
        "--vstep",
        required=True,
        type=float,
        help="Voltage step in V.",
    )

    parser.add_argument(
        "--swept-contact",
        choices=["anode", "cathode"],
        default="anode",
        help="Contact voltage varied by the sweep.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Base random seed. If omitted, run.seed from the YAML file is used.",
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
        default=None,
        help="Override run.threads for each simulation; otherwise use the YAML value.",
    )

    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse existing device_history.csv files when available.",
    )

    parser.add_argument(
        "--show",
        action="store_true",
        help="Show plots interactively.",
    )

    parser.add_argument(
        "--set",
        dest="config_overrides",
        action="append",
        default=[],
        metavar="PATH=VALUE",
        help=(
            "Additional YAML override applied to every run. Repeat for multiple "
            "settings. Sweep-owned settings take precedence."
        ),
    )

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not args.exe.is_file():
        raise FileNotFoundError(f"Executable not found: {args.exe}")

    if not os.access(args.exe, os.X_OK):
        raise PermissionError(f"Executable is not executable: {args.exe}")

    if not args.config.is_file():
        raise FileNotFoundError(f"AMC configuration not found: {args.config}")

    finite_values = {
        "--vmin": args.vmin,
        "--vmax": args.vmax,
        "--vstep": args.vstep,
        "--transient-fraction": args.transient_fraction,
    }
    for option, value in finite_values.items():
        if not math.isfinite(value):
            raise ValueError(f"{option} must be finite.")

    if args.vstep == 0.0:
        raise ValueError("--vstep must be non-zero.")

    if (args.vmax - args.vmin) * args.vstep < 0.0:
        raise ValueError("--vstep sign is inconsistent with --vmin and --vmax.")

    if not (0.0 <= args.transient_fraction < 1.0):
        raise ValueError("--transient-fraction must be in [0, 1).")

    if args.jobs <= 0:
        raise ValueError("--jobs must be positive.")

    if args.threads_per_run is not None and args.threads_per_run <= 0:
        raise ValueError("--threads-per-run must be positive.")

    if args.history_filename != "device_history.csv":
        raise ValueError(
            "--history-filename must be 'device_history.csv'; "
            "the simulator does not support a custom history filename."
        )

    for override in args.config_overrides:
        key, separator, value = override.partition("=")
        if not separator or not key or not value:
            raise ValueError(
                f"Invalid --set value '{override}'. Expected PATH=VALUE."
            )

        runner_owned_keys = {
            f"contacts.{args.swept_contact}_voltage_V",
            "run.output_directory",
        }
        if args.seed is not None:
            runner_owned_keys.add("run.seed")
        if args.threads_per_run is not None:
            runner_owned_keys.add("run.threads")

        if key in runner_owned_keys:
            raise ValueError(
                f"--set {key}=... conflicts with an IV runner-owned setting."
            )


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


def voltage_directory_name(contact: str, voltage: float) -> str:
    prefix = "Va" if contact == "anode" else "Vc"
    return f"{prefix}_{voltage:+.6e}_V".replace("+", "p").replace("-", "m")


def build_command(
    args: argparse.Namespace,
    voltage: float,
    run_dir: Path,
    seed: int | None,
) -> list[str]:
    command = [
        str(args.exe),
        "--config",
        str(args.config),
    ]

    for override in args.config_overrides:
        command.extend(["--set", override])

    swept_voltage_key = f"contacts.{args.swept_contact}_voltage_V"
    command.extend(["--set", f"{swept_voltage_key}={voltage}"])
    command.extend(["--set", f"run.output_directory={run_dir}"])

    if args.seed is not None:
        command.extend(["--set", f"run.seed={seed}"])

    if args.threads_per_run is not None:
        command.extend(["--set", f"run.threads={args.threads_per_run}"])

    return command


def run_one_voltage(
    args: argparse.Namespace,
    voltage: float,
    seed: int | None,
) -> Path:
    run_dir = args.outdir / voltage_directory_name(args.swept_contact, voltage)
    run_dir.mkdir(parents=True, exist_ok=True)

    history_file = run_dir / args.history_filename

    if args.resume and history_file.exists():
        print(
            f"Reusing {args.swept_contact} = {voltage:.6e} V",
            flush=True,
        )
        return history_file

    command = build_command(args, voltage, run_dir, seed)

    log_file = run_dir / "stdout.log"
    command_file = run_dir / "command.txt"

    command_file.write_text(shlex.join(command) + "\n", encoding="utf-8")

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
            f"Simulation failed for {args.swept_contact}={voltage:.6e} V. "
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
        "mean_current_A": mean_current,
        "std_current_A": std_current,
        "stderr_current_A": float(stderr_current),
        "mean_abs_current_A": abs(mean_current),
        "mean_current_electron_A": float(np.mean(current_e)),
        "mean_current_hole_A": float(np.mean(current_h)),
        "std_current_electron_A": float(np.std(current_e, ddof=1)),
        "std_current_hole_A": float(np.std(current_h, ddof=1)),
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
    optional_mean_final(result, df, steady, "anode_voltage_V")
    optional_mean_final(result, df, steady, "cathode_voltage_V")
    optional_mean_final(result, df, steady, "quench_bias_voltage_V")
    optional_mean_final(result, df, steady, "quench_device_current_A")
    optional_mean_final(result, df, steady, "quench_resistor_current_A")
    optional_mean_max(result, df, steady, "quench_voltage_drop_V")

    return result


def run_and_extract_one_voltage(
    args: argparse.Namespace,
    voltage: float,
    seed: int | None,
) -> dict[str, float | str]:
    print(f"Running {args.swept_contact} = {voltage:.6e} V", flush=True)

    started = time.perf_counter()

    history_file = run_one_voltage(args, voltage, seed)

    record = extract_iv_point(
        history_file,
        voltage,
        args.transient_fraction,
    )

    record["runtime_s"] = float(time.perf_counter() - started)
    if seed is not None:
        record["seed"] = seed

    print(f"Finished {args.swept_contact} = {voltage:.6e} V", flush=True)

    return record


def plot_iv_signed(df: pd.DataFrame, outdir: Path, show: bool) -> None:
    fig, ax = plt.subplots()

    ax.errorbar(
        df["voltage_V"],
        df["mean_current_A"],
        yerr=df["stderr_current_A"],
        marker="o",
        linestyle="-",
        capsize=3,
        label="total",
    )

    ax.plot(
        df["voltage_V"],
        df["mean_current_electron_A"],
        marker="s",
        linestyle="--",
        label="electron",
    )

    ax.plot(
        df["voltage_V"],
        df["mean_current_hole_A"],
        marker="^",
        linestyle="--",
        label="hole",
    )

    ax.set_xlabel("Swept contact voltage (V)")
    ax.set_ylabel("Mean Ramo current (A)")
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
    data = df[df["mean_abs_current_A"] > 0.0].copy()

    if data.empty:
        print("No positive absolute current values available for log plot.")
        return

    fig, ax = plt.subplots()

    ax.plot(
        data["voltage_V"],
        data["mean_abs_current_A"],
        marker="o",
        linestyle="-",
    )

    ax.set_yscale("log")
    ax.set_xlabel("Swept contact voltage (V)")
    ax.set_ylabel("|Mean Ramo current| (A)")
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

    ax.set_xlabel("Swept contact voltage (V)")
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

    ax.set_xlabel("Swept contact voltage (V)")
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

    ax.set_xlabel("Swept contact voltage (V)")
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
        stream.write(f"config = {args.config}\n")
        stream.write(f"outdir = {args.outdir}\n\n")

        stream.write(f"swept_contact = {args.swept_contact}\n")
        stream.write(f"vmin = {args.vmin:.8e}\n")
        stream.write(f"vmax = {args.vmax:.8e}\n")
        stream.write(f"vstep = {args.vstep:.8e}\n")
        stream.write("voltages = " + " ".join(f"{v:.8e}" for v in voltages) + "\n\n")

        stream.write(f"transient_fraction = {args.transient_fraction:.8e}\n")
        stream.write(f"jobs = {args.jobs}\n")
        stream.write(f"threads_per_run = {args.threads_per_run}\n")
        stream.write(f"resume = {args.resume}\n")
        stream.write(f"seed = {args.seed}\n")
        for override in args.config_overrides:
            stream.write(f"set = {override}\n")


def main() -> int:
    args = parse_args()
    validate_args(args)

    args.exe = args.exe.resolve()
    args.config = args.config.resolve()
    args.outdir = args.outdir.resolve()

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

    print(f"Swept contact: {args.swept_contact}", flush=True)
    print(
        f"Running with jobs={args.jobs}, threads_per_run="
        f"{args.threads_per_run if args.threads_per_run is not None else 'config'}",
        flush=True,
    )

    records: list[dict[str, float | str]] = []

    if args.jobs == 1:
        for index, voltage in enumerate(voltages):
            seed = args.seed + index if args.seed is not None else None
            records.append(
                run_and_extract_one_voltage(args, voltage, seed)
            )
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            futures = {
                executor.submit(
                    run_and_extract_one_voltage,
                    args,
                    voltage,
                    args.seed + index if args.seed is not None else None,
                ): voltage
                for index, voltage in enumerate(voltages)
            }

            for future in as_completed(futures):
                voltage = futures[future]

                try:
                    records.append(future.result())
                except Exception as exc:
                    raise RuntimeError(
                        f"Voltage run failed for {args.swept_contact} = "
                        f"{voltage:.6e} V"
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
