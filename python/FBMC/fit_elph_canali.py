#!/usr/bin/env python3
"""Fit high-energy FBMC electron-phonon parameters to the Canali velocity curve.

For each trial profile this driver:
  1. reconstructs a full-band rate CSV from a reusable kernel;
  2. evaluates intrinsic MRTA mobility to preserve the low-field fit;
  3. runs bulk FBMC at fixed electric fields and a fixed random seed;
  4. compares projected drift-speed magnitudes with the zero-doping Canali curve.

Use the low-field Arora fit as --base-params. The kernel must cover the maximum
FBMC carrier energy; a 0.3 eV low-field kernel is normally insufficient.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
MOBILITY_PATTERN = re.compile(r"(?:μ_iso|mu_iso)\s*=\s*([0-9.eE+-]+)\s*cm\^2/\(V[·*]s\)")
SUPPORTED_PARAMETERS = (
    "acoustic-a",
    "optical-a",
    "acoustic-end",
    "optical-end",
    "threshold",
)
KERNEL_HEADER = [
    "vertex_index",
    "local_band_index",
    "energy_eV",
    "kernel_ac_L_ab",
    "kernel_ac_T_ab",
    "kernel_op_L_ab",
    "kernel_op_T_ab",
    "kernel_ac_L_em",
    "kernel_ac_T_em",
    "kernel_op_L_em",
    "kernel_op_T_em",
    "transport_kernel_ac",
    "transport_kernel_op",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--elph-exe", type=Path, default=REPO_ROOT / "build/apps/elph.epm")
    parser.add_argument("--fbmc-exe", type=Path, default=REPO_ROOT / "build/apps/fbmc.epm")
    parser.add_argument("--mesh", required=True, type=Path)
    parser.add_argument("--kernel", required=True, type=Path, help="Reusable electron kernel at --temperature.")
    parser.add_argument(
        "--base-params",
        required=True,
        type=Path,
        help="Low-field fitted YAML, normally fit_elph_arora/best.yaml.",
    )
    parser.add_argument(
        "--arora-params",
        type=Path,
        default=REPO_ROOT / "data/materials/Si/admc/arora-canali.yaml",
    )
    parser.add_argument("--output-dir", type=Path, default=REPO_ROOT / "fit_elph_canali_output")
    parser.add_argument("--material", default="Si")
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument(
        "--fields",
        default="10000,30000,100000,300000",
        help="Positive field magnitudes in V/cm, comma-separated.",
    )
    parser.add_argument(
        "--fit",
        default="acoustic-end,optical-end,threshold",
        help="Comma-separated parameters: " + ",".join(SUPPORTED_PARAMETERS),
    )
    parser.add_argument("--ncbands", type=int, default=2)
    parser.add_argument("--nvbands", type=int, default=4)
    parser.add_argument("--nthreads", type=int, default=16)
    parser.add_argument("--npart", type=int, default=300)
    parser.add_argument("--time", type=float, default=20.0e-12)
    parser.add_argument("--warmup", type=float, default=0.3)
    parser.add_argument("--max-energy", type=float, default=2.5)
    parser.add_argument("--gamma-safety", type=float, default=1.3)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--bz-domain", choices=("full", "octant"), default="full")
    parser.add_argument("--velocity-weight", type=float, default=1.0)
    parser.add_argument("--mobility-weight", type=float, default=0.25)
    parser.add_argument("--regularization", type=float, default=1.0e-3)
    parser.add_argument("--max-strength-factor", type=float, default=30.0)
    parser.add_argument("--threshold-min", type=float, default=0.05)
    parser.add_argument("--threshold-max", type=float, default=None)
    parser.add_argument("--maxiter", type=int, default=20)
    parser.add_argument("--xtol", type=float, default=5.0e-3)
    parser.add_argument("--evaluate-only", action="store_true")
    return parser.parse_args()


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as stream:
        value = yaml.safe_load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a YAML mapping in {path}")
    return value


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_numeric_tree(value):
    if isinstance(value, dict):
        return {key: normalize_numeric_tree(item) for key, item in value.items()}
    if isinstance(value, list):
        return [normalize_numeric_tree(item) for item in value]
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return value
    return value


def validate_executable(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")
    if not os.access(path, os.X_OK):
        raise PermissionError(f"{label} is not executable: {path}")


def validate_kernel_header(path: Path) -> None:
    with path.open("r", encoding="utf-8", newline="") as stream:
        header = next(csv.reader(stream), None)
    if header != KERNEL_HEADER:
        raise ValueError(f"{path} does not have the expected electron-phonon kernel header")


def read_and_validate_kernel_metadata(args: argparse.Namespace, base_config: dict) -> dict | None:
    metadata_path = Path(str(args.kernel) + ".meta.yaml")
    if not metadata_path.is_file():
        print(
            "warning: kernel metadata sidecar is absent; mesh, temperature, dispersion, and energy coverage "
            "cannot be verified",
            flush=True,
        )
        return None
    metadata = load_yaml(metadata_path)
    expected = {
        "model": "electron_phonon_kernel",
        "material": args.material,
        "carrier": "electron",
        "n_conduction_bands": args.ncbands,
        "n_valence_bands": args.nvbands,
        "bz_domain": args.bz_domain,
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(f"Kernel metadata {key}={metadata.get(key)!r}; expected {value!r}")
    if not math.isclose(float(metadata["temperature_K"]), args.temperature, abs_tol=1.0e-9):
        raise ValueError("Kernel temperature does not match --temperature")
    kernel_window = float(metadata["energy_window_eV"])
    if args.max_energy > kernel_window + 1.0e-12:
        raise ValueError(
            f"--max-energy={args.max_energy:g} eV exceeds kernel window {kernel_window:g} eV; "
            "generate a larger high-field kernel"
        )
    if int(metadata["mesh_size_bytes"]) != args.mesh.stat().st_size:
        raise ValueError("Kernel mesh metadata does not match --mesh")
    if normalize_numeric_tree(metadata.get("dispersion")) != normalize_numeric_tree(base_config.get("dispersion")):
        raise ValueError("Kernel dispersion does not match --base-params")
    if not math.isclose(float(metadata["Radius-WS"]), float(base_config["Radius-WS"]), rel_tol=1.0e-12):
        raise ValueError("Kernel Radius-WS does not match --base-params")
    return metadata


def validate_args(args: argparse.Namespace) -> None:
    validate_executable(args.elph_exe, "elph executable")
    validate_executable(args.fbmc_exe, "fbmc executable")
    for path, label in (
        (args.mesh, "mesh"),
        (args.kernel, "kernel"),
        (args.base_params, "base profile"),
        (args.arora_params, "Arora-Canali profile"),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"{label} not found: {path}")
    validate_kernel_header(args.kernel)

    args.fields = sorted({float(value) for value in args.fields.split(",") if value.strip()})
    if not args.fields or any(not math.isfinite(value) or value <= 0.0 for value in args.fields):
        raise ValueError("--fields must contain positive finite field magnitudes")
    args.fit = [name.strip() for name in args.fit.split(",") if name.strip()]
    unknown = sorted(set(args.fit) - set(SUPPORTED_PARAMETERS))
    if not args.fit or unknown or len(args.fit) != len(set(args.fit)):
        raise ValueError(f"Invalid --fit list; unsupported entries: {', '.join(unknown)}")
    if args.ncbands <= 0 or args.nvbands <= 0 or args.nthreads <= 0 or args.npart <= 0:
        raise ValueError("Band counts, threads, and particle count must be positive")
    if args.temperature <= 0.0 or args.time <= 0.0 or args.max_energy <= 0.0:
        raise ValueError("Temperature, simulation time, and maximum energy must be positive")
    if not 0.0 <= args.warmup < 1.0:
        raise ValueError("--warmup must be in [0, 1)")
    if args.gamma_safety < 1.0 or args.seed < 0:
        raise ValueError("--gamma-safety must be >= 1 and --seed must be non-negative")
    if min(args.velocity_weight, args.mobility_weight, args.regularization) < 0.0:
        raise ValueError("Objective weights must be non-negative")
    if args.velocity_weight == 0.0 and args.mobility_weight == 0.0:
        raise ValueError("At least one data objective weight must be positive")
    if args.max_strength_factor <= 1.0:
        raise ValueError("--max-strength-factor must be greater than one")
    if args.threshold_min <= 0.0:
        raise ValueError("--threshold-min must be positive")
    if args.threshold_max is None:
        args.threshold_max = args.max_energy
    if args.threshold_max < args.threshold_min:
        raise ValueError("--threshold-max must be >= --threshold-min")


def arora_zero_doping_mobility(config: dict, temperature_K: float) -> float:
    reference_temperature = float(config["reference_temperature_K"])
    low_field = config["electron"]["low_field"]
    ratio = temperature_K / reference_temperature
    return float(low_field["mu_min_300_cm2_per_V_s"]) * ratio ** float(
        low_field["mu_min_temperature_exponent"]
    ) + float(low_field["mu_dop_300_cm2_per_V_s"]) * ratio ** float(
        low_field["mu_dop_temperature_exponent"]
    )


def canali_velocity_m_per_s(config: dict, temperature_K: float, field_V_per_cm: float) -> float:
    mobility_cm2_per_V_s = arora_zero_doping_mobility(config, temperature_K)
    reference_temperature = float(config["reference_temperature_K"])
    high_field = config["electron"]["high_field"]
    ratio = temperature_K / reference_temperature
    saturation_velocity_cm_per_s = float(high_field["saturation_velocity_300_cm_per_s"]) * ratio ** float(
        high_field["saturation_velocity_temperature_exponent"]
    )
    beta = float(high_field["beta_300"]) * ratio ** float(high_field["beta_temperature_exponent"])
    velocity_ratio = mobility_cm2_per_V_s * field_V_per_cm / saturation_velocity_cm_per_s
    velocity_cm_per_s = mobility_cm2_per_V_s * field_V_per_cm / (
        1.0 + velocity_ratio**beta
    ) ** (1.0 / beta)
    return velocity_cm_per_s * 1.0e-2


def profile_parameters(config: dict) -> dict[str, float]:
    electron = config["deformation-potential"]["electron"]
    threshold = float(electron["energy-threshold"])
    acoustic_a = float(electron["acoustic"]["A"])
    optical_a = float(electron["optic"]["A"])
    return {
        "acoustic-a": acoustic_a,
        "optical-a": optical_a,
        "acoustic-end": acoustic_a + float(electron["acoustic"]["B"]) * threshold,
        "optical-end": optical_a + float(electron["optic"]["B"]) * threshold,
        "threshold": threshold,
    }


def write_trial_profile(
    base_config: dict,
    baseline: dict[str, float],
    fitted_names: list[str],
    log_values: list[float],
    path: Path,
) -> dict[str, float]:
    values = dict(baseline)
    for name, log_value in zip(fitted_names, log_values, strict=True):
        values[name] = math.exp(float(log_value))
    if any(not math.isfinite(value) or value <= 0.0 for value in values.values()):
        raise ValueError("Trial parameter is non-positive or non-finite")

    config = copy.deepcopy(base_config)
    config["parameter_set"] = "canali-fit-trial"
    electron = config["deformation-potential"]["electron"]
    threshold = values["threshold"]
    electron["energy-threshold"] = threshold
    electron["acoustic"]["A"] = values["acoustic-a"]
    electron["acoustic"]["B"] = (values["acoustic-end"] - values["acoustic-a"]) / threshold
    electron["optic"]["A"] = values["optical-a"]
    electron["optic"]["B"] = (values["optical-end"] - values["optical-a"]) / threshold
    with path.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(config, stream, sort_keys=False)
    return values


def run_command(command: list[str], cwd: Path, log_path: Path) -> None:
    with log_path.open("w", encoding="utf-8") as stream:
        completed = subprocess.run(command, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT, check=False, text=True)
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {completed.returncode}; see {log_path}")


def reconstruct_rates_and_mobility(
    args: argparse.Namespace,
    profile: Path,
    evaluation_dir: Path,
) -> tuple[Path, float]:
    rates_file = evaluation_dir / "phonon_rates.csv"
    log_path = evaluation_dir / "elph.log"
    command = [
        str(args.elph_exe.resolve()),
        "--meshbandfile",
        str(args.mesh.resolve()),
        "--material",
        args.material,
        "--carrier",
        "electron",
        "--phonon-params-file",
        str(profile.resolve()),
        "--ncbands",
        str(args.ncbands),
        "--nvbands",
        str(args.nvbands),
        "--temperature",
        f"{args.temperature:.12g}",
        "--energy_window",
        f"{args.max_energy:.12g}",
        "--nthreads",
        str(args.nthreads),
        "--bz-domain",
        args.bz_domain,
        "--kernel-file",
        str(args.kernel.resolve()),
        "--export-rates",
        "--rates-out",
        str(rates_file.resolve()),
        "--outdir",
        str((evaluation_dir / "elph_output").resolve()),
    ]
    run_command(command, evaluation_dir, log_path)
    text = log_path.read_text(encoding="utf-8")
    match = MOBILITY_PATTERN.search(text)
    if not match:
        raise RuntimeError(f"Could not extract MRTA mobility; see {log_path}")
    mobility = float(match.group(1))
    if not rates_file.is_file() or not math.isfinite(mobility) or mobility <= 0.0:
        raise RuntimeError("Invalid reconstructed rate file or MRTA mobility")
    return rates_file, mobility


def read_fbmc_observables(path: Path) -> tuple[float, float]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"No FBMC observables in {path}")
    row = rows[-1]
    velocity_key = (
        "drift_velocity_parallel_m_per_s"
        if "drift_velocity_parallel_m_per_s" in row
        else "mean_velocity_x_m_per_s"
    )
    velocity = abs(float(row[velocity_key]))
    energy = float(row["mean_kinetic_energy_eV"])
    if not math.isfinite(velocity) or velocity <= 0.0 or not math.isfinite(energy):
        raise RuntimeError(f"Invalid FBMC observables in {path}")
    return velocity, energy


def run_fbmc(
    args: argparse.Namespace,
    profile: Path,
    rates_file: Path,
    field_V_per_cm: float,
    evaluation_dir: Path,
) -> tuple[float, float]:
    field_dir = evaluation_dir / f"field_{field_V_per_cm:.6g}_Vcm"
    field_dir.mkdir(parents=True, exist_ok=True)
    command = [
        str(args.fbmc_exe.resolve()),
        "--meshbandfile",
        str(args.mesh.resolve()),
        "--phononfile",
        str(rates_file.resolve()),
        "--material",
        args.material,
        "--carrier",
        "electron",
        "--phonon-params-file",
        str(profile.resolve()),
        "--ncbands",
        str(args.ncbands),
        "--nvbands",
        str(args.nvbands),
        "--temperature",
        f"{args.temperature:.12g}",
        "--npart",
        str(args.npart),
        "--nthreads",
        str(args.nthreads),
        "--maxenergy",
        f"{args.max_energy:.12g}",
        "--gamma-safety",
        f"{args.gamma_safety:.12g}",
        "--time",
        f"{args.time:.12g}",
        "--warmup",
        f"{args.warmup:.12g}",
        "--Ex",
        f"{field_V_per_cm:.12g}",
        "--seed",
        str(args.seed),
        "--bz-domain",
        args.bz_domain,
        "--skip-mesh-vtk",
        "--outdir",
        str(field_dir.resolve()),
    ]
    run_command(command, field_dir, field_dir / "fbmc.log")
    return read_fbmc_observables(field_dir / "observables.csv")


class Objective:
    def __init__(self, args: argparse.Namespace, base_config: dict, arora_config: dict) -> None:
        self.args = args
        self.base_config = base_config
        self.arora_config = arora_config
        self.baseline = profile_parameters(base_config)
        self.output_dir = args.output_dir.resolve()
        self.trial_profile = self.output_dir / "trial.yaml"
        self.best_profile = self.output_dir / "best.yaml"
        self.best_rates = self.output_dir / "best_rates.csv"
        self.history = self.output_dir / "history.csv"
        self.best_loss = math.inf
        self.evaluation = 0
        self.cache: dict[tuple[float, ...], float] = {}
        for path in (self.trial_profile, self.best_profile, self.best_rates):
            path.unlink(missing_ok=True)
        with self.history.open("w", newline="", encoding="utf-8") as stream:
            csv.writer(stream).writerow(
                [
                    "evaluation",
                    *SUPPORTED_PARAMETERS,
                    "field_V_per_cm",
                    "target_velocity_m_per_s",
                    "model_velocity_m_per_s",
                    "mean_energy_eV",
                    "velocity_log_residual",
                    "mrta_target_cm2_per_V_s",
                    "mrta_model_cm2_per_V_s",
                    "velocity_loss",
                    "mobility_loss",
                    "regularization_loss",
                    "total_loss",
                ]
            )

    def __call__(self, log_values) -> float:
        values = [float(value) for value in log_values]
        key = tuple(round(value, 12) for value in values)
        if key in self.cache:
            return self.cache[key]
        self.evaluation += 1
        evaluation_dir = self.output_dir / "runs" / f"eval_{self.evaluation:04d}"
        evaluation_dir.mkdir(parents=True, exist_ok=True)
        parameters = write_trial_profile(
            self.base_config,
            self.baseline,
            self.args.fit,
            values,
            self.trial_profile,
        )
        shutil.copy2(self.trial_profile, evaluation_dir / "trial.yaml")
        rows = []
        try:
            rates_file, mrta_model = reconstruct_rates_and_mobility(
                self.args, self.trial_profile, evaluation_dir
            )
            for field in self.args.fields:
                model_velocity, mean_energy = run_fbmc(
                    self.args, self.trial_profile, rates_file, field, evaluation_dir
                )
                target_velocity = canali_velocity_m_per_s(self.arora_config, self.args.temperature, field)
                residual = math.log(model_velocity / target_velocity)
                rows.append((field, target_velocity, model_velocity, mean_energy, residual))
        except (OSError, RuntimeError, ValueError) as error:
            (evaluation_dir / "failure.txt").write_text(str(error) + "\n", encoding="utf-8")
            print(f"evaluation {self.evaluation:04d} failed: {error}", flush=True)
            self.cache[key] = 1.0e12
            return 1.0e12

        velocity_loss = sum(row[4] ** 2 for row in rows) / len(rows)
        mrta_target = arora_zero_doping_mobility(self.arora_config, self.args.temperature)
        mobility_loss = math.log(mrta_model / mrta_target) ** 2
        regularization_loss = self.args.regularization * sum(
            (value - math.log(self.baseline[name])) ** 2
            for name, value in zip(self.args.fit, values, strict=True)
        )
        total_loss = (
            self.args.velocity_weight * velocity_loss
            + self.args.mobility_weight * mobility_loss
            + regularization_loss
        )

        with self.history.open("a", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            for field, target, model, energy, residual in rows:
                writer.writerow(
                    [
                        self.evaluation,
                        *[f"{parameters[name]:.12g}" for name in SUPPORTED_PARAMETERS],
                        f"{field:.12g}",
                        f"{target:.12g}",
                        f"{model:.12g}",
                        f"{energy:.12g}",
                        f"{residual:.12g}",
                        f"{mrta_target:.12g}",
                        f"{mrta_model:.12g}",
                        f"{velocity_loss:.12g}",
                        f"{mobility_loss:.12g}",
                        f"{regularization_loss:.12g}",
                        f"{total_loss:.12g}",
                    ]
                )

        if total_loss < self.best_loss:
            self.best_loss = total_loss
            shutil.copy2(self.trial_profile, self.best_profile)
            shutil.copy2(rates_file, self.best_rates)
        text = " ".join(f"{name}={parameters[name]:.5g}" for name in self.args.fit)
        print(
            f"evaluation {self.evaluation:04d} loss={total_loss:.6g} "
            f"velocity={velocity_loss:.4g} mobility={mobility_loss:.4g} {text}",
            flush=True,
        )
        self.cache[key] = total_loss
        return total_loss


def parameter_bounds(args: argparse.Namespace, baseline: dict[str, float]) -> list[tuple[float, float]]:
    factor_log = math.log(args.max_strength_factor)
    bounds = []
    for name in args.fit:
        if name == "threshold":
            bounds.append((math.log(args.threshold_min), math.log(args.threshold_max)))
        else:
            center = math.log(baseline[name])
            bounds.append((center - factor_log, center + factor_log))
    return bounds


def main() -> int:
    args = parse_args()
    validate_args(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    base_config = load_yaml(args.base_params)
    arora_config = load_yaml(args.arora_params)
    if base_config.get("material") != args.material or base_config.get("model") != "electron_phonon":
        raise ValueError("Invalid base electron-phonon profile")
    if arora_config.get("material") != args.material or arora_config.get("model") != "arora_canali_mobility":
        raise ValueError("Invalid Arora-Canali profile")
    metadata = read_and_validate_kernel_metadata(args, base_config)
    baseline = profile_parameters(base_config)
    if any(not math.isfinite(value) or value <= 0.0 for value in baseline.values()):
        raise ValueError("Base profile contains a non-positive fitted strength")
    if "threshold" in args.fit and not args.threshold_min <= baseline["threshold"] <= args.threshold_max:
        raise ValueError(
            f"Base threshold {baseline['threshold']:g} eV is outside the optimization bounds "
            f"[{args.threshold_min:g}, {args.threshold_max:g}] eV"
        )

    manifest = {
        "elph_executable": str(args.elph_exe.resolve()),
        "fbmc_executable": str(args.fbmc_exe.resolve()),
        "mesh": str(args.mesh.resolve()),
        "mesh_sha256": sha256_file(args.mesh),
        "kernel": str(args.kernel.resolve()),
        "kernel_sha256": sha256_file(args.kernel),
        "kernel_metadata": metadata,
        "base_params": str(args.base_params.resolve()),
        "base_params_sha256": sha256_file(args.base_params),
        "arora_params": str(args.arora_params.resolve()),
        "arora_params_sha256": sha256_file(args.arora_params),
        "temperature_K": args.temperature,
        "fields_V_per_cm": args.fields,
        "fit_parameters": args.fit,
        "npart": args.npart,
        "time_s": args.time,
        "warmup": args.warmup,
        "max_energy_eV": args.max_energy,
        "gamma_safety": args.gamma_safety,
        "seed": args.seed,
        "objective_weights": {
            "velocity": args.velocity_weight,
            "mobility": args.mobility_weight,
            "regularization": args.regularization,
        },
    }
    with (args.output_dir / "run_manifest.json").open("w", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")

    objective = Objective(args, base_config, arora_config)
    x0 = [math.log(baseline[name]) for name in args.fit]
    initial_loss = objective(x0)
    if args.evaluate_only:
        print(f"base loss={initial_loss:.6g}", flush=True)
        return 0

    try:
        from scipy.optimize import minimize
    except ImportError as error:
        raise RuntimeError("SciPy is required for high-field optimization") from error
    result = minimize(
        objective,
        x0,
        method="Powell",
        bounds=parameter_bounds(args, baseline),
        options={"maxiter": args.maxiter, "xtol": args.xtol, "ftol": 1.0e-4, "disp": True},
    )
    print(f"optimizer success={result.success} message={result.message}", flush=True)
    print(f"best loss={objective.best_loss:.6g}", flush=True)
    print(f"best profile: {objective.best_profile}", flush=True)
    print(f"best rates: {objective.best_rates}", flush=True)
    print(f"history: {objective.history}", flush=True)
    return 0 if math.isfinite(objective.best_loss) else 1


if __name__ == "__main__":
    raise SystemExit(main())
