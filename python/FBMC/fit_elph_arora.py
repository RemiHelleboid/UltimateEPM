#!/usr/bin/env python3
"""Fit FBMC electron-phonon strengths to zero-doping Arora mobility.

The expensive, deformation-potential-independent kernel must already exist for
each fitted temperature. This script writes temporary electron-phonon YAML
profiles, runs elph.epm with each kernel, extracts the isotropic MRTA mobility,
and minimizes logarithmic mobility error.

The default fit adjusts only the zero-energy acoustic and optical strengths.
High-energy slopes and the threshold are deliberately left unchanged because
low-field mobility does not identify them reliably.
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
from dataclasses import dataclass
from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
MOBILITY_PATTERN = re.compile(
    r"(?:μ_iso|mu_iso)\s*=\s*([0-9.eE+-]+)\s*cm\^2/\(V[·*]s\)"
)
SUPPORTED_PARAMETERS = (
    "acoustic-a",
    "optical-a",
    "acoustic-end",
    "optical-end",
)


@dataclass(frozen=True)
class KernelInput:
    temperature_K: float
    path: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=REPO_ROOT / "build/apps/elph.epm")
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument(
        "--kernel",
        action="append",
        required=True,
        metavar="T=PATH",
        help="Kernel CSV for one temperature. Repeat, e.g. --kernel 300=electron_kernels_300K.csv.",
    )
    parser.add_argument(
        "--base-params",
        type=Path,
        default=REPO_ROOT / "data/materials/Si/electron_phonon/remi-2026.yaml",
    )
    parser.add_argument(
        "--arora-params",
        type=Path,
        default=REPO_ROOT / "data/materials/Si/admc/arora-canali.yaml",
    )
    parser.add_argument("--output-dir", type=Path, default=REPO_ROOT / "fit_elph_arora_output")
    parser.add_argument("--material", default="Si")
    parser.add_argument("--ncbands", type=int, default=2)
    parser.add_argument(
        "--nvbands",
        type=int,
        default=4,
        help="Valence bands loaded for the intrinsic Fermi-level calculation.",
    )
    parser.add_argument("--nthreads", type=int, default=1)
    parser.add_argument("--energy-window", type=float, default=0.3)
    parser.add_argument("--bz-domain", choices=("full", "octant"), default="full")
    parser.add_argument(
        "--fit",
        default="acoustic-a,optical-a",
        help="Comma-separated parameters: " + ",".join(SUPPORTED_PARAMETERS),
    )
    parser.add_argument(
        "--regularization",
        type=float,
        default=1.0e-3,
        help="Penalty on squared log-change from remi-2026.",
    )
    parser.add_argument("--maxiter", type=int, default=40)
    parser.add_argument("--xtol", type=float, default=2.0e-3)
    parser.add_argument(
        "--max-strength-factor",
        type=float,
        default=100.0,
        help="Maximum multiplicative change allowed for each fitted positive strength.",
    )
    parser.add_argument(
        "--evaluate-only",
        action="store_true",
        help="Evaluate the base profile once without optimizing.",
    )
    return parser.parse_args()


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    if not isinstance(config, dict):
        raise ValueError(f"Expected a YAML mapping in {path}")
    return config


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_kernel_header(path: Path) -> None:
    expected = [
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
    with path.open("r", encoding="utf-8", newline="") as stream:
        header = next(csv.reader(stream), None)
    if header != expected:
        raise ValueError(f"{path} does not have the expected electron-phonon kernel header")


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


def validate_kernel_metadata(
    kernel: KernelInput,
    args: argparse.Namespace,
    base_config: dict,
) -> dict | None:
    metadata_path = Path(str(kernel.path) + ".meta.yaml")
    if not metadata_path.is_file():
        print(
            f"warning: no metadata sidecar for {kernel.path}; temperature, dispersion, and mesh compatibility "
            "cannot be verified",
            flush=True,
        )
        return None
    metadata = load_yaml(metadata_path)
    checks = {
        "model": "electron_phonon_kernel",
        "material": args.material,
        "carrier": "electron",
        "n_conduction_bands": args.ncbands,
        "n_valence_bands": args.nvbands,
        "bz_domain": args.bz_domain,
    }
    for key, expected in checks.items():
        if metadata.get(key) != expected:
            raise ValueError(
                f"Kernel metadata mismatch for {kernel.path}: {key}={metadata.get(key)!r}, expected {expected!r}"
            )
    if not math.isclose(float(metadata["temperature_K"]), kernel.temperature_K, rel_tol=0.0, abs_tol=1.0e-9):
        raise ValueError(f"Kernel temperature metadata does not match --kernel label for {kernel.path}")
    kernel_window_eV = float(metadata["energy_window_eV"])
    if kernel_window_eV + 1.0e-12 < args.energy_window:
        raise ValueError(
            f"Kernel energy window {kernel_window_eV:g} eV is smaller than the requested "
            f"MRTA window {args.energy_window:g} eV for {kernel.path}"
        )
    if int(metadata["mesh_size_bytes"]) != args.mesh.stat().st_size:
        raise ValueError(f"Kernel mesh size metadata does not match --mesh for {kernel.path}")
    if normalize_numeric_tree(metadata.get("dispersion")) != normalize_numeric_tree(base_config.get("dispersion")):
        raise ValueError(f"Kernel phonon dispersion does not match --base-params for {kernel.path}")
    if not math.isclose(
        float(metadata["Radius-WS"]),
        float(base_config["Radius-WS"]),
        rel_tol=1.0e-12,
        abs_tol=0.0,
    ):
        raise ValueError(f"Kernel Radius-WS does not match --base-params for {kernel.path}")
    return metadata


def parse_kernels(values: list[str]) -> list[KernelInput]:
    kernels: list[KernelInput] = []
    temperatures: set[float] = set()
    for value in values:
        try:
            temperature_text, path_text = value.split("=", 1)
            temperature_K = float(temperature_text)
        except ValueError as error:
            raise ValueError(f"Invalid --kernel value {value!r}; expected T=PATH") from error
        path = Path(path_text).expanduser().resolve()
        if not math.isfinite(temperature_K) or temperature_K <= 0.0:
            raise ValueError(f"Invalid kernel temperature: {temperature_text}")
        if temperature_K in temperatures:
            raise ValueError(f"Duplicate kernel temperature: {temperature_K:g} K")
        if not path.is_file():
            raise FileNotFoundError(f"Kernel file not found: {path}")
        validate_kernel_header(path)
        temperatures.add(temperature_K)
        kernels.append(KernelInput(temperature_K, path))
    return sorted(kernels, key=lambda item: item.temperature_K)


def validate_args(args: argparse.Namespace) -> None:
    if not args.exe.is_file():
        raise FileNotFoundError(f"Executable not found: {args.exe}")
    if not os.access(args.exe, os.X_OK):
        raise PermissionError(f"Executable is not executable: {args.exe}")
    if not args.mesh.is_file():
        raise FileNotFoundError(f"Mesh file not found: {args.mesh}")
    if not args.base_params.is_file():
        raise FileNotFoundError(f"Base electron-phonon profile not found: {args.base_params}")
    if not args.arora_params.is_file():
        raise FileNotFoundError(f"Arora profile not found: {args.arora_params}")
    if args.ncbands <= 0 or args.nvbands <= 0:
        raise ValueError("This fitter requires conduction and valence bands for the intrinsic Fermi-level calculation")
    if args.nthreads <= 0:
        raise ValueError("--nthreads must be positive")
    if args.energy_window <= 0.0:
        raise ValueError("--energy-window must be positive")
    if args.regularization < 0.0:
        raise ValueError("--regularization must be non-negative")
    if args.max_strength_factor <= 1.0 or not math.isfinite(args.max_strength_factor):
        raise ValueError("--max-strength-factor must be finite and greater than one")

    args.fit = [name.strip() for name in args.fit.split(",") if name.strip()]
    if not args.fit:
        raise ValueError("--fit must contain at least one parameter")
    unknown = sorted(set(args.fit) - set(SUPPORTED_PARAMETERS))
    if unknown:
        raise ValueError(f"Unsupported fit parameters: {', '.join(unknown)}")
    if len(args.fit) != len(set(args.fit)):
        raise ValueError("--fit contains duplicate parameters")


def arora_zero_doping_mobility(config: dict, temperature_K: float) -> float:
    reference_temperature = float(config["reference_temperature_K"])
    low_field = config["electron"]["low_field"]
    relative_temperature = temperature_K / reference_temperature
    mu_min = float(low_field["mu_min_300_cm2_per_V_s"]) * relative_temperature ** float(
        low_field["mu_min_temperature_exponent"]
    )
    mu_dop = float(low_field["mu_dop_300_cm2_per_V_s"]) * relative_temperature ** float(
        low_field["mu_dop_temperature_exponent"]
    )
    return mu_min + mu_dop


def base_strengths(config: dict) -> dict[str, float]:
    electron = config["deformation-potential"]["electron"]
    threshold = float(electron["energy-threshold"])
    acoustic_a = float(electron["acoustic"]["A"])
    optical_a = float(electron["optic"]["A"])
    return {
        "acoustic-a": acoustic_a,
        "optical-a": optical_a,
        "acoustic-end": acoustic_a + float(electron["acoustic"]["B"]) * threshold,
        "optical-end": optical_a + float(electron["optic"]["B"]) * threshold,
    }


def write_trial_profile(
    base_config: dict,
    baseline: dict[str, float],
    fitted_names: list[str],
    log_values: list[float],
    path: Path,
) -> dict[str, float]:
    fitted = {
        name: math.exp(float(log_value))
        for name, log_value in zip(fitted_names, log_values, strict=True)
    }

    config = copy.deepcopy(base_config)
    config["parameter_set"] = "arora-fit-trial"
    electron = config["deformation-potential"]["electron"]
    threshold = float(electron["energy-threshold"])
    if threshold <= 0.0:
        raise ValueError("Electron energy-threshold must be positive")

    acoustic_a = fitted.get("acoustic-a", baseline["acoustic-a"])
    optical_a = fitted.get("optical-a", baseline["optical-a"])
    electron["acoustic"]["A"] = acoustic_a
    electron["optic"]["A"] = optical_a
    if "acoustic-end" in fitted:
        electron["acoustic"]["B"] = (fitted["acoustic-end"] - acoustic_a) / threshold
    if "optical-end" in fitted:
        electron["optic"]["B"] = (fitted["optical-end"] - optical_a) / threshold

    strengths = {
        "acoustic-a": acoustic_a,
        "optical-a": optical_a,
        "acoustic-end": acoustic_a + float(electron["acoustic"]["B"]) * threshold,
        "optical-end": optical_a + float(electron["optic"]["B"]) * threshold,
    }
    if any(not math.isfinite(value) or value <= 0.0 for value in strengths.values()):
        raise ValueError("Trial deformation-potential strength is non-positive or non-finite")

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(config, stream, sort_keys=False)
    return strengths


def run_elph(
    args: argparse.Namespace,
    kernel: KernelInput,
    profile: Path,
    run_dir: Path,
    evaluation: int,
) -> float:
    run_dir.mkdir(parents=True, exist_ok=True)
    command = [
        str(args.exe.resolve()),
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
        f"{kernel.temperature_K:.12g}",
        "--energy_window",
        f"{args.energy_window:.12g}",
        "--nthreads",
        str(args.nthreads),
        "--bz-domain",
        args.bz_domain,
        "--kernel-file",
        str(kernel.path),
        "--skip-mesh-vtk",
        "--outdir",
        str(run_dir.resolve()),
    ]
    completed = subprocess.run(
        command,
        cwd=run_dir,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    log_path = run_dir / f"elph_eval_{evaluation:04d}.log"
    log_path.write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"elph.epm failed at {kernel.temperature_K:g} K; see {log_path}"
        )
    match = MOBILITY_PATTERN.search(completed.stdout)
    if not match:
        raise RuntimeError(
            f"Could not extract isotropic mobility at {kernel.temperature_K:g} K; "
            f"see {log_path}"
        )
    mobility = float(match.group(1))
    if not math.isfinite(mobility) or mobility <= 0.0:
        raise RuntimeError(f"Invalid mobility returned at {kernel.temperature_K:g} K: {mobility}")
    return mobility


class Objective:
    def __init__(
        self,
        args: argparse.Namespace,
        kernels: list[KernelInput],
        base_config: dict,
        arora_config: dict,
    ) -> None:
        self.args = args
        self.kernels = kernels
        self.base_config = base_config
        self.baseline = base_strengths(base_config)
        self.arora_config = arora_config
        self.output_dir = args.output_dir.resolve()
        self.trial_profile = self.output_dir / "trial.yaml"
        self.history_path = self.output_dir / "history.csv"
        self.best_profile = self.output_dir / "best.yaml"
        self.best_loss = math.inf
        self.evaluation = 0
        self.cache: dict[tuple[float, ...], float] = {}
        self.trial_profile.unlink(missing_ok=True)
        self.best_profile.unlink(missing_ok=True)

        fieldnames = [
            "evaluation",
            "temperature_K",
            *SUPPORTED_PARAMETERS,
            "target_mobility_cm2_per_V_s",
            "model_mobility_cm2_per_V_s",
            "log_residual",
            "data_loss",
            "regularization_loss",
            "total_loss",
        ]
        with self.history_path.open("w", newline="", encoding="utf-8") as stream:
            csv.DictWriter(stream, fieldnames=fieldnames).writeheader()

    def __call__(self, log_values) -> float:
        values = [float(value) for value in log_values]
        cache_key = tuple(round(value, 12) for value in values)
        if cache_key in self.cache:
            return self.cache[cache_key]

        self.evaluation += 1
        try:
            strengths = write_trial_profile(
                self.base_config,
                self.baseline,
                self.args.fit,
                values,
                self.trial_profile,
            )
        except (OverflowError, ValueError) as error:
            print(
                f"evaluation {self.evaluation:04d} rejected before simulation: {error}",
                flush=True,
            )
            loss = 1.0e12
            self.cache[cache_key] = loss
            return loss

        rows = []
        residuals = []
        try:
            for kernel in self.kernels:
                target = arora_zero_doping_mobility(self.arora_config, kernel.temperature_K)
                model = run_elph(
                    self.args,
                    kernel,
                    self.trial_profile,
                    self.output_dir / "runs" / f"T_{kernel.temperature_K:g}K",
                    self.evaluation,
                )
                residual = math.log(model / target)
                residuals.append(residual)
                rows.append((kernel.temperature_K, target, model, residual))
        except (RuntimeError, ValueError) as error:
            print(f"evaluation {self.evaluation:04d} failed: {error}", flush=True)
            loss = 1.0e12
            self.cache[cache_key] = loss
            return loss

        data_loss = sum(value * value for value in residuals) / len(residuals)
        regularization = self.args.regularization * sum(
            (value - math.log(self.baseline[name])) ** 2
            for name, value in zip(self.args.fit, values, strict=True)
        )
        loss = data_loss + regularization

        with self.history_path.open("a", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=[
                "evaluation",
                "temperature_K",
                *SUPPORTED_PARAMETERS,
                "target_mobility_cm2_per_V_s",
                "model_mobility_cm2_per_V_s",
                "log_residual",
                "data_loss",
                "regularization_loss",
                "total_loss",
            ])
            for temperature, target, model, residual in rows:
                writer.writerow(
                    {
                        "evaluation": self.evaluation,
                        "temperature_K": f"{temperature:.12g}",
                        **{name: f"{strengths[name]:.12g}" for name in SUPPORTED_PARAMETERS},
                        "target_mobility_cm2_per_V_s": f"{target:.12g}",
                        "model_mobility_cm2_per_V_s": f"{model:.12g}",
                        "log_residual": f"{residual:.12g}",
                        "data_loss": f"{data_loss:.12g}",
                        "regularization_loss": f"{regularization:.12g}",
                        "total_loss": f"{loss:.12g}",
                    }
                )

        if loss < self.best_loss:
            self.best_loss = loss
            shutil.copy2(self.trial_profile, self.best_profile)

        parameter_text = " ".join(f"{name}={strengths[name]:.6g}" for name in self.args.fit)
        print(
            f"evaluation {self.evaluation:04d} loss={loss:.6g} {parameter_text}",
            flush=True,
        )
        self.cache[cache_key] = loss
        return loss


def main() -> int:
    args = parse_args()
    validate_args(args)
    kernels = parse_kernels(args.kernel)
    if len(kernels) < len(args.fit):
        raise ValueError(
            f"{len(args.fit)} fitted parameters require at least {len(args.fit)} temperature kernels; "
            f"received {len(kernels)}"
        )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    base_config = load_yaml(args.base_params)
    arora_config = load_yaml(args.arora_params)
    if base_config.get("material") != args.material:
        raise ValueError("Base electron-phonon profile material does not match --material")
    if base_config.get("model") != "electron_phonon":
        raise ValueError("Base profile is not an electron_phonon model")
    if arora_config.get("material") != args.material:
        raise ValueError("Arora profile material does not match --material")
    if arora_config.get("model") != "arora_canali_mobility":
        raise ValueError("Arora profile has an unexpected model identifier")
    kernel_metadata = [validate_kernel_metadata(kernel, args, base_config) for kernel in kernels]

    manifest = {
        "executable": str(args.exe.resolve()),
        "mesh": str(args.mesh.resolve()),
        "mesh_sha256": sha256_file(args.mesh),
        "base_params": str(args.base_params.resolve()),
        "base_params_sha256": sha256_file(args.base_params),
        "arora_params": str(args.arora_params.resolve()),
        "arora_params_sha256": sha256_file(args.arora_params),
        "material": args.material,
        "ncbands": args.ncbands,
        "nvbands": args.nvbands,
        "energy_window_eV": args.energy_window,
        "bz_domain": args.bz_domain,
        "fit_parameters": args.fit,
        "regularization": args.regularization,
        "kernels": [
            {
                "temperature_K": kernel.temperature_K,
                "path": str(kernel.path),
                "sha256": sha256_file(kernel.path),
                "metadata": metadata,
            }
            for kernel, metadata in zip(kernels, kernel_metadata, strict=True)
        ],
    }
    with (args.output_dir / "run_manifest.json").open("w", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")

    baseline = base_strengths(base_config)
    for name, value in baseline.items():
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"Base strength {name} must be positive; got {value}")
    x0 = [math.log(baseline[name]) for name in args.fit]

    objective = Objective(args, kernels, base_config, arora_config)
    initial_loss = objective(x0)
    if args.evaluate_only:
        print(f"base loss={initial_loss:.6g}", flush=True)
        return 0

    try:
        from scipy.optimize import minimize
    except ImportError as error:
        raise RuntimeError("SciPy is required for optimization; use --evaluate-only without it") from error

    print(f"Optionally fitting parameters: {', '.join(args.fit)}", flush=True)
    result = minimize(
        objective,
        x0,
        method="Powell",
        bounds=[
            (
                value - math.log(args.max_strength_factor),
                value + math.log(args.max_strength_factor),
            )
            for value in x0
        ],
        options={"maxiter": args.maxiter, "xtol": args.xtol, "ftol": 1.0e-5, "disp": True},
    )

    print(f"optimizer success={result.success} message={result.message}", flush=True)
    print(f"best loss={objective.best_loss:.6g}", flush=True)
    print(f"best profile: {objective.best_profile}", flush=True)
    print(f"history: {objective.history_path}", flush=True)
    return 0 if math.isfinite(objective.best_loss) else 1


if __name__ == "__main__":
    raise SystemExit(main())
