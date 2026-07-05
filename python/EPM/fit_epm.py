#!/usr/bin/env python3
"""Fit local EPM form factors by driving the C++ EPM analysis apps.

The script edits a generated EPM YAML parameter set, runs:
  - epm_band_edges for gap and valley-position targets
  - epm_valley_fit for effective masses and Kane non-parabolicity
  - optionally epsilon.epm for the optical dielectric function

Then it minimizes a weighted normalized least-squares score.
"""

from __future__ import annotations

import argparse
import copy
import csv
import math
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_EPSILON_REFERENCE = REPO_ROOT / "examples" / "references" / "DielectricFunction_Si_300K.csv"


@dataclass(frozen=True)
class Target:
    value: float
    scale: float
    weight: float = 1.0


@dataclass(frozen=True)
class FitParameter:
    label: str
    section: str
    keys: tuple[str, ...]


NONLOCAL_TIED_ALIASES: dict[str, tuple[str, str]] = {
    "alpha_0": ("alpha_0_cation", "alpha_0_anion"),
    "beta_0": ("beta_0_cation", "beta_0_anion"),
    "A2": ("A2_cation", "A2_anion"),
    "R0": ("R0_cation", "R0_anion"),
    "R2": ("R2_cation", "R2_anion"),
}
NONLOCAL_TIED_ALIASES_LOWER: dict[str, tuple[str, str]] = {
    key.lower(): value for key, value in NONLOCAL_TIED_ALIASES.items()
}


DEFAULT_TARGETS: dict[str, Target] = {
    "indirect_gap_eV": Target(1.12, 0.02, 4.0),
    "gamma_gap_eV": Target(3.40, 0.05, 1.0),
    "delta_k": Target(0.85, 0.01, 2.0),
    "x_cbm_rel_vbm_eV": Target(1.20, 0.05, 0.5),
    "l_cbm_rel_vbm_eV": Target(2.10, 0.08, 0.5),
    "mt_m0": Target(0.1905, 0.01, 2.0),
    "ml_m0": Target(0.9163, 0.03, 2.0),
    "alpha_eV_inv": Target(0.5, 0.08, 1.0),
}

PENALTY_SCORE = 1.0e12


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--material", default="Si")
    parser.add_argument("--base-set", default="local-remi-2026", help="Named repository EPM parameter set.")
    parser.add_argument(
        "--base-file",
        default="",
        help="External EPM YAML file to use as the starting parameter set. Overrides --base-set.",
    )
    parser.add_argument("--work-set", default="local-fit-working")
    parser.add_argument("--build-dir", default=str(REPO_ROOT / "build"))
    parser.add_argument(
        "--output-dir",
        default="fit_epm_output",
        help="Directory for fit logs and spectra. Relative paths are resolved from the launch directory.",
    )
    parser.add_argument("--maxiter", type=int, default=60)
    parser.add_argument("--method", choices=["auto", "powell", "nelder-mead", "coordinate"], default="auto")
    parser.add_argument(
        "--fit-params",
        default="pseudo.V3S,pseudo.V8S,pseudo.V11S",
        help=(
            "Comma-separated fitted parameters. Examples: "
            "pseudo.V3S,pseudo.V8S,pseudo.V11S,nonlocal.alpha_0,nonlocal.beta_0,nonlocal.R0"
        ),
    )
    parser.add_argument(
        "--param-bounds",
        default="",
        help=(
            "Comma-separated bound overrides as name:min:max, e.g. "
            "pseudo.V3S:-0.5:0,nonlocal.alpha_0:0:2,nonlocal.beta_0:0:2"
        ),
    )
    parser.add_argument(
        "--no-default-bounds",
        action="store_true",
        help="Disable built-in parameter bounds. Explicit --param-bounds still apply.",
    )
    parser.add_argument(
        "--nonlocal",
        dest="enable_nonlocal",
        action="store_true",
        help="Enable nonlocal EPM corrections in epm_band_edges, epm_valley_fit, and epsilon.epm.",
    )
    parser.add_argument(
        "--metrics",
        choices=["bands", "epsilon", "both"],
        default="bands",
        help="Which objective terms to activate.",
    )
    parser.add_argument("--nthreads", type=int, default=2)
    parser.add_argument("--nearest-neighbors", type=int, default=10)
    parser.add_argument("--delta-samples", type=int, default=401)
    parser.add_argument("--band", type=int, default=4)
    parser.add_argument("--k0", default="0.85,0,0")
    parser.add_argument("--mass-radius", type=float, default=0.01)
    parser.add_argument("--mass-shells", type=int, default=10)
    parser.add_argument("--alpha-radius", type=float, default=0.05)
    parser.add_argument("--alpha-shells", type=int, default=30)
    parser.add_argument("--alpha-max-energy", type=float, default=1.0e-2)
    parser.add_argument("--epsilon-reference", default=str(DEFAULT_EPSILON_REFERENCE), help=argparse.SUPPRESS)
    parser.add_argument(
        "--epsilon-component",
        choices=["real", "imag", "both"],
        default="both",
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--epsilon-weight", type=float, default=1.0)
    parser.add_argument("--epsilon-static-target", type=float, default=11.7, help="Target epsilon1(E=0).")
    parser.add_argument("--epsilon-argmax-target", type=float, default=3.3, help="Target energy in eV of max epsilon1.")
    parser.add_argument("--epsilon-argmin-target", type=float, default=4.5, help="Target energy in eV of min epsilon1.")
    parser.add_argument("--epsilon-static-scale", type=float, default=1.0, help="Normalization scale for epsilon1(E=0).")
    parser.add_argument("--epsilon-extrema-scale", type=float, default=0.2, help="Normalization scale in eV for extrema positions.")
    parser.add_argument("--epsilon-real-weight", type=float, default=1.0, help=argparse.SUPPRESS)
    parser.add_argument("--epsilon-imag-weight", type=float, default=1.0, help=argparse.SUPPRESS)
    parser.add_argument("--epsilon-real-scale", type=float, default=10.0, help=argparse.SUPPRESS)
    parser.add_argument("--epsilon-imag-scale", type=float, default=10.0, help=argparse.SUPPRESS)
    parser.add_argument("--epsilon-bands", type=int, default=16)
    parser.add_argument(
        "--epsilon-bz-sampling",
        choices=["full", "q100-octant", "fcc-ibz", "1", "8", "48"],
        default="full",
        help="BZ sampling mode passed to epsilon.epm.",
    )
    parser.add_argument("--epsilon-nkx", type=int, default=8)
    parser.add_argument("--epsilon-nky", type=int, default=8)
    parser.add_argument("--epsilon-nkz", type=int, default=8)
    parser.add_argument("--epsilon-q", default="1e-3")
    parser.add_argument("--epsilon-direction", default="1,0,0")
    parser.add_argument("--epsilon-emin", type=float, default=0.0)
    parser.add_argument("--epsilon-emax", type=float, default=6.0)
    parser.add_argument("--epsilon-estep", type=float, default=0.1)
    parser.add_argument("--epsilon-eta", type=float, default=0.15)
    parser.add_argument("--epsilon-nearest-neighbors", type=int, default=10)
    parser.add_argument(
        "--epsilon-mpi-ranks",
        type=int,
        default=1,
        help="Number of MPI ranks for epsilon.epm. Use 1 for serial.",
    )
    parser.add_argument(
        "--epsilon-mpi-runner",
        default="auto",
        help="MPI launcher for epsilon.epm. 'auto' prefers /usr/bin/mpirun when available.",
    )
    parser.add_argument(
        "--epsilon-mpi-extra-args",
        default="",
        help="Extra MPI launcher arguments, parsed like a shell string.",
    )
    parser.add_argument("--initial-step", type=float, default=0.01, help="Initial coordinate-search step in Rydberg.")
    parser.add_argument("--min-step", type=float, default=2.0e-4, help="Minimum coordinate-search step in Rydberg.")
    parser.add_argument("--dry-run", action="store_true", help="Evaluate the initial parameter set once and exit.")
    return parser.parse_args()


def uses_band_metric(args: argparse.Namespace) -> bool:
    return args.metrics in ("bands", "both")


def uses_epsilon_metric(args: argparse.Namespace) -> bool:
    return args.metrics in ("epsilon", "both")


def parameter_file(material: str, parameter_set: str) -> Path:
    return REPO_ROOT / "data" / "materials" / material / "epm" / f"{parameter_set}.yaml"


def resolve_from_launch_dir(path_like: str, launch_dir: Path) -> Path:
    path = Path(path_like).expanduser()
    if path.is_absolute():
        return path
    return launch_dir / path


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as stream:
        return yaml.safe_load(stream)


def validate_epm_config(config: dict, path: Path, material: str) -> None:
    if config.get("material") != material:
        raise ValueError(f"{path} has material={config.get('material')!r}, expected {material!r}")
    if config.get("model") != "epm":
        raise ValueError(f"{path} has model={config.get('model')!r}, expected 'epm'")
    if "pseudo-potential-parameters" not in config:
        raise ValueError(f"{path} does not define pseudo-potential-parameters")


def parse_fit_parameter(token: str) -> FitParameter:
    raw = token.strip()
    if not raw:
        raise ValueError("empty fit parameter")
    if "." in raw:
        namespace, key = raw.split(".", 1)
    else:
        namespace, key = "pseudo", raw
    namespace = namespace.strip().lower()
    key = key.strip()
    if not key:
        raise ValueError(f"invalid fit parameter '{token}'")

    if namespace in ("pseudo", "pseudopotential", "local"):
        return FitParameter(f"pseudo.{key}", "pseudo-potential-parameters", (key,))
    if namespace in ("nonlocal", "non-local"):
        keys = NONLOCAL_TIED_ALIASES.get(key, NONLOCAL_TIED_ALIASES_LOWER.get(key.lower(), (key,)))
        return FitParameter(f"nonlocal.{key}", "non-local-parameters", tuple(keys))
    if namespace in ("soc", "spinorbit", "spin-orbit"):
        return FitParameter(f"spinorbit.{key}", "spin-orbit-parameters", (key,))
    raise ValueError(f"unknown fit parameter namespace '{namespace}' in '{token}'")


def parse_fit_parameters(value: str) -> list[FitParameter]:
    params = [parse_fit_parameter(token) for token in value.split(",") if token.strip()]
    if not params:
        raise ValueError("--fit-params must contain at least one parameter")
    labels = [param.label for param in params]
    if len(labels) != len(set(labels)):
        raise ValueError(f"--fit-params contains duplicate entries: {value}")
    return params


def get_config_section(config: dict, section: str) -> dict:
    if section not in config or config[section] is None:
        raise KeyError(f"parameter file does not define '{section}'")
    return config[section]


def initial_parameter_values(config: dict, fit_params: list[FitParameter]) -> list[float]:
    values: list[float] = []
    for param in fit_params:
        section = get_config_section(config, param.section)
        missing = [key for key in param.keys if key not in section]
        if missing:
            raise KeyError(f"{param.label} references missing YAML key(s): {', '.join(missing)}")
        tied_values = [float(section[key]) for key in param.keys]
        if max(tied_values) - min(tied_values) > 1.0e-12:
            print(
                f"warning: tied parameter {param.label} has unequal initial values {tied_values}; "
                f"using {tied_values[0]} and tying them during the fit",
                flush=True,
            )
        values.append(tied_values[0])
    return values


def default_bound_for_param(param: FitParameter, initial_value: float) -> tuple[float, float]:
    if param.section == "pseudo-potential-parameters":
        return (-1.0, 1.0)
    if param.section == "non-local-parameters":
        key_lower = param.label.lower()
        if ".r0" in key_lower or ".r2" in key_lower:
            return (0.0, 5.0)
        if ".alpha_0" in key_lower or ".beta_0" in key_lower:
            return (0.0, 2.0)
        if ".a2" in key_lower:
            return (-2.0, 2.0)
        return (-5.0, 5.0)
    span = max(1.0, 5.0 * abs(initial_value))
    return (initial_value - span, initial_value + span)


def parse_bound_overrides(value: str) -> dict[str, tuple[float, float]]:
    overrides: dict[str, tuple[float, float]] = {}
    if not value.strip():
        return overrides
    for token in value.split(","):
        token = token.strip()
        if not token:
            continue
        parts = token.split(":")
        if len(parts) != 3:
            raise ValueError(f"invalid bound override '{token}', expected name:min:max")
        label = parse_fit_parameter(parts[0]).label
        lower = float(parts[1])
        upper = float(parts[2])
        if not (math.isfinite(lower) and math.isfinite(upper) and lower < upper):
            raise ValueError(f"invalid finite bounds for {label}: {lower}, {upper}")
        overrides[label] = (lower, upper)
    return overrides


def parameter_bounds(args: argparse.Namespace, fit_params: list[FitParameter], x0: list[float]) -> list[tuple[float, float]]:
    overrides = parse_bound_overrides(args.param_bounds)
    known_labels = {param.label for param in fit_params}
    unknown = sorted(set(overrides) - known_labels)
    if unknown:
        raise ValueError(f"--param-bounds references parameter(s) not in --fit-params: {', '.join(unknown)}")

    bounds: list[tuple[float, float]] = []
    for param, initial_value in zip(fit_params, x0):
        lower, upper = (-math.inf, math.inf) if args.no_default_bounds else default_bound_for_param(param, initial_value)
        if param.label in overrides:
            lower, upper = overrides[param.label]
        if not (lower <= initial_value <= upper):
            raise ValueError(f"initial value {initial_value} for {param.label} is outside bounds [{lower}, {upper}]")
        bounds.append((lower, upper))
    return bounds


def finite_params(params: list[float]) -> bool:
    return all(math.isfinite(value) for value in params)


def within_bounds(params: list[float], bounds: list[tuple[float, float]]) -> bool:
    return all(lower <= value <= upper for value, (lower, upper) in zip(params, bounds))


def clip_to_bounds(params: list[float], bounds: list[tuple[float, float]]) -> list[float]:
    clipped: list[float] = []
    for value, (lower, upper) in zip(params, bounds):
        if math.isfinite(lower):
            value = max(value, lower)
        if math.isfinite(upper):
            value = min(value, upper)
        clipped.append(value)
    return clipped


def write_working_yaml(
    base_config: dict,
    material: str,
    work_set: str,
    fit_params: list[FitParameter],
    params: list[float],
) -> Path:
    if len(fit_params) != len(params):
        raise ValueError(f"got {len(params)} values for {len(fit_params)} fit parameters")
    config = copy.deepcopy(base_config)
    config["parameter_set"] = work_set
    for spec, value in zip(fit_params, params):
        section = get_config_section(config, spec.section)
        for key in spec.keys:
            if key not in section:
                raise KeyError(f"{spec.label} references missing YAML key '{key}'")
            section[key] = float(value)

    path = parameter_file(material, work_set)
    with path.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(config, stream, sort_keys=False)
    return path


def format_params(fit_params: list[FitParameter], params: list[float]) -> str:
    return " ".join(f"{spec.label}={value:.7g}" for spec, value in zip(fit_params, params))


def read_quantity_csv(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as stream:
        reader = csv.reader(stream)
        header = next(reader, None)
        if header != ["quantity", "value"]:
            raise RuntimeError(f"{path} does not look like a quantity,value CSV")
        for row in reader:
            if len(row) >= 2:
                values[row[0]] = row[1]
    return values


def read_epsilon_csv(path: Path, energy_key: str, real_key: str, imag_key: str) -> tuple[list[float], list[float], list[float]]:
    energies: list[float] = []
    real: list[float] = []
    imag: list[float] = []
    with path.open("r", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        for row in reader:
            energies.append(float(row[energy_key]))
            real.append(float(row[real_key]))
            imag.append(float(row[imag_key]))
    if len(energies) < 2:
        raise RuntimeError(f"{path} does not contain enough dielectric samples")
    return energies, real, imag


def interpolate_linear(xs: list[float], ys: list[float], x: float) -> float:
    if x < xs[0] or x > xs[-1]:
        raise ValueError(f"{x} outside interpolation range [{xs[0]}, {xs[-1]}]")
    lo = 0
    hi = len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    if xs[lo] == x or hi == lo:
        return ys[lo]
    dx = xs[hi] - xs[lo]
    if dx <= 0.0:
        raise RuntimeError("dielectric energy grid must be strictly increasing")
    t = (x - xs[lo]) / dx
    return (1.0 - t) * ys[lo] + t * ys[hi]


def extremum_energy(energies: list[float], values: list[float], find_max: bool) -> float:
    if len(energies) != len(values) or len(energies) < 3:
        raise RuntimeError("need at least three epsilon samples to locate an extremum")
    index = max(range(len(values)), key=values.__getitem__) if find_max else min(range(len(values)), key=values.__getitem__)
    if index == 0 or index == len(values) - 1:
        return energies[index]

    x0, x1, x2 = energies[index - 1], energies[index], energies[index + 1]
    y0, y1, y2 = values[index - 1], values[index], values[index + 1]
    if not (math.isfinite(x0) and math.isfinite(x1) and math.isfinite(x2)):
        return energies[index]
    if abs((x1 - x0) - (x2 - x1)) > 1.0e-8:
        return energies[index]

    denominator = y0 - 2.0 * y1 + y2
    if abs(denominator) < 1.0e-14:
        return energies[index]
    step = x1 - x0
    x_vertex = x1 + 0.5 * step * (y0 - y2) / denominator
    if x0 <= x_vertex <= x2:
        return x_vertex
    return energies[index]


def score_epsilon_features(
    args: argparse.Namespace,
    model_csv: Path,
) -> dict[str, float]:
    model_energy, model_real, model_imag = read_epsilon_csv(
        model_csv,
        "Energy (eV)",
        "EpsilonReal",
        "EpsilonImaginary",
    )
    del model_imag
    if any(not math.isfinite(value) for value in model_energy + model_real):
        raise RuntimeError(f"{model_csv} contains non-finite epsilon samples")

    epsilon0 = interpolate_linear(model_energy, model_real, 0.0)
    argmax_eV = extremum_energy(model_energy, model_real, find_max=True)
    argmin_eV = extremum_energy(model_energy, model_real, find_max=False)

    static_residual = (epsilon0 - args.epsilon_static_target) / args.epsilon_static_scale
    argmax_residual = (argmax_eV - args.epsilon_argmax_target) / args.epsilon_extrema_scale
    argmin_residual = (argmin_eV - args.epsilon_argmin_target) / args.epsilon_extrema_scale
    score = args.epsilon_weight * (static_residual * static_residual +
                                   argmax_residual * argmax_residual +
                                   argmin_residual * argmin_residual)

    return {
        "epsilon_score": score,
        "epsilon0": epsilon0,
        "epsilon_argmax_eV": argmax_eV,
        "epsilon_argmin_eV": argmin_eV,
    }


def as_float(values: dict[str, str], key: str) -> float:
    return float(values[key])


def run_command(cmd: list[str], cwd: Path) -> None:
    result = subprocess.run(cmd, cwd=cwd, check=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if result.returncode != 0:
        command = " ".join(cmd)
        raise RuntimeError(f"command failed with exit code {result.returncode}: {command}\n{result.stdout}")


def resolve_mpi_runner(runner: str) -> str:
    if runner != "auto":
        return runner
    system_mpirun = Path("/usr/bin/mpirun")
    if system_mpirun.exists():
        return str(system_mpirun)
    return "mpirun"


def epsilon_command(args: argparse.Namespace, work_set: str, epsilon_prefix: Path) -> list[str]:
    build_apps = Path(args.build_dir) / "apps"
    command = [
        str(build_apps / "epsilon.epm"),
        "--material",
        args.material,
        "--epm-set",
        work_set,
        "--bands",
        str(args.epsilon_bands),
        "--nearest-neighbors",
        str(args.epsilon_nearest_neighbors),
        "--Nkx",
        str(args.epsilon_nkx),
        "--Nky",
        str(args.epsilon_nky),
        "--Nkz",
        str(args.epsilon_nkz),
        "--mode",
        "q-list",
        "--q-values",
        args.epsilon_q,
        "--direction",
        args.epsilon_direction,
        "--emin",
        str(args.epsilon_emin),
        "--emax",
        str(args.epsilon_emax),
        "--estep",
        str(args.epsilon_estep),
        "--eta",
        str(args.epsilon_eta),
        "--bz-sampling",
        args.epsilon_bz_sampling,
        "--out",
        str(epsilon_prefix),
    ]
    if args.enable_nonlocal:
        command.append("--nonlocal")
    if args.epsilon_mpi_ranks <= 1:
        return command

    mpi_command = [resolve_mpi_runner(args.epsilon_mpi_runner), "-np", str(args.epsilon_mpi_ranks)]
    if args.epsilon_mpi_extra_args:
        mpi_command.extend(shlex.split(args.epsilon_mpi_extra_args))
    return mpi_command + command


def measure_bands(
    args: argparse.Namespace,
    work_set: str,
    output_dir: Path,
    iteration: int,
) -> dict[str, float]:
    edges_csv = output_dir / f"edges_{iteration:04d}.csv"
    valley_csv = output_dir / f"valley_{iteration:04d}.csv"
    build_apps = Path(args.build_dir) / "apps"

    edges_command = [
        str(build_apps / "epm_band_edges"),
        "--material",
        args.material,
        "--epm-set",
        work_set,
        "--nthreads",
        str(args.nthreads),
        "--nearestNeighbors",
        str(args.nearest_neighbors),
        "--delta-samples",
        str(args.delta_samples),
        "--out",
        str(edges_csv),
    ]
    if args.enable_nonlocal:
        edges_command.append("--nonlocal-correction")
    run_command(edges_command, REPO_ROOT)

    edges = read_quantity_csv(edges_csv)
    k_delta = as_float(edges, "delta_cbm_kx_reduced")

    valley_command = [
        str(build_apps / "epm_valley_fit"),
        "--material",
        args.material,
        "--epm-set",
        work_set,
        "--band",
        str(args.band),
        "--k0",
        f"{k_delta},0,0",
        "--mass-radius",
        str(args.mass_radius),
        "--mass-shells",
        str(args.mass_shells),
        "--radius",
        str(args.alpha_radius),
        "--shells",
        str(args.alpha_shells),
        "--alpha-max-energy",
        str(args.alpha_max_energy),
        "--nthreads",
        str(args.nthreads),
        "--nearestNeighbors",
        str(args.nearest_neighbors),
        "--out",
        str(valley_csv),
    ]
    if args.enable_nonlocal:
        valley_command.append("--nonlocal-correction")
    run_command(valley_command, REPO_ROOT)

    valley = read_quantity_csv(valley_csv)
    masses = sorted([as_float(valley, "m1_m0"), as_float(valley, "m2_m0"), as_float(valley, "m3_m0")])
    return {
        "indirect_gap_eV": as_float(edges, "indirect_gap_eV"),
        "gamma_gap_eV": as_float(edges, "gamma_gap_eV"),
        "delta_k": k_delta,
        "x_cbm_rel_vbm_eV": as_float(edges, "x_cbm_rel_vbm_eV"),
        "l_cbm_rel_vbm_eV": as_float(edges, "l_cbm_rel_vbm_eV"),
        "mt_m0": 0.5 * (masses[0] + masses[1]),
        "ml_m0": masses[2],
        "alpha_eV_inv": as_float(valley, "non_parabolicity_eV_inv"),
        "mass_rms_error_meV": as_float(valley, "mass_rms_error_meV"),
        "alpha_rms_error_meV": as_float(valley, "alpha_rms_error_meV"),
    }


def measure_epsilon(
    args: argparse.Namespace,
    work_set: str,
    output_dir: Path,
    iteration: int,
) -> dict[str, float]:
    epsilon_prefix = output_dir / f"epsilon_{iteration:04d}"
    run_command(epsilon_command(args, work_set, epsilon_prefix), REPO_ROOT)

    spectra = sorted(
        path
        for path in output_dir.glob(f"epsilon_{iteration:04d}_*.csv")
        if not path.name.endswith("_kpoints.csv")
    )
    if not spectra:
        raise RuntimeError(f"epsilon.epm did not write a spectrum for prefix {epsilon_prefix}")
    return score_epsilon_features(args, spectra[0])


def measure(
    args: argparse.Namespace,
    work_set: str,
    output_dir: Path,
    iteration: int,
) -> dict[str, float]:
    observables: dict[str, float] = {}
    if uses_band_metric(args):
        observables.update(measure_bands(args, work_set, output_dir, iteration))
    if uses_epsilon_metric(args):
        observables.update(measure_epsilon(args, work_set, output_dir, iteration))
    return observables


def score_observables(args: argparse.Namespace, observables: dict[str, float]) -> float:
    score = 0.0
    if uses_band_metric(args):
        for key, target in DEFAULT_TARGETS.items():
            residual = (observables[key] - target.value) / target.scale
            score += target.weight * residual * residual
        if observables["alpha_rms_error_meV"] > 5.0:
            score += ((observables["alpha_rms_error_meV"] - 5.0) / 5.0) ** 2
    if uses_epsilon_metric(args):
        score += observables["epsilon_score"]
    return score


def observable_columns(args: argparse.Namespace) -> list[str]:
    columns: list[str] = []
    if uses_band_metric(args):
        columns.extend(DEFAULT_TARGETS.keys())
        columns.extend(["mass_rms_error_meV", "alpha_rms_error_meV"])
    if uses_epsilon_metric(args):
        columns.extend(["epsilon_score", "epsilon0", "epsilon_argmax_eV", "epsilon_argmin_eV"])
    return columns


class Objective:
    def __init__(
        self,
        args: argparse.Namespace,
        base_config: dict,
        output_dir: Path,
        fit_params: list[FitParameter],
        bounds: list[tuple[float, float]],
    ) -> None:
        self.args = args
        self.base_config = base_config
        self.output_dir = output_dir
        self.fit_params = fit_params
        self.bounds = bounds
        self.iteration = 0
        self.best_score = math.inf
        self.best_params: list[float] | None = None
        self.log_path = output_dir / "iterations.csv"
        self.observable_columns = observable_columns(args)
        with self.log_path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                [
                    "iteration",
                    "score",
                    *(param.label for param in self.fit_params),
                    *self.observable_columns,
                ]
            )

    def __call__(self, params_like) -> float:
        params = [float(x) for x in params_like]
        iteration = self.iteration
        self.iteration += 1
        if not finite_params(params):
            print(f"iteration {iteration:04d} rejected non-finite params: {params}", flush=True)
            observables = {key: math.nan for key in self.observable_columns}
            self.write_log_row(iteration, PENALTY_SCORE, params, observables)
            return PENALTY_SCORE
        if not within_bounds(params, self.bounds):
            print(f"iteration {iteration:04d} rejected out-of-bounds params: {format_params(self.fit_params, params)}", flush=True)
            observables = {key: math.nan for key in self.observable_columns}
            self.write_log_row(iteration, PENALTY_SCORE, params, observables)
            return PENALTY_SCORE
        write_working_yaml(self.base_config, self.args.material, self.args.work_set, self.fit_params, params)
        try:
            observables = measure(self.args, self.args.work_set, self.output_dir, iteration)
            score = score_observables(self.args, observables)
            if not math.isfinite(score):
                raise RuntimeError(f"non-finite objective score: {score}")
        except Exception as exc:
            print(f"iteration {iteration:04d} failed for {params}: {exc}", flush=True)
            score = PENALTY_SCORE
            observables = {key: math.nan for key in self.observable_columns}

        self.write_log_row(iteration, score, params, observables)

        if score < self.best_score and score < PENALTY_SCORE:
            self.best_score = score
            self.best_params = params
            shutil.copy2(parameter_file(self.args.material, self.args.work_set), self.output_dir / "best.yaml")

        print(
            f"iter {iteration:04d} score={score:.6g} {format_params(self.fit_params, params)}",
            flush=True,
        )
        return score

    def write_log_row(
        self,
        iteration: int,
        score: float,
        params: list[float],
        observables: dict[str, float],
    ) -> None:
        with self.log_path.open("a", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                [
                    iteration,
                    score,
                    *params,
                    *(observables[key] for key in self.observable_columns),
                ]
            )


def coordinate_search(
    objective: Callable[[list[float]], float],
    x0: list[float],
    maxiter: int,
    step0: float,
    min_step: float,
    bounds: list[tuple[float, float]],
):
    x = clip_to_bounds(list(x0), bounds)
    best = objective(x)
    step = step0
    iterations = 1
    while iterations < maxiter and step >= min_step:
        improved = False
        for dim in range(len(x)):
            for direction in (1.0, -1.0):
                trial = list(x)
                trial[dim] += direction * step
                trial = clip_to_bounds(trial, bounds)
                value = objective(trial)
                iterations += 1
                if value < best:
                    x = trial
                    best = value
                    improved = True
                if iterations >= maxiter:
                    break
            if iterations >= maxiter:
                break
        if not improved:
            step *= 0.5
    return x, best


def run_scipy(objective: Objective, x0: list[float], args: argparse.Namespace) -> None:
    try:
        from scipy.optimize import Bounds, minimize
    except ImportError:
        print("SciPy not found; falling back to coordinate search.", flush=True)
        coordinate_search(objective, x0, args.maxiter, args.initial_step, args.min_step, objective.bounds)
        return

    method = "Powell" if args.method in ("auto", "powell") else "Nelder-Mead"
    options = {"maxiter": args.maxiter, "disp": True}
    if method == "Powell":
        options["xtol"] = args.min_step
        options["ftol"] = 1.0e-4
    else:
        options["xatol"] = args.min_step
        options["fatol"] = 1.0e-4
    lower = [bound[0] for bound in objective.bounds]
    upper = [bound[1] for bound in objective.bounds]
    minimize(objective, x0, method=method, bounds=Bounds(lower, upper), options=options)


def main() -> int:
    launch_dir = Path.cwd().resolve()
    args = parse_args()
    if args.nthreads <= 0:
        raise SystemExit("--nthreads must be positive")
    if uses_epsilon_metric(args):
        if args.epsilon_bands <= 0:
            raise SystemExit("--epsilon-bands must be positive")
        if args.epsilon_mpi_ranks <= 0:
            raise SystemExit("--epsilon-mpi-ranks must be positive")
        if args.epsilon_nkx <= 0 or args.epsilon_nky <= 0 or args.epsilon_nkz <= 0:
            raise SystemExit("--epsilon-nkx/--epsilon-nky/--epsilon-nkz must be positive")
        if args.epsilon_estep <= 0.0:
            raise SystemExit("--epsilon-estep must be positive")
        if args.epsilon_eta <= 0.0:
            raise SystemExit("--epsilon-eta must be positive")
        if args.epsilon_static_scale <= 0.0 or args.epsilon_extrema_scale <= 0.0:
            raise SystemExit("--epsilon-static-scale/--epsilon-extrema-scale must be positive")
        for name in ("epsilon_static_target", "epsilon_argmax_target", "epsilon_argmin_target"):
            if not math.isfinite(getattr(args, name)):
                raise SystemExit(f"--{name.replace('_', '-')} must be finite")

    base_path = resolve_from_launch_dir(args.base_file, launch_dir) if args.base_file else parameter_file(args.material, args.base_set)
    if not base_path.exists():
        raise SystemExit(f"base parameter file does not exist: {base_path}")
    base_config = load_yaml(base_path)
    try:
        validate_epm_config(base_config, base_path, args.material)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    try:
        fit_params = parse_fit_parameters(args.fit_params)
        if any(param.section == "non-local-parameters" for param in fit_params) and not args.enable_nonlocal:
            raise ValueError("--fit-params includes nonlocal.* entries; add --nonlocal so they affect the Hamiltonian")
        if args.enable_nonlocal:
            get_config_section(base_config, "non-local-parameters")
        x0 = initial_parameter_values(base_config, fit_params)
        bounds = parameter_bounds(args, fit_params, x0)
    except (KeyError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc

    output_dir = resolve_from_launch_dir(args.output_dir, launch_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"output dir: {output_dir}", flush=True)
    print(f"base file: {base_path}", flush=True)
    print(f"working set: {args.work_set}", flush=True)
    print(f"fit params: {', '.join(param.label for param in fit_params)}", flush=True)
    print(
        "bounds: "
        + ", ".join(f"{param.label}=[{lower:g},{upper:g}]" for param, (lower, upper) in zip(fit_params, bounds)),
        flush=True,
    )
    print(f"nonlocal corrections: {int(args.enable_nonlocal)}", flush=True)

    objective = Objective(args, base_config, output_dir, fit_params, bounds)
    if args.dry_run:
        objective(x0)
    elif args.method == "coordinate":
        coordinate_search(objective, x0, args.maxiter, args.initial_step, args.min_step, bounds)
    else:
        run_scipy(objective, x0, args)

    if objective.best_params is not None:
        print(
            f"best score={objective.best_score:.6g} {format_params(fit_params, objective.best_params)}",
            flush=True,
        )
        print(f"log: {objective.log_path}", flush=True)
        print(f"best yaml: {output_dir / 'best.yaml'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
