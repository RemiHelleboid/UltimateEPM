#!/usr/bin/env python3
"""Fit local EPM form factors by driving the C++ EPM analysis apps.

The script edits a generated EPM YAML parameter set, runs:
  - epm_band_edges for gap and valley-position targets
  - epm_valley_fit for effective masses and Kane non-parabolicity

Then it minimizes a weighted normalized least-squares score.
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class Target:
    value: float
    scale: float
    weight: float = 1.0


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--material", default="Si")
    parser.add_argument("--base-set", default="local-remi")
    parser.add_argument("--work-set", default="local-fit-working")
    parser.add_argument("--build-dir", default=str(REPO_ROOT / "build"))
    parser.add_argument("--output-dir", default=str(REPO_ROOT / "fit_epm_output"))
    parser.add_argument("--maxiter", type=int, default=60)
    parser.add_argument("--method", choices=["auto", "powell", "nelder-mead", "coordinate"], default="auto")
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
    parser.add_argument("--initial-step", type=float, default=0.01, help="Initial coordinate-search step in Rydberg.")
    parser.add_argument("--min-step", type=float, default=2.0e-4, help="Minimum coordinate-search step in Rydberg.")
    parser.add_argument("--dry-run", action="store_true", help="Evaluate the initial parameter set once and exit.")
    return parser.parse_args()


def parameter_file(material: str, parameter_set: str) -> Path:
    return REPO_ROOT / "data" / "materials" / material / "epm" / f"{parameter_set}.yaml"


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as stream:
        return yaml.safe_load(stream)


def write_working_yaml(base_config: dict, material: str, work_set: str, params: list[float]) -> Path:
    config = dict(base_config)
    config["parameter_set"] = work_set
    pseudo = dict(config["pseudo-potential-parameters"])
    pseudo["V3S"] = float(params[0])
    pseudo["V8S"] = float(params[1])
    pseudo["V11S"] = float(params[2])
    config["pseudo-potential-parameters"] = pseudo

    path = parameter_file(material, work_set)
    with path.open("w", encoding="utf-8") as stream:
        yaml.safe_dump(config, stream, sort_keys=False)
    return path


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


def as_float(values: dict[str, str], key: str) -> float:
    return float(values[key])


def run_command(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, cwd=cwd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)


def measure(
    args: argparse.Namespace,
    work_set: str,
    output_dir: Path,
    iteration: int,
) -> dict[str, float]:
    edges_csv = output_dir / f"edges_{iteration:04d}.csv"
    valley_csv = output_dir / f"valley_{iteration:04d}.csv"
    build_apps = Path(args.build_dir) / "apps"

    run_command(
        [
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
        ],
        REPO_ROOT,
    )

    edges = read_quantity_csv(edges_csv)
    k_delta = as_float(edges, "delta_cbm_kx_reduced")

    run_command(
        [
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
        ],
        REPO_ROOT,
    )

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


def score_observables(observables: dict[str, float]) -> float:
    score = 0.0
    for key, target in DEFAULT_TARGETS.items():
        residual = (observables[key] - target.value) / target.scale
        score += target.weight * residual * residual
    if observables["alpha_rms_error_meV"] > 5.0:
        score += ((observables["alpha_rms_error_meV"] - 5.0) / 5.0) ** 2
    return score


class Objective:
    def __init__(self, args: argparse.Namespace, base_config: dict, output_dir: Path) -> None:
        self.args = args
        self.base_config = base_config
        self.output_dir = output_dir
        self.iteration = 0
        self.best_score = math.inf
        self.best_params: list[float] | None = None
        self.log_path = output_dir / "iterations.csv"
        with self.log_path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                [
                    "iteration",
                    "score",
                    "V3S",
                    "V8S",
                    "V11S",
                    *DEFAULT_TARGETS.keys(),
                    "mass_rms_error_meV",
                    "alpha_rms_error_meV",
                ]
            )

    def __call__(self, params_like) -> float:
        params = [float(x) for x in params_like]
        iteration = self.iteration
        self.iteration += 1
        write_working_yaml(self.base_config, self.args.material, self.args.work_set, params)
        try:
            observables = measure(self.args, self.args.work_set, self.output_dir, iteration)
            score = score_observables(observables)
        except Exception as exc:
            print(f"iteration {iteration:04d} failed for {params}: {exc}", flush=True)
            score = 1.0e12
            observables = {key: math.nan for key in DEFAULT_TARGETS}
            observables["mass_rms_error_meV"] = math.nan
            observables["alpha_rms_error_meV"] = math.nan

        with self.log_path.open("a", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                [
                    iteration,
                    score,
                    params[0],
                    params[1],
                    params[2],
                    *(observables[key] for key in DEFAULT_TARGETS),
                    observables["mass_rms_error_meV"],
                    observables["alpha_rms_error_meV"],
                ]
            )

        if score < self.best_score:
            self.best_score = score
            self.best_params = params
            shutil.copy2(parameter_file(self.args.material, self.args.work_set), self.output_dir / "best.yaml")

        print(
            f"iter {iteration:04d} score={score:.6g} "
            f"V3S={params[0]:.7g} V8S={params[1]:.7g} V11S={params[2]:.7g}",
            flush=True,
        )
        return score


def coordinate_search(objective: Callable[[list[float]], float], x0: list[float], maxiter: int, step0: float, min_step: float):
    x = list(x0)
    best = objective(x)
    step = step0
    iterations = 1
    while iterations < maxiter and step >= min_step:
        improved = False
        for dim in range(len(x)):
            for direction in (1.0, -1.0):
                trial = list(x)
                trial[dim] += direction * step
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
        from scipy.optimize import minimize
    except ImportError:
        print("SciPy not found; falling back to coordinate search.", flush=True)
        coordinate_search(objective, x0, args.maxiter, args.initial_step, args.min_step)
        return

    method = "Powell" if args.method in ("auto", "powell") else "Nelder-Mead"
    options = {"maxiter": args.maxiter, "disp": True}
    if method == "Powell":
        options["xtol"] = args.min_step
        options["ftol"] = 1.0e-4
    else:
        options["xatol"] = args.min_step
        options["fatol"] = 1.0e-4
    minimize(objective, x0, method=method, options=options)


def main() -> int:
    args = parse_args()
    if args.nthreads <= 0:
        raise SystemExit("--nthreads must be positive")

    base_path = parameter_file(args.material, args.base_set)
    base_config = load_yaml(base_path)
    pseudo = base_config["pseudo-potential-parameters"]
    x0 = [float(pseudo["V3S"]), float(pseudo["V8S"]), float(pseudo["V11S"])]

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    objective = Objective(args, base_config, output_dir)
    if args.dry_run:
        objective(x0)
    elif args.method == "coordinate":
        coordinate_search(objective, x0, args.maxiter, args.initial_step, args.min_step)
    else:
        run_scipy(objective, x0, args)

    if objective.best_params is not None:
        print(
            f"best score={objective.best_score:.6g} "
            f"V3S={objective.best_params[0]:.7g} "
            f"V8S={objective.best_params[1]:.7g} "
            f"V11S={objective.best_params[2]:.7g}",
            flush=True,
        )
        print(f"log: {objective.log_path}", flush=True)
        print(f"best yaml: {output_dir / 'best.yaml'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
