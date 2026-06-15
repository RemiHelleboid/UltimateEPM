#!/usr/bin/env python3

"""Run the two-pass BZ meshing and band-energy refinement loop."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def positive_float(value: str) -> float:
    parsed = float(value)
    if parsed <= 0.0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def fraction(value: str) -> float:
    parsed = float(value)
    if not 0.0 < parsed <= 1.0:
        raise argparse.ArgumentTypeError("value must be in (0, 1]")
    return parsed


def nonnegative_float(value: str) -> float:
    parsed = float(value)
    if parsed < 0.0:
        raise argparse.ArgumentTypeError("value must be non-negative")
    return parsed


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[2]

    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Iteratively generate an FCC BZ mesh, compute bands in the IBZ, "
            "and remesh from the ranked energy-spread refinement map."
        )
    )
    parser.add_argument(
        "--mesher",
        type=Path,
        default=repo_root / "build/apps/bz_meshing.epm",
        help="Path to bz_meshing.epm.",
    )
    parser.add_argument(
        "--bands",
        type=Path,
        default=repo_root / "build/apps/BandsOnBZ",
        help="Path to BandsOnBZ.",
    )
    parser.add_argument(
        "--epm-set",
        default="local-cohen",
        help="Named EPM parameter set.",
    )
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory.")
    parser.add_argument("--prefix", default="bz", help="Iteration filename prefix.")
    parser.add_argument("--material", default="Si", help="Material symbol.")
    parser.add_argument(
        "--iterations",
        type=positive_int,
        default=3,
        help="Maximum number of band-evaluation iterations, including iteration 0.",
    )
    parser.add_argument("--threads", type=positive_int, default=1, help="Band solver threads.")
    parser.add_argument("--valence-bands", type=int, default=4, help="Top valence bands to diagnose/export.")
    parser.add_argument("--conduction-bands", type=int, default=8, help="Lowest conduction bands to diagnose/export.")
    parser.add_argument("--nearest-neighbors", type=positive_int, default=10, help="EPM nearest neighbors.")
    parser.add_argument(
        "--energy-target",
        type=positive_float,
        default=0.020,
        help="Energy-spread target in eV.",
    )
    parser.add_argument(
        "--max-refinement-points",
        type=positive_int,
        default=5000,
        help="Absolute maximum refinement points exported per iteration.",
    )
    parser.add_argument(
        "--max-refinement-fraction",
        type=fraction,
        default=0.10,
        help="Maximum new refinement points as a fraction of the current IBZ nodes.",
    )
    parser.add_argument(
        "--min-refinement-points",
        type=positive_int,
        default=100,
        help="Minimum batch allowance before the fractional cap becomes active.",
    )
    parser.add_argument(
        "--max-new-refinement-points",
        type=positive_int,
        default=1000,
        help="Hard cap on new refinement points added in one iteration.",
    )
    parser.add_argument(
        "--max-cumulative-refinement-points",
        type=positive_int,
        default=20000,
        help="Maximum refinement points retained across all iterations.",
    )
    parser.add_argument(
        "--retained-history-factor",
        type=nonnegative_float,
        default=1.0,
        help=(
            "Maximum older refinement points retained per newly selected point. "
            "Current ranked points always take priority."
        ),
    )
    parser.add_argument(
        "--adaptive-radius-factor",
        type=positive_float,
        default=4.0,
        help="Adaptive field radius divided by target mesh size.",
    )
    parser.add_argument(
        "--mesh-size",
        type=positive_float,
        default=0.10,
        help="Coarse uniform mesh size used as the background throughout the loop.",
    )
    parser.add_argument("--nonlocal-correction", action="store_true", help="Enable the EPM non-local correction.")
    parser.add_argument("--soc", action="store_true", help="Enable spin-orbit coupling.")
    parser.add_argument(
        "--max-full-nodes",
        type=positive_int,
        default=1_000_000,
        help="Stop before band evaluation if a generated mesh exceeds this node count.",
    )
    parser.add_argument(
        "--max-full-tets",
        type=positive_int,
        default=6_000_000,
        help="Stop before band evaluation if a generated mesh exceeds this tetrahedron count.",
    )
    parser.add_argument(
        "--max-mesh-growth-factor",
        type=positive_float,
        default=3.0,
        help="Stop before band evaluation if tetrahedra grow by more than this factor in one iteration.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse complete iteration outputs already present in the output directory.",
    )
    return parser.parse_args()


def checked_path(path: Path, description: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{description} does not exist: {resolved}")
    return resolved


def run_logged(
    command: list[str],
    cwd: Path,
    log_path: Path,
    *,
    input_digest: str | None = None,
    check: bool = True,
) -> int:
    printable = " ".join(command)
    print(f"\n$ {printable}", flush=True)

    with log_path.open("w", encoding="utf-8") as log:
        log.write(f"$ {printable}\n\n")
        if input_digest is not None:
            log.write(f"refinement_map_sha256={input_digest}\n\n")
        process = subprocess.Popen(
            command,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            log.write(line)

        return_code = process.wait()

    if check and return_code != 0:
        raise subprocess.CalledProcessError(return_code, command)
    return return_code


def read_metadata(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            key, separator, value = line.partition("=")
            if separator:
                values[key.strip()] = value.strip()
    return values


def read_adaptive_summary(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        row = next(reader, None)
        if row is None:
            raise ValueError(f"adaptive summary has no data row: {path}")
        return row


def log_matches_command(path: Path, command: list[str], input_digest: str | None = None) -> bool:
    if not path.is_file():
        return False
    with path.open(encoding="utf-8") as stream:
        first_line = stream.readline().rstrip("\n")
        if first_line != f"$ {' '.join(command)}":
            return False
        if input_digest is None:
            return True
        return any(line.rstrip("\n") == f"refinement_map_sha256={input_digest}" for line in stream)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_refinement_map(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames:
            raise ValueError(f"refinement map has no CSV header: {path}")
        return reader.fieldnames, list(reader)


def write_refinement_map(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def refinement_point_count(path: Path) -> int:
    _, rows = read_refinement_map(path)
    return len(rows)


def truncate_refinement_map(path: Path, limit: int) -> int:
    fieldnames, rows = read_refinement_map(path)
    selected = rows[:limit]
    write_refinement_map(path, fieldnames, selected)
    return len(selected)


def refinement_backoff_path(path: Path) -> Path:
    return path.with_suffix(path.suffix + ".backoff.json")


def read_persisted_refinement_limit(path: Path, candidate_digest: str) -> int | None:
    if not path.is_file():
        return None
    state = json.loads(path.read_text(encoding="utf-8"))
    if state.get("candidate_sha256") != candidate_digest:
        return None
    value = int(state["limit"])
    if value <= 0:
        raise ValueError(f"invalid persisted refinement limit in {path}")
    return value


def persist_refinement_limit(
    path: Path,
    limit: int,
    candidate_digest: str,
    accepted_digest: str,
) -> None:
    state = {
        "candidate_sha256": candidate_digest,
        "accepted_sha256": accepted_digest,
        "limit": limit,
    }
    path.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")


def next_backoff_limit(current: int, current_batch: int, minimum: int) -> int | None:
    if current <= 1:
        return None
    if current > current_batch:
        return max(1, current_batch)
    if current > minimum:
        return max(minimum, current // 2)
    return max(1, current // 2)


def select_refinement_points(source: Path, destination: Path, limit: int) -> tuple[int, int]:
    fieldnames, rows = read_refinement_map(source)
    selected = rows[:limit]
    write_refinement_map(destination, fieldnames, selected)

    return len(rows), len(selected)


def merge_refinement_points(
    previous: Path | None,
    selected: Path,
    destination: Path,
    max_points: int,
) -> tuple[int, int]:
    selected_fieldnames, selected_rows = read_refinement_map(selected)
    previous_rows: list[dict[str, str]] = []
    if previous is not None:
        previous_fieldnames, previous_rows = read_refinement_map(previous)
        if previous_fieldnames != selected_fieldnames:
            raise ValueError("refinement map CSV headers do not match")

    rows = selected_rows + previous_rows
    fieldnames = selected_fieldnames
    if not selected_rows:
        write_refinement_map(destination, fieldnames, [])
        return len(rows), 0

    # Current energy errors are authoritative. Older points only fill the
    # remaining budget and cannot displace a newly selected point.
    accepted: list[dict[str, str]] = []
    accepted_is_current: list[bool] = []
    spatial_bins: dict[tuple[int, int, int], list[int]] = {}
    bin_size = max(float(row["target_h"]) for row in rows)
    if not math.isfinite(bin_size) or bin_size <= 0.0:
        raise ValueError("refinement map contains an invalid target_h")
    bin_size *= 0.25

    for row_index, row in enumerate(rows):
        is_current = row_index < len(selected_rows)
        point = (float(row["x"]), float(row["y"]), float(row["z"]))
        target_h = float(row["target_h"])
        key = tuple(math.floor(coordinate / bin_size) for coordinate in point)
        merge_index: int | None = None

        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for index in spatial_bins.get((key[0] + dx, key[1] + dy, key[2] + dz), []):
                        other = accepted[index]
                        distance = math.sqrt(
                            (point[0] - float(other["x"])) ** 2
                            + (point[1] - float(other["y"])) ** 2
                            + (point[2] - float(other["z"])) ** 2
                        )
                        merge_radius = 0.25 * min(target_h, float(other["target_h"]))
                        if distance <= merge_radius:
                            merge_index = index
                            break
                    if merge_index is not None:
                        break
                if merge_index is not None:
                    break
            if merge_index is not None:
                break

        if merge_index is not None:
            if accepted_is_current[merge_index] and not is_current:
                continue
            retained = accepted[merge_index]
            if target_h < float(retained["target_h"]):
                retained["target_h"] = row["target_h"]
            if float(row["error_eV"]) > float(retained["error_eV"]):
                retained["error_eV"] = row["error_eV"]
                retained["band"] = row["band"]
                retained["tetra"] = row["tetra"]
            continue

        if len(accepted) >= max_points:
            continue

        spatial_bins.setdefault(key, []).append(len(accepted))
        accepted.append(row)
        accepted_is_current.append(is_current)

    write_refinement_map(destination, fieldnames, accepted)

    return len(rows), len(accepted)


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "iteration",
        "mesh",
        "full_nodes",
        "full_tetrahedra",
        "ibz_nodes",
        "ibz_tetrahedra",
        "refinement_point_limit",
        "refinement_points_found",
        "refinement_points_exported",
        "refinement_points_used",
        "cumulative_refinement_points",
        "violating_volume_fraction",
        "p50_error_eV",
        "p90_error_eV",
        "p99_error_eV",
        "max_error_eV",
        "converged",
    ]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()

    if args.valence_bands < 0 or args.conduction_bands < 0:
        raise ValueError("band counts cannot be negative")
    if args.valence_bands + args.conduction_bands == 0:
        raise ValueError("at least one band must be selected")
    if args.adaptive_radius_factor <= 1.0:
        raise ValueError("--adaptive-radius-factor must be greater than 1")

    mesher = checked_path(args.mesher, "mesher executable")
    bands = checked_path(args.bands, "band executable")
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, object]] = []
    previous_cumulative_map: Path | None = None
    previous_full_tets: int | None = None
    converged = False

    for iteration in range(args.iterations):
        stem = f"{args.prefix}_iter{iteration:02d}"
        mesh_name = f"{stem}.msh"
        metadata_name = f"{stem}_metadata.txt"
        kstar_name = f"{stem}_kstar_ibz_to_bz.txt"
        refinement_name = f"{stem}_refinement.csv"
        refinement_summary_name = f"{stem}_refinement_summary.csv"
        selected_refinement_name = f"{stem}_refinement_selected.csv"
        cumulative_refinement_name = f"{stem}_refinement_cumulative.csv"
        bands_name = f"{stem}_bands.msh"

        mesh_path = outdir / mesh_name
        metadata_path = outdir / metadata_name
        kstar_path = outdir / kstar_name
        refinement_path = outdir / refinement_name
        refinement_summary_path = outdir / refinement_summary_name
        selected_refinement_path = outdir / selected_refinement_name
        cumulative_refinement_path = outdir / cumulative_refinement_name
        bands_path = outdir / bands_name
        mesher_log_path = outdir / f"{stem}_mesher.log"

        max_ibz_nodes = max(1, args.max_full_nodes // 48)
        max_ibz_tets = max(1, args.max_full_tets // 48)
        if iteration >= 2 and previous_full_tets is not None:
            growth_full_tet_cap = int(args.max_mesh_growth_factor * previous_full_tets)
            max_ibz_tets = min(max_ibz_tets, max(1, growth_full_tet_cap // 48))

        mesh_command = [
            str(mesher),
            "--mesh",
            str(args.mesh_size),
            "--nogui",
            "--outfile",
            mesh_name,
        ]
        if previous_cumulative_map is None:
            mesh_command.append("--uniform")
        else:
            mesh_command.extend(
                [
                    "--delta-h",
                    str(args.mesh_size),
                    "--no-tube",
                    "--no-L",
                    "--refinement-map",
                    previous_cumulative_map.name,
                    "--adaptive-radius-factor",
                    str(args.adaptive_radius_factor),
                ]
            )
        mesh_command.extend(
            [
                "--max-ibz-nodes",
                str(max_ibz_nodes),
                "--max-ibz-tets",
                str(max_ibz_tets),
            ]
        )

        mesh_was_generated = False
        while True:
            refinement_digest = (
                file_sha256(previous_cumulative_map) if previous_cumulative_map is not None else None
            )
            mesh_complete = (
                mesh_path.is_file()
                and metadata_path.is_file()
                and kstar_path.is_file()
                and log_matches_command(mesher_log_path, mesh_command, refinement_digest)
            )
            if args.resume and mesh_complete:
                print(f"\nReusing mesh iteration {iteration}: {mesh_path}")
                break

            mesh_was_generated = True
            mesh_path.unlink(missing_ok=True)
            metadata_path.unlink(missing_ok=True)
            kstar_path.unlink(missing_ok=True)
            return_code = run_logged(
                mesh_command,
                outdir,
                mesher_log_path,
                input_digest=refinement_digest,
                check=False,
            )
            if return_code == 0:
                break
            if return_code != 3 or previous_cumulative_map is None:
                raise subprocess.CalledProcessError(return_code, mesh_command)

            current_count = refinement_point_count(previous_cumulative_map)
            current_batch = int(summary_rows[-1]["refinement_points_used"])
            reduced_limit = next_backoff_limit(
                current_count,
                current_batch,
                args.min_refinement_points,
            )
            if reduced_limit is None or reduced_limit >= current_count:
                print(
                    f"Stopping: the mesh exceeds its IBZ safety cap even with "
                    f"{current_count} refinement point."
                )
                return 2

            backoff_path = refinement_backoff_path(previous_cumulative_map)
            current_digest = file_sha256(previous_cumulative_map)
            candidate_digest = current_digest
            if backoff_path.is_file():
                state = json.loads(backoff_path.read_text(encoding="utf-8"))
                if state.get("accepted_sha256") == current_digest:
                    candidate_digest = state.get("candidate_sha256", current_digest)
            retained_count = truncate_refinement_map(previous_cumulative_map, reduced_limit)
            persist_refinement_limit(
                backoff_path,
                retained_count,
                candidate_digest,
                file_sha256(previous_cumulative_map),
            )
            summary_rows[-1]["cumulative_refinement_points"] = retained_count
            write_summary(outdir / f"{args.prefix}_adaptive_summary.csv", summary_rows)
            print(
                f"Retrying mesh iteration {iteration}: reducing the refinement map "
                f"from {current_count} to {retained_count} points."
            )

        metadata = read_metadata(metadata_path)
        full_nodes = int(metadata["full_nodes"])
        full_tets = int(metadata["full_tetrahedra"])
        ibz_nodes = int(metadata["ibz_nodes"])
        ibz_tets = int(metadata["ibz_tetrahedra"])
        if full_nodes > args.max_full_nodes or full_tets > args.max_full_tets:
            print(
                f"Stopping: {mesh_name} has {full_nodes} nodes and {full_tets} tetrahedra, "
                "which exceeds the configured safety cap despite the IBZ pre-check."
            )
            return 2
        if iteration >= 2 and previous_full_tets is not None and full_tets > args.max_mesh_growth_factor * previous_full_tets:
            print(
                f"Stopping: {mesh_name} grew from {previous_full_tets} to {full_tets} tetrahedra "
                f"({full_tets / previous_full_tets:.2f}x), above the "
                f"{args.max_mesh_growth_factor:.2f}x growth cap despite the IBZ pre-check."
            )
            return 2

        refinement_point_limit = min(
            args.max_refinement_points,
            args.max_new_refinement_points,
            max(args.min_refinement_points, int(args.max_refinement_fraction * ibz_nodes)),
        )
        print(
            f"Adaptive point budget: {refinement_point_limit} "
            f"(IBZ nodes={ibz_nodes}, fraction={args.max_refinement_fraction:.3f}, "
            f"hard cap={args.max_new_refinement_points})."
        )

        band_command = [
            str(bands),
            "--meshfile",
            mesh_name,
            "--material",
            args.material,
            "--epm-set",
            args.epm_set,
            "--IrrWedge",
            "--nvbands",
            str(args.valence_bands),
            "--ncbands",
            str(args.conduction_bands),
            "--nearestNeighbors",
            str(args.nearest_neighbors),
            "--nthreads",
            str(args.threads),
            "--adaptive-energy-target",
            str(args.energy_target),
            "--adaptive-max-points",
            str(args.max_refinement_points),
            "--adaptive-background-h",
            str(args.mesh_size),
            "--refinement-map",
            refinement_name,
            "--outfile",
            bands_name,
        ]
        if args.nonlocal_correction:
            band_command.append("--nonlocal-correction")
        if args.soc:
            band_command.append("--soc")

        bands_log_path = outdir / f"{stem}_bands.log"
        bands_complete = (
            refinement_path.is_file()
            and refinement_summary_path.is_file()
            and bands_path.is_file()
            and log_matches_command(bands_log_path, band_command)
            and not mesh_was_generated
        )
        if not (args.resume and bands_complete):
            bands_path.unlink(missing_ok=True)
            run_logged(band_command, outdir, bands_log_path)
        else:
            print(f"Reusing band results for iteration {iteration}: {bands_path}")

        refinement_points_exported, refinement_points_used = select_refinement_points(
            refinement_path,
            selected_refinement_path,
            refinement_point_limit,
        )
        adaptive_summary = read_adaptive_summary(refinement_summary_path)
        refinement_points_found = int(adaptive_summary["violating_tetrahedra"])
        violating_volume_fraction = float(adaptive_summary["violating_volume_fraction"])
        p50_error = float(adaptive_summary["p50_error_eV"])
        p90_error = float(adaptive_summary["p90_error_eV"])
        p99_error = float(adaptive_summary["p99_error_eV"])
        max_error = float(adaptive_summary["max_error_eV"])
        iteration_converged = refinement_points_found == 0
        cumulative_limit = min(
            args.max_cumulative_refinement_points,
            refinement_points_used + int(args.retained_history_factor * refinement_points_used),
        )
        _, cumulative_refinement_points = merge_refinement_points(
            previous_cumulative_map,
            selected_refinement_path,
            cumulative_refinement_path,
            cumulative_limit,
        )
        backoff_path = refinement_backoff_path(cumulative_refinement_path)
        candidate_digest = file_sha256(cumulative_refinement_path)
        persisted_limit = (
            read_persisted_refinement_limit(backoff_path, candidate_digest) if args.resume else None
        )
        if persisted_limit is not None and persisted_limit < cumulative_refinement_points:
            cumulative_refinement_points = truncate_refinement_map(
                cumulative_refinement_path,
                persisted_limit,
            )
        elif not args.resume:
            backoff_path.unlink(missing_ok=True)
        summary_rows.append(
            {
                "iteration": iteration,
                "mesh": mesh_name,
                "full_nodes": full_nodes,
                "full_tetrahedra": full_tets,
                "ibz_nodes": ibz_nodes,
                "ibz_tetrahedra": ibz_tets,
                "refinement_point_limit": refinement_point_limit,
                "refinement_points_found": refinement_points_found,
                "refinement_points_exported": refinement_points_exported,
                "refinement_points_used": refinement_points_used,
                "cumulative_refinement_points": cumulative_refinement_points,
                "violating_volume_fraction": violating_volume_fraction,
                "p50_error_eV": p50_error,
                "p90_error_eV": p90_error,
                "p99_error_eV": p99_error,
                "max_error_eV": max_error,
                "converged": iteration_converged,
            }
        )
        write_summary(outdir / f"{args.prefix}_adaptive_summary.csv", summary_rows)

        print(
            f"Iteration {iteration}: {full_nodes} full-BZ nodes, {full_tets} tetrahedra, "
            f"{refinement_points_found}/{ibz_tets} violating IBZ tetrahedra "
            f"({100.0 * violating_volume_fraction:.2f}% of IBZ volume), "
            f"{refinement_points_exported} candidates exported, "
            f"{refinement_points_used} new points selected, "
            f"{cumulative_refinement_points} cumulative points retained; "
            f"p50={p50_error:.6f} eV, max={max_error:.6f} eV."
        )
        if iteration_converged:
            print("Converged: no tetrahedra exceeded the energy-spread target.")
            converged = True
            break

        previous_cumulative_map = cumulative_refinement_path
        previous_full_tets = full_tets

    summary_path = outdir / f"{args.prefix}_adaptive_summary.csv"
    if converged:
        print(f"\nAdaptive run converged. Summary: {summary_path}")
        return 0

    last = summary_rows[-1]
    print(
        f"\nStopped without convergence after {args.iterations} iterations: "
        f"{last['refinement_points_found']} IBZ tetrahedra still exceed "
        f"{args.energy_target:.6f} eV, covering "
        f"{100.0 * float(last['violating_volume_fraction']):.2f}% of the IBZ volume. "
        f"Summary: {summary_path}"
    )
    return 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        print(f"Error: {error}", file=sys.stderr)
        raise SystemExit(1)
