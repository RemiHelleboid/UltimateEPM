#!/usr/bin/env python3
"""Recursively render every matching PBMC mesh/particle frame as a PNG."""

from __future__ import annotations

import argparse
import logging
import os
import sys
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm, ListedColormap, Normalize

from plot_animation import (
    TimeStep,
    configure_axes,
    load_mesh_frame,
    load_particle_frame,
    zoomed_limits,
)


LOGGER = logging.getLogger("plot_all_frames")
MESH_PREFIX = "mesh_"
PARTICLE_PREFIX = "particles_"


@dataclass(frozen=True)
class FramePair:
    trajectory_dir: Path
    step: str
    mesh_path: Path
    particle_path: Path


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Recursively find trajectory/mesh/mesh_<step>.vtu and "
            "trajectory/particles/particles_<step>.vtp pairs and render each pair "
            "as a PNG."
        )
    )
    parser.add_argument("root", type=Path, help="Root directory to search recursively.")
    parser.add_argument(
        "--output-dir-name",
        default="png",
        help="Output directory created inside each trajectory directory (default: png).",
    )
    parser.add_argument("--field", default="PoissonSolution", help="Mesh point-data field.")
    parser.add_argument(
        "--particle-field",
        default="particle_type",
        help="Particle point-data field containing integer particle types.",
    )
    parser.add_argument("--vmin", type=float, default=0.0, help="Field color minimum.")
    parser.add_argument("--vmax", type=float, default=1.0, help="Field color maximum.")
    parser.add_argument("--dpi", type=int, default=300, help="PNG resolution in DPI.")
    parser.add_argument(
        "--particle-size",
        type=float,
        default=4.0,
        help="Particle marker area in points squared.",
    )
    parser.add_argument(
        "--plane", choices=("xy", "xz", "yz"), default="xy", help="Projection plane."
    )
    parser.add_argument(
        "--zoom",
        type=float,
        default=1.0,
        help="Fraction of the mesh bounds shown around the center.",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Replace PNG files that already exist."
    )
    parser.add_argument(
        "--fail-fast", action="store_true", help="Stop at the first frame that cannot be rendered."
    )
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    return parser.parse_args()


def step_from_name(path: Path, prefix: str, suffix: str) -> str | None:
    name = path.name
    if not name.startswith(prefix) or not name.endswith(suffix):
        return None
    step = name[len(prefix) : -len(suffix)]
    return step or None


def discover_frame_pairs(root: Path) -> tuple[list[FramePair], list[Path]]:
    pairs: list[FramePair] = []
    unmatched: list[Path] = []

    for mesh_dir in sorted(path for path in root.rglob("mesh") if path.is_dir()):
        trajectory_dir = mesh_dir.parent
        particle_dir = trajectory_dir / "particles"
        if trajectory_dir.name != "trajectory" or not particle_dir.is_dir():
            continue

        meshes: dict[str, Path] = {}
        particles: dict[str, Path] = {}
        for path in mesh_dir.glob("mesh_*.vtu"):
            step = step_from_name(path, MESH_PREFIX, ".vtu")
            if step is not None:
                meshes[step] = path
        for path in particle_dir.glob("particles_*.vtp"):
            step = step_from_name(path, PARTICLE_PREFIX, ".vtp")
            if step is not None:
                particles[step] = path

        common_steps = sorted(meshes.keys() & particles.keys())
        pairs.extend(
            FramePair(trajectory_dir, step, meshes[step], particles[step])
            for step in common_steps
        )
        unmatched.extend(meshes[step] for step in sorted(meshes.keys() - particles.keys()))
        unmatched.extend(particles[step] for step in sorted(particles.keys() - meshes.keys()))

    return pairs, unmatched


def render_pair(pair: FramePair, output_path: Path, args: argparse.Namespace) -> None:
    mesh = load_mesh_frame(TimeStep(0.0, pair.mesh_path), args.field, args.plane)
    particles = load_particle_frame(
        TimeStep(0.0, pair.particle_path), args.particle_field, args.plane
    )

    field_norm = Normalize(vmin=args.vmin, vmax=args.vmax, clip=True)
    particle_cmap = ListedColormap(
        np.array([[0.25, 0.0, 1.0, 1.0], [1.0, 0.0, 0.0, 1.0]])
    )
    particle_norm = BoundaryNorm([-0.5, 0.5, 1.5], particle_cmap.N, clip=True)

    figure, ax = plt.subplots(figsize=(11.0, 7.0), constrained_layout=True)
    try:
        figure.patch.set_facecolor("black")
        configure_axes(ax, args.plane)
        x_limits, y_limits = zoomed_limits(mesh.coordinates, args.zoom)
        ax.set_xlim(*x_limits)
        ax.set_ylim(*y_limits)

        triangulation = mtri.Triangulation(
            mesh.coordinates[:, 0], mesh.coordinates[:, 1], mesh.triangles
        )
        ax.tripcolor(
            triangulation,
            mesh.values,
            shading="gouraud",
            cmap="jet",
            norm=field_norm,
            zorder=1,
        )
        ax.scatter(
            particles.coordinates[:, 0],
            particles.coordinates[:, 1],
            c=particles.particle_types,
            cmap=particle_cmap,
            norm=particle_norm,
            s=args.particle_size,
            linewidths=0.1,
            edgecolors="black",
            zorder=2,
        )
        ax.set_title(f"iteration: {pair.step}", color="white", fontsize=18)

        field_mappable = ScalarMappable(norm=field_norm, cmap="jet")
        field_mappable.set_array([])
        field_bar = figure.colorbar(field_mappable, ax=ax, location="bottom", pad=0.02)
        field_bar.set_label(args.field, color="white")
        field_bar.ax.tick_params(colors="white")
        field_bar.outline.set_edgecolor("white")

        output_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_path, dpi=args.dpi, facecolor=figure.get_facecolor())
    finally:
        plt.close(figure)


def main() -> int:
    args = parse_arguments()
    logging.basicConfig(
        level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s"
    )

    root = args.root.expanduser().resolve()
    if not root.is_dir():
        LOGGER.error("Search root is not a directory: %s", root)
        return 1
    output_dir_name = Path(args.output_dir_name)
    if (
        not args.output_dir_name
        or output_dir_name.is_absolute()
        or ".." in output_dir_name.parts
    ):
        LOGGER.error("--output-dir-name must stay inside each trajectory directory.")
        return 1

    pairs, unmatched = discover_frame_pairs(root)
    for path in unmatched:
        LOGGER.warning("No matching mesh/particle file for %s", path)
    if not pairs:
        LOGGER.error("No matching mesh/particle pairs found below %s", root)
        return 1

    LOGGER.info("Found %d matching frame pair(s).", len(pairs))
    rendered = 0
    skipped = 0
    failed = 0
    for index, pair in enumerate(pairs, start=1):
        output_path = (
            pair.trajectory_dir / args.output_dir_name / f"frame_{pair.step}.png"
        )
        if output_path.exists() and not args.overwrite:
            LOGGER.info("Skipping existing %s", output_path)
            skipped += 1
            continue

        try:
            render_pair(pair, output_path, args)
            rendered += 1
            LOGGER.info("Rendered %d/%d: %s", index, len(pairs), output_path)
        except (FileNotFoundError, KeyError, RuntimeError, ValueError) as error:
            failed += 1
            LOGGER.error("Could not render %s: %s", pair.mesh_path, error)
            if args.fail_fast:
                break

    LOGGER.info(
        "Finished: %d rendered, %d skipped, %d failed, %d unmatched.",
        rendered,
        skipped,
        failed,
        len(unmatched),
    )
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
