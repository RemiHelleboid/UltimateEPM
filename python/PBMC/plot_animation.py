#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from matplotlib.colors import BoundaryNorm, ListedColormap, Normalize
from matplotlib.cm import ScalarMappable

LOGGER = logging.getLogger("animate_scene")


@dataclass(frozen=True)
class TimeStep:
    time: float
    path: Path


@dataclass
class MeshFrame:
    coordinates: np.ndarray
    triangles: np.ndarray
    values: np.ndarray


@dataclass
class ParticleFrame:
    coordinates: np.ndarray
    particle_types: np.ndarray


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render a ParaView PVD time series with Matplotlib. The mesh is colored "
            "by a point-data field and particles are overlaid as a scatter plot."
        )
    )
    parser.add_argument(
        "scene_dir",
        type=Path,
        help="Directory containing mesh/mesh.pvd and particles/particles.pvd.",
    )
    parser.add_argument(
        "--mesh-pvd",
        type=Path,
        default=Path("mesh/mesh.pvd"),
        help="Mesh PVD path, relative to scene_dir unless absolute.",
    )
    parser.add_argument(
        "--particles-pvd",
        type=Path,
        default=Path("particles/particles.pvd"),
        help="Particle PVD path, relative to scene_dir unless absolute.",
    )
    parser.add_argument(
        "--field",
        default="PoissonSolution",
        help="Mesh point-data field to display.",
    )
    parser.add_argument(
        "--particle-field",
        default="particle_type",
        help="Particle point-data field containing integer particle types.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("animation.mp4"),
        help="Output .mp4 or .gif file.",
    )
    parser.add_argument("--fps", type=float, default=30.0, help="Output frames per second.")
    parser.add_argument("--dpi", type=int, default=300, help="Output resolution in DPI.")
    parser.add_argument(
        "--interval",
        type=float,
        default=None,
        help="Interactive frame interval in milliseconds. Defaults to 1000/fps.",
    )
    parser.add_argument("--vmin", type=float, default=0.0, help="Field color minimum.")
    parser.add_argument("--vmax", type=float, default=1.0, help="Field color maximum.")
    parser.add_argument(
        "--particle-size",
        type=float,
        default=14.0,
        help="Particle marker area in points squared.",
    )
    parser.add_argument(
        "--plane",
        choices=("xy", "xz", "yz"),
        default="xy",
        help="Coordinate plane used for the 2D projection.",
    )
    parser.add_argument(
        "--zoom",
        type=float,
        default=1.0,
        help="Fraction of the mesh bounds shown around the center.",
    )
    parser.add_argument(
        "--time-tolerance",
        type=float,
        default=None,
        help="Absolute tolerance used to match particle timesteps to mesh timesteps.",
    )
    parser.add_argument(
        "--min-iteration",
        type=int,
        default=0,
        help="First PVD iteration to animate, using zero-based indexing (default: 0).",
    )
    parser.add_argument(
        "--max-iteration",
        type=int,
        default=None,
        help="Last PVD iteration to animate, inclusive (default: final iteration).",
    )
    parser.add_argument(
        "--iterstep",
        type=int,
        default=1,
        help="Iteration stride between animated frames (default: 1).",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Open an interactive window after saving the animation.",
    )
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    return parser.parse_args()


def resolve_path(base_directory: Path, path: Path) -> Path:
    return path if path.is_absolute() else base_directory / path


def select_iterations(
    entries: list[TimeStep],
    min_iteration: int,
    max_iteration: int | None,
    iterstep: int,
) -> tuple[list[int], list[TimeStep]]:
    if min_iteration < 0:
        raise ValueError("--min-iteration must be greater than or equal to zero.")
    if max_iteration is not None and max_iteration < 0:
        raise ValueError("--max-iteration must be greater than or equal to zero.")
    if iterstep <= 0:
        raise ValueError("--iterstep must be greater than zero.")
    if min_iteration >= len(entries):
        raise ValueError(
            f"--min-iteration={min_iteration} is outside the available range "
            f"0..{len(entries) - 1}."
        )

    last_iteration = len(entries) - 1 if max_iteration is None else max_iteration
    if last_iteration >= len(entries):
        raise ValueError(
            f"--max-iteration={last_iteration} is outside the available range "
            f"0..{len(entries) - 1}."
        )
    if last_iteration < min_iteration:
        raise ValueError("--max-iteration must be greater than or equal to --min-iteration.")

    indices = list(range(min_iteration, last_iteration + 1, iterstep))
    return indices, [entries[index] for index in indices]


def read_pvd(path: Path) -> list[TimeStep]:
    if not path.is_file():
        raise FileNotFoundError(f"PVD file does not exist: {path}")

    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as error:
        raise ValueError(f"Invalid PVD XML in {path}: {error}") from error

    dataset_elements = root.findall(".//DataSet")
    if not dataset_elements:
        dataset_elements = root.findall(".//{*}DataSet")

    entries: list[TimeStep] = []
    for element in dataset_elements:
        file_name = element.get("file")
        if not file_name:
            continue

        time_text = element.get("timestep", "0")
        try:
            time = float(time_text)
        except ValueError as error:
            raise ValueError(
                f"Invalid timestep {time_text!r} in PVD file {path}."
            ) from error

        data_path = (path.parent / unquote(file_name)).resolve()
        entries.append(TimeStep(time=time, path=data_path))

    if not entries:
        raise ValueError(f"No DataSet entries were found in {path}.")

    entries.sort(key=lambda entry: entry.time)
    duplicate_times = [
        entries[index].time
        for index in range(1, len(entries))
        if entries[index].time == entries[index - 1].time
    ]
    if duplicate_times:
        raise ValueError(
            f"{path} contains multiple DataSet entries for the same timestep. "
            "This script expects one data file per timestep."
        )

    return entries


def create_vtk_reader(path: Path):
    readers = {
        ".vtp": vtk.vtkXMLPolyDataReader,
        ".vtu": vtk.vtkXMLUnstructuredGridReader,
        ".pvtp": vtk.vtkXMLPPolyDataReader,
        ".pvtu": vtk.vtkXMLPUnstructuredGridReader,
        ".vtk": vtk.vtkDataSetReader,
    }

    reader_type = readers.get(path.suffix.lower())
    if reader_type is None:
        supported = ", ".join(sorted(readers))
        raise ValueError(
            f"Unsupported VTK data file {path}. Supported extensions: {supported}."
        )

    reader = reader_type()
    reader.SetFileName(str(path))
    return reader


def read_vtk_dataset(path: Path):
    if not path.is_file():
        raise FileNotFoundError(f"Referenced data file does not exist: {path}")

    reader = create_vtk_reader(path)
    try:
        reader.Update()
    except Exception as error:
        raise RuntimeError(f"Could not read VTK data file {path}: {error}") from error

    error_code = reader.GetErrorCode() if hasattr(reader, "GetErrorCode") else 0
    if error_code:
        error_name = vtk.vtkErrorCode.GetStringFromErrorCode(error_code)
        raise RuntimeError(f"Could not read VTK data file {path}: {error_name}")

    output = reader.GetOutputDataObject(0)
    dataset = vtk.vtkDataSet.SafeDownCast(output)
    if dataset is None:
        output_type = output.GetClassName() if output is not None else "none"
        raise RuntimeError(
            f"Could not read VTK data file {path}: expected vtkDataSet, got {output_type}."
        )

    return dataset


def dataset_points(dataset, path: Path) -> np.ndarray:
    points = dataset.GetPoints()
    if points is None:
        if dataset.GetNumberOfPoints() == 0:
            return np.empty((0, 3), dtype=float)
        raise ValueError(f"Dataset contains no point coordinates: {path}")

    coordinates = np.asarray(vtk_to_numpy(points.GetData()), dtype=float)
    if coordinates.ndim != 2:
        raise ValueError(
            f"Expected a two-dimensional point array in {path}, got {coordinates.shape}."
        )
    return coordinates


def available_point_data_names(dataset) -> list[str]:
    point_data = dataset.GetPointData()
    return [
        point_data.GetArrayName(index) or f"array_{index}"
        for index in range(point_data.GetNumberOfArrays())
    ]


def point_data_array(dataset, field_name: str, path: Path) -> np.ndarray:
    point_data = dataset.GetPointData()
    vtk_array = point_data.GetArray(field_name)

    if vtk_array is None:
        if dataset.GetNumberOfPoints() == 0:
            return np.empty(0, dtype=float)
        available = ", ".join(sorted(available_point_data_names(dataset))) or "none"
        raise KeyError(
            f"Point-data field {field_name!r} is missing from {path}. "
            f"Available fields: {available}."
        )

    values = np.asarray(vtk_to_numpy(vtk_array)).squeeze()
    if values.ndim != 1:
        raise ValueError(
            f"Point-data field {field_name!r} in {path} is not scalar; "
            f"shape={values.shape}."
        )
    if values.shape[0] != dataset.GetNumberOfPoints():
        raise ValueError(
            f"Point-data field {field_name!r} in {path} has {values.shape[0]} values "
            f"for {dataset.GetNumberOfPoints()} points."
        )
    return values


def triangulated_surface(dataset, path: Path):
    surface_filter = vtk.vtkDataSetSurfaceFilter()
    surface_filter.SetInputData(dataset)

    triangle_filter = vtk.vtkTriangleFilter()
    triangle_filter.SetInputConnection(surface_filter.GetOutputPort())
    triangle_filter.PassVertsOff()
    triangle_filter.PassLinesOff()
    triangle_filter.Update()

    surface = triangle_filter.GetOutput()
    if surface is None or surface.GetNumberOfPolys() == 0:
        raise ValueError(f"The mesh contains no polygonal surface cells: {path}")
    return surface


def surface_triangles(surface, path: Path) -> np.ndarray:
    polygons = surface.GetPolys()
    triangles: list[tuple[int, int, int]] = []
    point_ids = vtk.vtkIdList()

    polygons.InitTraversal()
    while polygons.GetNextCell(point_ids):
        if point_ids.GetNumberOfIds() != 3:
            raise ValueError(
                f"Triangle conversion produced a non-triangular cell in {path}."
            )
        triangles.append(
            (
                point_ids.GetId(0),
                point_ids.GetId(1),
                point_ids.GetId(2),
            )
        )

    if not triangles:
        raise ValueError(f"The mesh contains no triangles after conversion: {path}")

    return np.ascontiguousarray(triangles, dtype=np.int64)

def project_points(points: np.ndarray, plane: str) -> np.ndarray:
    if points.ndim != 2 or points.shape[1] < 2:
        raise ValueError(f"Expected an N x 2 or N x 3 point array, got {points.shape}.")

    if plane == "xy":
        axes = (0, 1)
    elif plane == "xz":
        if points.shape[1] < 3:
            raise ValueError("The xz projection requires three-dimensional point coordinates.")
        axes = (0, 2)
    else:
        if points.shape[1] < 3:
            raise ValueError("The yz projection requires three-dimensional point coordinates.")
        axes = (1, 2)

    return np.asarray(points[:, axes], dtype=float)


def load_mesh_frame(entry: TimeStep, field_name: str, plane: str) -> MeshFrame:
    dataset = read_vtk_dataset(entry.path)
    surface = triangulated_surface(dataset, entry.path)
    return MeshFrame(
        coordinates=project_points(dataset_points(surface, entry.path), plane),
        triangles=surface_triangles(surface, entry.path),
        values=np.asarray(point_data_array(surface, field_name, entry.path), dtype=float),
    )


def load_particle_frame(
    entry: TimeStep,
    field_name: str,
    plane: str,
) -> ParticleFrame:
    dataset = read_vtk_dataset(entry.path)
    coordinates = project_points(dataset_points(dataset, entry.path), plane)
    particle_types = np.asarray(
        point_data_array(dataset, field_name, entry.path), dtype=np.int64
    )
    return ParticleFrame(coordinates=coordinates, particle_types=particle_types)

def closest_timestep(
    entries: list[TimeStep],
    target_time: float,
    absolute_tolerance: float | None,
) -> TimeStep:
    times = np.fromiter((entry.time for entry in entries), dtype=float)
    index = int(np.argmin(np.abs(times - target_time)))
    entry = entries[index]
    distance = abs(entry.time - target_time)
    tolerance = (
        absolute_tolerance
        if absolute_tolerance is not None
        else 1.0e-9 * max(1.0, abs(target_time))
    )

    if distance > tolerance:
        raise ValueError(
            "No particle timestep matches mesh time "
            f"{target_time:.16g}; nearest particle time is {entry.time:.16g}. "
            "Use --time-tolerance to permit a larger mismatch."
        )

    return entry


def zoomed_limits(coordinates: np.ndarray, zoom: float) -> tuple[tuple[float, float], tuple[float, float]]:
    if not 0.0 < zoom:
        raise ValueError("--zoom must be greater than zero.")

    minimum = np.nanmin(coordinates, axis=0)
    maximum = np.nanmax(coordinates, axis=0)
    center = 0.5 * (minimum + maximum)
    half_range = 0.5 * (maximum - minimum) * zoom

    for axis in range(2):
        if not np.isfinite(half_range[axis]) or half_range[axis] <= 0.0:
            half_range[axis] = 0.5

    return (
        (center[0] - half_range[0], center[0] + half_range[0]),
        (center[1] - half_range[1], center[1] + half_range[1]),
    )


def configure_axes(ax: plt.Axes, plane: str) -> None:
    labels = {
        "xy": ("x ($\\mu$m)", "y ($\\mu$m)"),
        "xz": ("x ($\\mu$m)", "z ($\\mu$m)"),
        "yz": ("y ($\\mu$m)", "z ($\\mu$m)"),
    }
    x_label, y_label = labels[plane]

    ax.set_facecolor("black")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(x_label, color="white")
    ax.set_ylabel(y_label, color="white")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_color("white")


def create_animation(args: argparse.Namespace) -> tuple[plt.Figure, animation.FuncAnimation]:
    scene_dir = args.scene_dir.expanduser().resolve()
    mesh_pvd = resolve_path(scene_dir, args.mesh_pvd).resolve()
    particles_pvd = resolve_path(scene_dir, args.particles_pvd).resolve()

    all_mesh_steps = read_pvd(mesh_pvd)
    particle_steps = read_pvd(particles_pvd)
    mesh_iterations, mesh_steps = select_iterations(
        all_mesh_steps,
        args.min_iteration,
        args.max_iteration,
        args.iterstep,
    )

    LOGGER.info(
        "Animating %d iteration(s): %d through %d with step %d",
        len(mesh_steps),
        mesh_iterations[0],
        mesh_iterations[-1],
        args.iterstep,
    )

    first_mesh = load_mesh_frame(mesh_steps[0], args.field, args.plane)
    first_particle_entry = closest_timestep(
        particle_steps, mesh_steps[0].time, args.time_tolerance
    )
    first_particles = load_particle_frame(
        first_particle_entry, args.particle_field, args.plane
    )

    field_norm = Normalize(vmin=args.vmin, vmax=args.vmax, clip=True)
    particle_cmap = ListedColormap(
        np.array(
            [
                [0.25, 0.0, 1.0, 1.0],
                [1.0, 0.0, 0.0, 1.0],
            ]
        )
    )
    particle_norm = BoundaryNorm([-0.5, 0.5, 1.5], particle_cmap.N, clip=True)

    figure, ax = plt.subplots(figsize=(11.0, 7.0), constrained_layout=True)
    figure.patch.set_facecolor("black")
    configure_axes(ax, args.plane)

    x_limits, y_limits = zoomed_limits(first_mesh.coordinates, args.zoom)
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)

    triangulation = mtri.Triangulation(
        first_mesh.coordinates[:, 0],
        first_mesh.coordinates[:, 1],
        first_mesh.triangles,
    )
    mesh_artist = ax.tripcolor(
        triangulation,
        first_mesh.values,
        shading="gouraud",
        cmap="jet",
        norm=field_norm,
        zorder=1,
    )
    particle_artist = ax.scatter(
        first_particles.coordinates[:, 0],
        first_particles.coordinates[:, 1],
        c=first_particles.particle_types,
        cmap=particle_cmap,
        norm=particle_norm,
        s=args.particle_size,
        linewidths=0.0,
        zorder=2,
    )
    time_artist = ax.text(
        0.5,
        1.1,
        f"time: {mesh_steps[0].time:.3e} s",
        transform=ax.transAxes,
        ha="center",
        va="top",
        color="white",
        fontsize=18,
        zorder=3,
    )

    field_mappable = ScalarMappable(norm=field_norm, cmap="jet")
    field_mappable.set_array([])
    field_bar = figure.colorbar(field_mappable, ax=ax, location="bottom", pad=0.02)
    field_bar.set_label(args.field, color="white")
    field_bar.ax.tick_params(colors="white")
    field_bar.outline.set_edgecolor("white")

    particle_mappable = ScalarMappable(norm=particle_norm, cmap=particle_cmap)
    particle_mappable.set_array([])
    # particle_bar = figure.colorbar(
    #     particle_mappable,
    #     ax=ax,
    #     location="right",
    #     pad=0.10,
    #     ticks=[0, 1],
    # )
    # particle_bar.set_label(args.particle_field, color="white")
    # particle_bar.set_ticklabels(["electron", "hole"])
    # particle_bar.ax.tick_params(colors="white")
    # particle_bar.outline.set_edgecolor("white")

    state: dict[str, object] = {
        "mesh_artist": mesh_artist,
        "coordinates": first_mesh.coordinates,
        "triangles": first_mesh.triangles,
    }

    def update(frame_index: int) -> tuple[object, ...]:
        mesh_entry = mesh_steps[frame_index]
        particle_entry = closest_timestep(
            particle_steps, mesh_entry.time, args.time_tolerance
        )
        mesh_frame = load_mesh_frame(mesh_entry, args.field, args.plane)
        particle_frame = load_particle_frame(
            particle_entry, args.particle_field, args.plane
        )

        geometry_is_unchanged = (
            np.array_equal(state["coordinates"], mesh_frame.coordinates)
            and np.array_equal(state["triangles"], mesh_frame.triangles)
        )

        current_mesh_artist = state["mesh_artist"]
        if geometry_is_unchanged:
            current_mesh_artist.set_array(mesh_frame.values)
        else:
            current_mesh_artist.remove()
            current_triangulation = mtri.Triangulation(
                mesh_frame.coordinates[:, 0],
                mesh_frame.coordinates[:, 1],
                mesh_frame.triangles,
            )
            current_mesh_artist = ax.tripcolor(
                current_triangulation,
                mesh_frame.values,
                shading="gouraud",
                cmap="jet",
                norm=field_norm,
                zorder=1,
            )
            state["mesh_artist"] = current_mesh_artist
            state["coordinates"] = mesh_frame.coordinates
            state["triangles"] = mesh_frame.triangles

        particle_artist.set_offsets(particle_frame.coordinates)
        particle_artist.set_array(particle_frame.particle_types.astype(float, copy=False))
        time_artist.set_text(f"time: {mesh_entry.time:.3e} s")

        LOGGER.info(
            "Rendered frame %d/%d from iteration %d at t=%.6g s",
            frame_index + 1,
            len(mesh_steps),
            mesh_iterations[frame_index],
            mesh_entry.time,
        )
        return current_mesh_artist, particle_artist, time_artist

    frame_interval = args.interval if args.interval is not None else 1000.0 / args.fps
    scene_animation = animation.FuncAnimation(
        figure,
        update,
        frames=len(mesh_steps),
        interval=frame_interval,
        repeat=False,
        blit=False,
        cache_frame_data=False,
    )
    return figure, scene_animation


def save_animation(
    figure: plt.Figure,
    scene_animation: animation.FuncAnimation,
    output_path: Path,
    fps: float,
    dpi: int,
) -> None:
    output_path = output_path.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    suffix = output_path.suffix.lower()

    if suffix == ".mp4":
        if not animation.writers.is_available("ffmpeg"):
            raise RuntimeError(
                "Matplotlib cannot find ffmpeg. Install ffmpeg or choose a .gif output."
            )
        writer = animation.FFMpegWriter(
            fps=fps,
            codec="h264",
            bitrate=-1,
            extra_args=["-pix_fmt", "yuv420p"],
            metadata={"title": "ADMC simulation"},
        )
    elif suffix == ".gif":
        if not animation.writers.is_available("pillow"):
            raise RuntimeError(
                "The Pillow animation writer is unavailable. Install pillow or choose .mp4."
            )
        writer = animation.PillowWriter(fps=fps)
    else:
        raise ValueError("--output must have a .mp4 or .gif extension.")

    LOGGER.info("Saving %s", output_path)
    scene_animation.save(output_path, writer=writer, dpi=dpi)
    LOGGER.info("Saved %s", output_path)


def main() -> int:
    args = parse_arguments()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s: %(message)s",
    )

    try:
        figure, scene_animation = create_animation(args)
        save_animation(figure, scene_animation, args.output, args.fps, args.dpi)
        if args.show:
            plt.show()
        else:
            plt.close(figure)
    except (FileNotFoundError, KeyError, RuntimeError, ValueError) as error:
        LOGGER.error("%s", error)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())