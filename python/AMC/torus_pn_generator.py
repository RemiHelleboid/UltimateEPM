#!/usr/bin/env python3

import argparse
import importlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d


FULL_CIRCLE = 2.0 * np.pi
gmsh = None


def require_gmsh():
    global gmsh

    if gmsh is None:
        try:
            gmsh = importlib.import_module("gmsh")
        except ImportError as exc:
            raise RuntimeError(
                "The gmsh Python module is required to generate the mesh. "
                "Install it with: python -m pip install gmsh"
            ) from exc

    return gmsh


def smooth_log(x, values, smoothing_length, floor):
    values = np.clip(values, floor, None)

    if smoothing_length <= 0.0:
        return values

    dx = x[1] - x[0]
    sigma = smoothing_length / dx

    log_values = np.log10(values)
    smoothed_log_values = gaussian_filter1d(log_values, sigma)

    return 10.0**smoothed_log_values


def raw_acceptor_concentration(
    s,
    length,
    contact_doping_length,
    contact_doping_level,
    peak_p_level,
    diffusion_length,
):
    value = peak_p_level * np.exp(
        -(((s + contact_doping_length) - length) ** 2) / diffusion_length**2
    )

    value = np.where(s > length - contact_doping_length, contact_doping_level, value)

    return value


def raw_donor_concentration(
    s,
    contact_doping_length,
    contact_doping_level,
    peak_n_level,
    diffusion_length,
):
    value = peak_n_level * np.exp(
        -((s - contact_doping_length) ** 2) / diffusion_length**2
    )

    value = np.where(s < contact_doping_length, contact_doping_level, value)

    return value


def compute_raw_doping_profile(
    s,
    length,
    contact_doping_length,
    contact_doping_level,
    peak_n_level,
    peak_p_level,
    diffusion_length,
    min_doping,
    max_doping,
):
    donor = raw_donor_concentration(
        s=s,
        contact_doping_length=contact_doping_length,
        contact_doping_level=contact_doping_level,
        peak_n_level=peak_n_level,
        diffusion_length=diffusion_length,
    )

    acceptor = raw_acceptor_concentration(
        s=s,
        length=length,
        contact_doping_length=contact_doping_length,
        contact_doping_level=contact_doping_level,
        peak_p_level=peak_p_level,
        diffusion_length=diffusion_length,
    )

    donor = np.clip(donor, min_doping, max_doping)
    acceptor = np.clip(acceptor, min_doping, max_doping)

    return acceptor, donor


def compute_doping_profile(
    s,
    length,
    contact_doping_length,
    contact_doping_level,
    peak_n_level,
    peak_p_level,
    diffusion_length,
    min_doping,
    max_doping,
    apply_smoothing,
    smoothing_length,
):
    s = np.asarray(s, dtype=float)
    s = np.clip(s, 0.0, length)

    if not apply_smoothing:
        acceptor, donor = compute_raw_doping_profile(
            s=s,
            length=length,
            contact_doping_length=contact_doping_length,
            contact_doping_level=contact_doping_level,
            peak_n_level=peak_n_level,
            peak_p_level=peak_p_level,
            diffusion_length=diffusion_length,
            min_doping=min_doping,
            max_doping=max_doping,
        )

        return acceptor, donor, donor - acceptor

    n_grid = max(1000, s.size)
    s_grid = np.linspace(0.0, length, n_grid)

    acceptor_grid, donor_grid = compute_raw_doping_profile(
        s=s_grid,
        length=length,
        contact_doping_length=contact_doping_length,
        contact_doping_level=contact_doping_level,
        peak_n_level=peak_n_level,
        peak_p_level=peak_p_level,
        diffusion_length=diffusion_length,
        min_doping=min_doping,
        max_doping=max_doping,
    )

    donor_grid = smooth_log(s_grid, donor_grid, smoothing_length, min_doping)
    acceptor_grid = smooth_log(s_grid, acceptor_grid, smoothing_length, min_doping)

    donor = np.interp(s, s_grid, donor_grid)
    acceptor = np.interp(s, s_grid, acceptor_grid)

    return acceptor, donor, donor - acceptor


def add_physical_group(dim, tags, name):
    physical_tag = gmsh.model.addPhysicalGroup(dim, tags)
    gmsh.model.setPhysicalName(dim, physical_tag, name)

    return physical_tag


def get_entity_center(dim, tag):
    return np.array(gmsh.model.occ.getCenterOfMass(dim, tag), dtype=float)


def find_surface_near_point(point, scale, label):
    point = np.asarray(point, dtype=float)
    tolerance = 1.0e-7 * max(scale, 1.0)

    matches = []
    distances = []

    for dim, tag in gmsh.model.getEntities(2):
        center = get_entity_center(dim, tag)
        distance = float(np.linalg.norm(center - point))
        distances.append((tag, distance, center))

        if distance <= tolerance:
            matches.append(tag)

    if len(matches) == 1:
        return matches

    distances.sort(key=lambda item: item[1])
    nearest = "; ".join(
        f"tag={tag}, distance={distance:.6e}, center=({center[0]:.6e}, {center[1]:.6e}, {center[2]:.6e})"
        for tag, distance, center in distances[:5]
    )

    raise RuntimeError(
        f"Expected one {label} surface near ({point[0]:.6e}, {point[1]:.6e}, {point[2]:.6e}), "
        f"found {len(matches)}. Nearest surfaces: {nearest}"
    )


def create_half_torus_geometry(major_radius, minor_radius, arc_angle):
    bulk_tag = gmsh.model.occ.addTorus(
        0.0,
        0.0,
        0.0,
        major_radius,
        minor_radius,
        -1,
        arc_angle,
    )

    gmsh.model.occ.synchronize()

    add_physical_group(3, [bulk_tag], "Silicon_1")

    start_point = np.array([major_radius, 0.0, 0.0])
    end_point = np.array([
        major_radius * np.cos(arc_angle),
        major_radius * np.sin(arc_angle),
        0.0,
    ])

    scale = major_radius + minor_radius

    cathode_surfaces = find_surface_near_point(start_point, scale, "cathode")
    anode_surfaces = find_surface_near_point(end_point, scale, "anode")

    add_physical_group(2, cathode_surfaces, "cathode")
    add_physical_group(2, anode_surfaces, "anode")


def arc_coordinate_from_torus_nodes(node_coordinates, major_radius, arc_angle):
    x = node_coordinates[0::3]
    y = node_coordinates[1::3]

    theta = np.mod(np.arctan2(y, x), FULL_CIRCLE)
    theta = np.clip(theta, 0.0, arc_angle)

    return major_radius * theta


def add_node_view(model_name, mesh_file, view_name, node_tags, values):
    view_tag = gmsh.view.add(view_name)

    gmsh.view.addHomogeneousModelData(
        view_tag,
        0,
        model_name,
        "NodeData",
        node_tags,
        values,
    )

    gmsh.view.write(view_tag, str(mesh_file), True)


def generate_mesh(
    mesh_file,
    major_radius,
    minor_radius,
    arc_angle,
    h_min,
    h_max,
    contact_doping_length,
    contact_doping_level,
    peak_n_level,
    peak_p_level,
    diffusion_length,
    min_doping,
    max_doping,
    apply_smoothing,
    smoothing_length,
):
    require_gmsh()

    model_name = "Half_Torus_PN_Junction"
    length = major_radius * arc_angle

    gmsh.initialize()

    try:
        gmsh.model.add(model_name)

        gmsh.option.setNumber("Mesh.MeshSizeMin", h_min)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h_max)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)

        create_half_torus_geometry(
            major_radius=major_radius,
            minor_radius=minor_radius,
            arc_angle=arc_angle,
        )

        gmsh.model.mesh.generate(3)

        node_tags, node_coordinates, _ = gmsh.model.mesh.getNodes()
        s = arc_coordinate_from_torus_nodes(
            node_coordinates=node_coordinates,
            major_radius=major_radius,
            arc_angle=arc_angle,
        )

        acceptor, donor, net_doping = compute_doping_profile(
            s=s,
            length=length,
            contact_doping_length=contact_doping_length,
            contact_doping_level=contact_doping_level,
            peak_n_level=peak_n_level,
            peak_p_level=peak_p_level,
            diffusion_length=diffusion_length,
            min_doping=min_doping,
            max_doping=max_doping,
            apply_smoothing=apply_smoothing,
            smoothing_length=smoothing_length,
        )

        gmsh.write(str(mesh_file))

        add_node_view(model_name, mesh_file, "DonorConcentration", node_tags, donor)
        add_node_view(model_name, mesh_file, "AcceptorConcentration", node_tags, acceptor)
        add_node_view(model_name, mesh_file, "DopingConcentration", node_tags, net_doping)

    finally:
        gmsh.finalize()


def export_profile_plot(
    output_prefix,
    length,
    contact_doping_length,
    contact_doping_level,
    peak_n_level,
    peak_p_level,
    diffusion_length,
    min_doping,
    max_doping,
    apply_smoothing,
    smoothing_length,
    show_plot,
):
    s = np.linspace(0.0, length, 1000)

    acceptor, donor, net_doping = compute_doping_profile(
        s=s,
        length=length,
        contact_doping_length=contact_doping_length,
        contact_doping_level=contact_doping_level,
        peak_n_level=peak_n_level,
        peak_p_level=peak_p_level,
        diffusion_length=diffusion_length,
        min_doping=min_doping,
        max_doping=max_doping,
        apply_smoothing=apply_smoothing,
        smoothing_length=smoothing_length,
    )

    profile_file = output_prefix.with_suffix(".profile.csv")

    np.savetxt(
        profile_file,
        np.column_stack([s, acceptor, donor, net_doping]),
        delimiter=",",
        header="arc_length,AcceptorConcentration,DonorConcentration,DopingConcentration",
        comments="",
    )

    fig, ax = plt.subplots()

    ax.plot(s, acceptor, label="Acceptor")
    ax.plot(s, donor, label="Donor")
    ax.plot(s, np.abs(net_doping), label="|Net doping|", linestyle="--")

    ax.set_xlabel("arc length s (µm)")
    ax.set_ylabel("Concentration (cm⁻³)")
    ax.set_yscale("log")
    ax.set_ylim(min_doping, 2.0 * max(peak_n_level, peak_p_level, contact_doping_level))
    ax.set_title(output_prefix.name)
    ax.grid(True, which="both")
    ax.legend()

    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".doping_profile.png"), dpi=200)

    if show_plot:
        plt.show()

    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a 3D half-torus PN junction mesh with curved doping profiles."
    )

    parser.add_argument(
        "mesh_name",
        help="Output mesh filename, for example half_torus_diode.msh.",
    )

    parser.add_argument(
        "--major-radius",
        type=float,
        default=1.0,
        help="Major radius of the torus centerline in µm.",
    )

    parser.add_argument(
        "--minor-radius",
        type=float,
        default=0.15,
        help="Minor radius of the torus tube in µm.",
    )

    parser.add_argument(
        "--angle-deg",
        type=float,
        default=180.0,
        help="Angular span of the torus in degrees. Use 180 for a half-torus.",
    )

    parser.add_argument(
        "--hmin",
        type=float,
        default=0.01,
        help="Minimum mesh size in µm.",
    )

    parser.add_argument(
        "--hmax",
        type=float,
        default=0.03,
        help="Maximum mesh size in µm.",
    )

    parser.add_argument(
        "--contact-doping-length",
        type=float,
        default=0.25,
        help="Highly doped contact length measured along the torus arc in µm.",
    )

    parser.add_argument(
        "--contact-doping-level",
        type=float,
        default=2.0e18,
        help="Contact doping concentration in cm^-3.",
    )

    parser.add_argument(
        "--n-level",
        type=float,
        default=2.0e18,
        help="Peak donor concentration in cm^-3.",
    )

    parser.add_argument(
        "--p-level",
        type=float,
        default=2.0e18,
        help="Peak acceptor concentration in cm^-3.",
    )

    parser.add_argument(
        "--diffusion-length",
        type=float,
        default=0.12,
        help="Gaussian diffusion length along the torus arc in µm.",
    )

    parser.add_argument(
        "--min-doping",
        type=float,
        default=1.0e11,
        help="Minimum clipped doping concentration in cm^-3.",
    )

    parser.add_argument(
        "--max-doping",
        type=float,
        default=1.0e19,
        help="Maximum clipped doping concentration in cm^-3.",
    )

    parser.add_argument(
        "--smooth",
        action="store_true",
        help="Apply log-scale Gaussian smoothing to the doping profile.",
    )

    parser.add_argument(
        "--smoothing-length",
        type=float,
        default=0.01,
        help="Smoothing length along the torus arc in µm.",
    )

    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Do not create the doping profile plot.",
    )

    parser.add_argument(
        "--show",
        action="store_true",
        help="Show the plot interactively.",
    )

    return parser.parse_args()


def validate_args(args):
    if args.major_radius <= 0.0:
        raise ValueError("--major-radius must be positive.")

    if args.minor_radius <= 0.0:
        raise ValueError("--minor-radius must be positive.")

    if args.minor_radius >= args.major_radius:
        raise ValueError("--minor-radius must be smaller than --major-radius.")

    if args.angle_deg <= 0.0:
        raise ValueError("--angle-deg must be positive.")

    if args.angle_deg >= 360.0:
        raise ValueError("--angle-deg must be smaller than 360.")

    if args.hmin <= 0.0:
        raise ValueError("--hmin must be positive.")

    if args.hmax <= 0.0:
        raise ValueError("--hmax must be positive.")

    if args.hmin > args.hmax:
        raise ValueError("--hmin cannot be larger than --hmax.")

    length = args.major_radius * np.deg2rad(args.angle_deg)

    if args.contact_doping_length <= 0.0:
        raise ValueError("--contact-doping-length must be positive.")

    if args.contact_doping_length >= length:
        raise ValueError("--contact-doping-length must be smaller than the torus arc length.")

    if args.contact_doping_level <= 0.0:
        raise ValueError("--contact-doping-level must be positive.")

    if args.n_level <= 0.0:
        raise ValueError("--n-level must be positive.")

    if args.p_level <= 0.0:
        raise ValueError("--p-level must be positive.")

    if args.diffusion_length <= 0.0:
        raise ValueError("--diffusion-length must be positive.")

    if args.min_doping <= 0.0:
        raise ValueError("--min-doping must be positive.")

    if args.max_doping <= 0.0:
        raise ValueError("--max-doping must be positive.")

    if args.min_doping > args.max_doping:
        raise ValueError("--min-doping cannot be larger than --max-doping.")

    if args.smoothing_length <= 0.0:
        raise ValueError("--smoothing-length must be positive.")


def main():
    args = parse_args()
    validate_args(args)

    mesh_file = Path(args.mesh_name)

    if mesh_file.suffix != ".msh":
        mesh_file = mesh_file.with_suffix(".msh")

    output_prefix = mesh_file.with_suffix("")
    arc_angle = np.deg2rad(args.angle_deg)
    length = args.major_radius * arc_angle

    if not args.no_plot:
        export_profile_plot(
            output_prefix=output_prefix,
            length=length,
            contact_doping_length=args.contact_doping_length,
            contact_doping_level=args.contact_doping_level,
            peak_n_level=args.n_level,
            peak_p_level=args.p_level,
            diffusion_length=args.diffusion_length,
            min_doping=args.min_doping,
            max_doping=args.max_doping,
            apply_smoothing=args.smooth,
            smoothing_length=args.smoothing_length,
            show_plot=args.show,
        )

    generate_mesh(
        mesh_file=mesh_file,
        major_radius=args.major_radius,
        minor_radius=args.minor_radius,
        arc_angle=arc_angle,
        h_min=args.hmin,
        h_max=args.hmax,
        contact_doping_length=args.contact_doping_length,
        contact_doping_level=args.contact_doping_level,
        peak_n_level=args.n_level,
        peak_p_level=args.p_level,
        diffusion_length=args.diffusion_length,
        min_doping=args.min_doping,
        max_doping=args.max_doping,
        apply_smoothing=args.smooth,
        smoothing_length=args.smoothing_length,
    )

    print(f"Wrote {mesh_file}")

    if not args.no_plot:
        print(f"Wrote {output_prefix.with_suffix('.profile.csv')}")
        print(f"Wrote {output_prefix.with_suffix('.doping_profile.png')}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError) as exc:
        raise SystemExit(f"error: {exc}") from None