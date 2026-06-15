#!/usr/bin/env python3

import argparse
from pathlib import Path

import gmsh
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d


def smooth_log(x, values, smoothing_length):
    values = np.clip(values, 1.0e11, None)

    if smoothing_length <= 0.0:
        return values

    dx = x[1] - x[0]
    sigma = smoothing_length / dx

    log_values = np.log10(values)
    smoothed_log_values = gaussian_filter1d(log_values, sigma)

    return 10.0**smoothed_log_values


def raw_acceptor_concentration(
    x,
    length,
    contact_doping_length,
    contact_doping_level,
    peak_p_level,
    diffusion_length,
):
    value = peak_p_level * np.exp(
        -(((x + contact_doping_length) - length) ** 2) / diffusion_length**2
    )

    value = np.where(x > length - contact_doping_length, contact_doping_level, value)

    return value


def raw_donor_concentration(
    x,
    contact_doping_length,
    contact_doping_level,
    peak_n_level,
    diffusion_length,
):
    value = peak_n_level * np.exp(
        -((x - contact_doping_length) ** 2) / diffusion_length**2
    )

    value = np.where(x < contact_doping_length, contact_doping_level, value)

    return value


def compute_raw_doping_profile(
    x,
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
        x=x,
        contact_doping_length=contact_doping_length,
        contact_doping_level=contact_doping_level,
        peak_n_level=peak_n_level,
        diffusion_length=diffusion_length,
    )

    acceptor = raw_acceptor_concentration(
        x=x,
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
    x,
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
    x = np.asarray(x, dtype=float)

    if not apply_smoothing:
        acceptor, donor = compute_raw_doping_profile(
            x=x,
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

    n_grid = max(1000, x.size)
    x_grid = np.linspace(0.0, length, n_grid)

    acceptor_grid, donor_grid = compute_raw_doping_profile(
        x=x_grid,
        length=length,
        contact_doping_length=contact_doping_length,
        contact_doping_level=contact_doping_level,
        peak_n_level=peak_n_level,
        peak_p_level=peak_p_level,
        diffusion_length=diffusion_length,
        min_doping=min_doping,
        max_doping=max_doping,
    )

    donor_grid = smooth_log(x_grid, donor_grid, smoothing_length)
    acceptor_grid = smooth_log(x_grid, acceptor_grid, smoothing_length)

    donor = np.interp(x, x_grid, donor_grid)
    acceptor = np.interp(x, x_grid, acceptor_grid)

    net_doping = donor - acceptor

    return acceptor, donor, net_doping


def add_physical_group(dim, tags, name):
    physical_tag = gmsh.model.addPhysicalGroup(dim, tags)
    gmsh.model.setPhysicalName(dim, physical_tag, name)

    return physical_tag


def find_entities_at_x(entity_dim, x_target, length, width):
    tolerance = 1.0e-8 * max(length, width, 1.0)

    matching_tags = []

    for dim, tag in gmsh.model.getEntities(entity_dim):
        center_x, _, _ = gmsh.model.occ.getCenterOfMass(dim, tag)

        if abs(center_x - x_target) <= tolerance:
            matching_tags.append(tag)

    return matching_tags


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


def create_geometry(dimension, length, width):
    if dimension == 3:
        bulk_dim = 3
        contact_dim = 2

        bulk_tag = gmsh.model.occ.addBox(
            0.0,
            0.0,
            0.0,
            length,
            width,
            width,
        )

    elif dimension == 2:
        bulk_dim = 2
        contact_dim = 1

        bulk_tag = gmsh.model.occ.addRectangle(
            0.0,
            0.0,
            0.0,
            length,
            width,
        )

    else:
        raise ValueError("dimension must be 2 or 3.")

    gmsh.model.occ.synchronize()

    add_physical_group(bulk_dim, [bulk_tag], "Silicon_1")

    cathode_entities = find_entities_at_x(
        entity_dim=contact_dim,
        x_target=0.0,
        length=length,
        width=width,
    )

    anode_entities = find_entities_at_x(
        entity_dim=contact_dim,
        x_target=length,
        length=length,
        width=width,
    )

    if len(anode_entities) != 1:
        raise RuntimeError(f"Expected one anode entity, found {len(anode_entities)}.")

    if len(cathode_entities) != 1:
        raise RuntimeError(f"Expected one cathode entity, found {len(cathode_entities)}.")

    add_physical_group(contact_dim, cathode_entities, "cathode")
    add_physical_group(contact_dim, anode_entities, "anode")


def generate_mesh(
    mesh_file,
    dimension,
    length,
    width,
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
    model_name = f"PN_Diode_{dimension}D"

    gmsh.initialize()

    try:
        gmsh.model.add(model_name)

        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.MeshSizeMin", h_min)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h_max)

        create_geometry(
            dimension=dimension,
            length=length,
            width=width,
        )

        gmsh.model.mesh.generate(dimension)

        node_tags, node_coordinates, _ = gmsh.model.mesh.getNodes()

        x = node_coordinates[0::3]

        acceptor, donor, net_doping = compute_doping_profile(
            x=x,
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
    x = np.linspace(0.0, length, 1000)

    acceptor, donor, net_doping = compute_doping_profile(
        x=x,
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
        np.column_stack([x, acceptor, donor, net_doping]),
        delimiter=",",
        header="x,AcceptorConcentration,DonorConcentration,DopingConcentration",
        comments="",
    )

    fig, ax = plt.subplots()

    ax.plot(x, acceptor, label="Acceptor")
    ax.plot(x, donor, label="Donor")
    ax.plot(x, np.abs(net_doping), label="|Net doping|", linestyle="--")

    ax.set_xlabel("x (µm)")
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
        description="Generate a simple 2D or 3D PN diode mesh with doping profiles."
    )

    parser.add_argument(
        "mesh_name",
        help="Output mesh filename, for example diode.msh.",
    )

    parser.add_argument(
        "--dimension",
        type=int,
        choices=[2, 3],
        default=3,
        help="Mesh dimension: 2 or 3.",
    )

    parser.add_argument(
        "--length",
        type=float,
        default=2.0,
        help="Device length in µm.",
    )

    parser.add_argument(
        "--width",
        type=float,
        default=0.3,
        help="Device width in µm. In 3D, the device is width × width in transverse directions.",
    )

    parser.add_argument(
        "--hmin",
        type=float,
        default=0.005,
        help="Minimum mesh size in µm.",
    )

    parser.add_argument(
        "--hmax",
        type=float,
        default=0.01,
        help="Maximum mesh size in µm.",
    )

    parser.add_argument(
        "--contact-doping-length",
        type=float,
        default=0.70,
        help="Length of the highly doped contact regions in µm.",
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
        default=0.15,
        help="Gaussian diffusion length in µm.",
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
        help="Smoothing length in µm.",
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
    if args.length <= 0.0:
        raise ValueError("--length must be positive.")

    if args.width <= 0.0:
        raise ValueError("--width must be positive.")

    if args.hmin <= 0.0:
        raise ValueError("--hmin must be positive.")

    if args.hmax <= 0.0:
        raise ValueError("--hmax must be positive.")

    if args.hmin > args.hmax:
        raise ValueError("--hmin cannot be larger than --hmax.")

    if args.contact_doping_length <= 0.0:
        raise ValueError("--contact-doping-length must be positive.")

    if args.contact_doping_length >= args.length:
        raise ValueError("--contact-doping-length must be smaller than --length.")

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

    if not args.no_plot:
        export_profile_plot(
            output_prefix=output_prefix,
            length=args.length,
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
        dimension=args.dimension,
        length=args.length,
        width=args.width,
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
    main()