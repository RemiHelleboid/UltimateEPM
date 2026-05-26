#!/usr/bin/env python3

import argparse
from pathlib import Path

import gmsh
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d


def smooth_log(x, values, smoothing_length):
    values = np.clip(values, 1.0e11, None)

    dx = x[1] - x[0]
    sigma = smoothing_length / dx

    log_values = np.log10(values)
    smoothed_log_values = gaussian_filter1d(log_values, sigma)

    return 10.0**smoothed_log_values


def acceptor_concentration(
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

    if x > length - contact_doping_length:
        value = contact_doping_level

    return value


def donor_concentration(
    x,
    contact_doping_length,
    contact_doping_level,
    peak_n_level,
    diffusion_length,
):
    value = peak_n_level * np.exp(
        -((x - contact_doping_length) ** 2) / diffusion_length**2
    )

    if x < contact_doping_length:
        value = contact_doping_level

    return value


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
    donor = np.array(
        [
            donor_concentration(
                xi,
                contact_doping_length,
                contact_doping_level,
                peak_n_level,
                diffusion_length,
            )
            for xi in x
        ]
    )

    acceptor = np.array(
        [
            acceptor_concentration(
                xi,
                length,
                contact_doping_length,
                contact_doping_level,
                peak_p_level,
                diffusion_length,
            )
            for xi in x
        ]
    )

    donor = np.clip(donor, min_doping, max_doping)
    acceptor = np.clip(acceptor, min_doping, max_doping)

    if apply_smoothing:
        donor = smooth_log(x, donor, smoothing_length)
        acceptor = smooth_log(x, acceptor, smoothing_length)

    net_doping = donor - acceptor

    return acceptor, donor, net_doping


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
    model_name = "PN_Diode"

    gmsh.initialize()
    try:
        gmsh.model.add(model_name)

        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.MeshSizeMin", h_min)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h_max)

        box_tag = gmsh.model.occ.addBox(0.0, 0.0, 0.0, length, width, width)
        gmsh.model.occ.synchronize()

        gmsh.model.addPhysicalGroup(3, [box_tag], -1, "Silicon_1")

        eps = 1.0e-6

        anode_surfaces = gmsh.model.getEntitiesInBoundingBox(
            -eps,
            -eps,
            -eps,
            eps,
            width + eps,
            width + eps,
            dim=2,
        )

        cathode_surfaces = gmsh.model.getEntitiesInBoundingBox(
            length - eps,
            -eps,
            -eps,
            length + eps,
            width + eps,
            width + eps,
            dim=2,
        )

        if len(anode_surfaces) != 1:
            raise RuntimeError(f"Expected one anode surface, found {len(anode_surfaces)}.")

        if len(cathode_surfaces) != 1:
            raise RuntimeError(
                f"Expected one cathode surface, found {len(cathode_surfaces)}."
            )

        gmsh.model.addPhysicalGroup(2, [anode_surfaces[0][1]], -1, "anode")
        gmsh.model.addPhysicalGroup(2, [cathode_surfaces[0][1]], -1, "cathode")

        gmsh.model.mesh.generate(3)

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
        description="Generate a simple 3D PN diode mesh with doping profiles."
    )

    parser.add_argument("mesh_name", help="Output mesh filename, for example diode.msh.")

    parser.add_argument("--length", type=float, default=2.0, help="Device length in µm.")
    parser.add_argument("--width", type=float, default=0.3, help="Device width in µm.")

    parser.add_argument("--hmin", type=float, default=0.01, help="Minimum mesh size in µm.")
    parser.add_argument("--hmax", type=float, default=0.015, help="Maximum mesh size in µm.")

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


def main():
    args = parse_args()

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