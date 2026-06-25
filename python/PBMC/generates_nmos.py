#!/usr/bin/env python3

import argparse
from pathlib import Path

import gmsh
import matplotlib.pyplot as plt
import numpy as np


SILICON_PHYSICAL_NAME = "Silicon_1"
OXIDE_PHYSICAL_NAME = "SiO2_1"


def add_physical_group(dim, tags, name):
    physical_tag = gmsh.model.addPhysicalGroup(dim, sorted(set(tags)))
    gmsh.model.setPhysicalName(dim, physical_tag, name)
    return physical_tag


def smooth_step(value, edge, smoothing_length, direction):
    value = np.asarray(value, dtype=float)

    if smoothing_length <= 0.0:
        if direction > 0:
            return np.where(value >= edge, 1.0, 0.0)
        return np.where(value <= edge, 1.0, 0.0)

    argument = (value - edge) / smoothing_length

    if direction > 0:
        return 0.5 * (1.0 + np.tanh(argument))
    return 0.5 * (1.0 - np.tanh(argument))


def curve_bbox(curve_tag):
    xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(1, curve_tag)
    return xmin, ymin, xmax, ymax


def interval_overlaps(a_min, a_max, b_min, b_max, tolerance):
    return min(a_max, b_max) >= max(a_min, b_min) - tolerance


def compute_doping_profile(
    x,
    y,
    length,
    silicon_thickness,
    source_length,
    drain_length,
    n_plus_level,
    p_body_level,
    background_doping,
    junction_depth,
    lateral_smoothing_length,
    vertical_smoothing_length,
    min_doping,
    max_doping,
):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    silicon_mask = y <= silicon_thickness + 1.0e-12 * max(length, silicon_thickness, 1.0)

    donor = np.zeros_like(x, dtype=float)
    acceptor = np.zeros_like(x, dtype=float)

    donor[silicon_mask] = background_doping
    acceptor[silicon_mask] = p_body_level

    source_mask = smooth_step(x, source_length, lateral_smoothing_length, -1)
    drain_mask = smooth_step(x, length - drain_length, lateral_smoothing_length, 1)
    surface_mask = smooth_step(
        y,
        silicon_thickness - junction_depth,
        vertical_smoothing_length,
        1,
    )

    implant_mask = np.clip((source_mask + drain_mask) * surface_mask, 0.0, 1.0)
    donor += n_plus_level * implant_mask * silicon_mask

    donor = np.where(silicon_mask, np.clip(donor, min_doping, max_doping), 0.0)
    acceptor = np.where(silicon_mask, np.clip(acceptor, min_doping, max_doping), 0.0)

    return acceptor, donor, donor - acceptor


def compute_temperature_profile(x, length, temperature_source, temperature_drain):
    x = np.asarray(x, dtype=float)
    return temperature_source + (temperature_drain - temperature_source) * x / length


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


def create_geometry(length, silicon_thickness, oxide_thickness, source_length, drain_length):
    gate_start = source_length
    gate_end = length - drain_length
    gate_length = gate_end - gate_start

    silicon = gmsh.model.occ.addRectangle(
        0.0,
        0.0,
        0.0,
        length,
        silicon_thickness,
    )
    oxide = gmsh.model.occ.addRectangle(
        gate_start,
        silicon_thickness,
        0.0,
        gate_length,
        oxide_thickness,
    )

    gmsh.model.occ.fragment([(2, silicon)], [(2, oxide)])
    gmsh.model.occ.synchronize()

    tolerance = 1.0e-6 * max(length, silicon_thickness + oxide_thickness, 1.0)
    oxide_top = silicon_thickness + oxide_thickness

    silicon_surfaces = []
    oxide_surfaces = []

    for _, tag in gmsh.model.getEntities(2):
        _, ymin, _, _, ymax, _ = gmsh.model.getBoundingBox(2, tag)
        center_y = 0.5 * (ymin + ymax)

        if center_y < silicon_thickness:
            silicon_surfaces.append(tag)
        else:
            oxide_surfaces.append(tag)

    if not silicon_surfaces:
        raise RuntimeError("No silicon surface found.")
    if not oxide_surfaces:
        raise RuntimeError("No SiO2 surface found.")

    add_physical_group(2, silicon_surfaces, SILICON_PHYSICAL_NAME)
    add_physical_group(2, oxide_surfaces, OXIDE_PHYSICAL_NAME)

    def is_horizontal_at(curve_tag, y_target):
        xmin, ymin, xmax, ymax = curve_bbox(curve_tag)
        center_y = 0.5 * (ymin + ymax)
        return (xmax - xmin > tolerance) and (abs(center_y - y_target) <= tolerance)

    all_curves = [tag for _, tag in gmsh.model.getEntities(1)]

    top_silicon_curves = [tag for tag in all_curves if is_horizontal_at(tag, silicon_thickness)]
    bottom_silicon_curves = [tag for tag in all_curves if is_horizontal_at(tag, 0.0)]
    top_oxide_curves = [tag for tag in all_curves if is_horizontal_at(tag, oxide_top)]

    source_curves = []
    drain_curves = []

    for curve_tag in top_silicon_curves:
        xmin, _, xmax, _ = curve_bbox(curve_tag)
        center_x = 0.5 * (xmin + xmax)

        if interval_overlaps(xmin, xmax, 0.0, source_length, tolerance) and center_x <= source_length + tolerance:
            source_curves.append(curve_tag)
        elif interval_overlaps(xmin, xmax, length - drain_length, length, tolerance) and center_x >= length - drain_length - tolerance:
            drain_curves.append(curve_tag)

    if not source_curves:
        raise RuntimeError("No source boundary found.")
    if not drain_curves:
        raise RuntimeError("No drain boundary found.")
    if not top_oxide_curves:
        raise RuntimeError("No gate boundary found.")
    if not bottom_silicon_curves:
        raise RuntimeError("No body boundary found.")

    add_physical_group(1, source_curves, "source")
    add_physical_group(1, drain_curves, "drain")
    add_physical_group(1, top_oxide_curves, "gate")
    add_physical_group(1, bottom_silicon_curves, "body")

    return gate_start, gate_end


def set_mesh_fields(
    length,
    silicon_thickness,
    oxide_thickness,
    source_length,
    drain_length,
    h_min,
    h_max,
    interface_refinement_length,
    junction_refinement_length,
    gate_edge_refinement_length,
    junction_depth,
):
    gate_start = source_length
    gate_end = length - drain_length
    oxide_top = silicon_thickness + oxide_thickness

    fields = []

    def add_box(x_min, x_max, y_min, y_max):
        field = gmsh.model.mesh.field.add("Box")
        gmsh.model.mesh.field.setNumber(field, "VIn", h_min)
        gmsh.model.mesh.field.setNumber(field, "VOut", h_max)
        gmsh.model.mesh.field.setNumber(field, "XMin", max(0.0, x_min))
        gmsh.model.mesh.field.setNumber(field, "XMax", min(length, x_max))
        gmsh.model.mesh.field.setNumber(field, "YMin", max(0.0, y_min))
        gmsh.model.mesh.field.setNumber(field, "YMax", min(oxide_top, y_max))
        fields.append(field)

    add_box(
        0.0,
        length,
        silicon_thickness - interface_refinement_length,
        oxide_top,
    )

    for x0 in (gate_start, gate_end):
        add_box(
            x0 - gate_edge_refinement_length,
            x0 + gate_edge_refinement_length,
            silicon_thickness - gate_edge_refinement_length,
            oxide_top,
        )

    for x0 in (source_length, length - drain_length):
        add_box(
            x0 - junction_refinement_length,
            x0 + junction_refinement_length,
            silicon_thickness - junction_depth - junction_refinement_length,
            silicon_thickness,
        )

    implant_transition_y = silicon_thickness - junction_depth
    implant_refinement_height = max(junction_refinement_length, 3.0 * h_min)

    add_box(
        0.0,
        source_length + junction_refinement_length,
        implant_transition_y - implant_refinement_height,
        silicon_thickness,
    )
    add_box(
        length - drain_length - junction_refinement_length,
        length,
        implant_transition_y - implant_refinement_height,
        silicon_thickness,
    )

    minimum_field = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(minimum_field, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(minimum_field)


def generate_mesh(
    mesh_file,
    length,
    silicon_thickness,
    oxide_thickness,
    source_length,
    drain_length,
    h_min,
    h_max,
    interface_refinement_length,
    junction_refinement_length,
    gate_edge_refinement_length,
    n_plus_level,
    p_body_level,
    background_doping,
    junction_depth,
    lateral_smoothing_length,
    vertical_smoothing_length,
    min_doping,
    max_doping,
    temperature_source,
    temperature_drain,
):
    model_name = "Simple_NMOS_2D"

    gmsh.initialize()

    try:
        gmsh.model.add(model_name)

        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.MeshSizeMin", h_min)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h_max)

        create_geometry(
            length=length,
            silicon_thickness=silicon_thickness,
            oxide_thickness=oxide_thickness,
            source_length=source_length,
            drain_length=drain_length,
        )

        set_mesh_fields(
            length=length,
            silicon_thickness=silicon_thickness,
            oxide_thickness=oxide_thickness,
            source_length=source_length,
            drain_length=drain_length,
            h_min=h_min,
            h_max=h_max,
            interface_refinement_length=interface_refinement_length,
            junction_refinement_length=junction_refinement_length,
            gate_edge_refinement_length=gate_edge_refinement_length,
            junction_depth=junction_depth,
        )

        gmsh.model.mesh.generate(2)

        node_tags, node_coordinates, _ = gmsh.model.mesh.getNodes()
        x = node_coordinates[0::3]
        y = node_coordinates[1::3]

        acceptor, donor, net_doping = compute_doping_profile(
            x=x,
            y=y,
            length=length,
            silicon_thickness=silicon_thickness,
            source_length=source_length,
            drain_length=drain_length,
            n_plus_level=n_plus_level,
            p_body_level=p_body_level,
            background_doping=background_doping,
            junction_depth=junction_depth,
            lateral_smoothing_length=lateral_smoothing_length,
            vertical_smoothing_length=vertical_smoothing_length,
            min_doping=min_doping,
            max_doping=max_doping,
        )
        temperature = compute_temperature_profile(
            x=x,
            length=length,
            temperature_source=temperature_source,
            temperature_drain=temperature_drain,
        )

        gmsh.write(str(mesh_file))

        add_node_view(model_name, mesh_file, "DonorConcentration", node_tags, donor)
        add_node_view(model_name, mesh_file, "AcceptorConcentration", node_tags, acceptor)
        add_node_view(model_name, mesh_file, "DopingConcentration", node_tags, net_doping)
        add_node_view(model_name, mesh_file, "Temperature", node_tags, temperature)

    finally:
        gmsh.finalize()


def export_profile_plot(
    output_prefix,
    length,
    silicon_thickness,
    source_length,
    drain_length,
    n_plus_level,
    p_body_level,
    background_doping,
    junction_depth,
    lateral_smoothing_length,
    vertical_smoothing_length,
    min_doping,
    max_doping,
    show_plot,
):
    x = np.linspace(0.0, length, 1200)
    y_surface = np.full_like(x, silicon_thickness)
    y_mid = np.full_like(x, 0.5 * silicon_thickness)

    acceptor_surface, donor_surface, net_surface = compute_doping_profile(
        x=x,
        y=y_surface,
        length=length,
        silicon_thickness=silicon_thickness,
        source_length=source_length,
        drain_length=drain_length,
        n_plus_level=n_plus_level,
        p_body_level=p_body_level,
        background_doping=background_doping,
        junction_depth=junction_depth,
        lateral_smoothing_length=lateral_smoothing_length,
        vertical_smoothing_length=vertical_smoothing_length,
        min_doping=min_doping,
        max_doping=max_doping,
    )
    acceptor_mid, donor_mid, net_mid = compute_doping_profile(
        x=x,
        y=y_mid,
        length=length,
        silicon_thickness=silicon_thickness,
        source_length=source_length,
        drain_length=drain_length,
        n_plus_level=n_plus_level,
        p_body_level=p_body_level,
        background_doping=background_doping,
        junction_depth=junction_depth,
        lateral_smoothing_length=lateral_smoothing_length,
        vertical_smoothing_length=vertical_smoothing_length,
        min_doping=min_doping,
        max_doping=max_doping,
    )

    profile_file = output_prefix.with_suffix(".channel_profile.csv")
    np.savetxt(
        profile_file,
        np.column_stack(
            [
                x,
                acceptor_surface,
                donor_surface,
                net_surface,
                acceptor_mid,
                donor_mid,
                net_mid,
            ]
        ),
        delimiter=",",
        header=(
            "x,SurfaceAcceptorConcentration,SurfaceDonorConcentration,"
            "SurfaceDopingConcentration,MidAcceptorConcentration,"
            "MidDonorConcentration,MidDopingConcentration"
        ),
        comments="",
    )

    fig, ax = plt.subplots()
    ax.plot(x, np.clip(donor_surface, min_doping, None), label="Donor at Si surface")
    ax.plot(x, np.clip(acceptor_surface, min_doping, None), label="Acceptor at Si surface")
    ax.plot(x, np.clip(np.abs(net_surface), min_doping, None), label="|Net| at Si surface", linestyle="--")
    ax.plot(x, np.clip(np.abs(net_mid), min_doping, None), label="|Net| at mid Si", linestyle=":")

    ax.set_xlabel("x (µm)")
    ax.set_ylabel("Concentration (cm⁻³)")
    ax.set_yscale("log")
    ax.set_ylim(min_doping, 2.0 * max(n_plus_level, p_body_level, background_doping))
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
        description="Generate a simple 2D nMOS mesh with Silicon_1 and SiO2_1 physical groups."
    )

    parser.add_argument("mesh_name", help="Output mesh filename, for example nmos.msh.")

    parser.add_argument("--length", type=float, default=1.0, help="Total device length in µm.")
    parser.add_argument("--silicon-thickness", type=float, default=0.2, help="Silicon thickness in µm.")
    parser.add_argument("--oxide-thickness", type=float, default=0.01, help="Gate SiO2 thickness in µm.")
    parser.add_argument("--source-length", type=float, default=0.2, help="Source contact and implant length in µm.")
    parser.add_argument("--drain-length", type=float, default=0.2, help="Drain contact and implant length in µm.")

    parser.add_argument("--hmin", type=float, default=0.002, help="Minimum mesh size in µm.")
    parser.add_argument("--hmax", type=float, default=0.02, help="Maximum mesh size in µm.")
    parser.add_argument(
        "--interface-refinement-length",
        type=float,
        default=0.025,
        help="Vertical refinement distance around the Si/SiO2 interface in µm.",
    )
    parser.add_argument(
        "--junction-refinement-length",
        type=float,
        default=0.025,
        help="Refinement half-width around source/channel and drain/channel junctions in µm.",
    )
    parser.add_argument(
        "--gate-edge-refinement-length",
        type=float,
        default=0.015,
        help="Refinement half-width around the two gate edges in µm.",
    )

    parser.add_argument("--n-plus-level", type=float, default=1.0e19, help="Peak N+ source/drain donor concentration in cm^-3.")
    parser.add_argument("--p-body-level", type=float, default=1.0e16, help="P-body acceptor concentration in cm^-3.")
    parser.add_argument("--background-doping", type=float, default=1.0e11, help="Background donor concentration in silicon in cm^-3.")
    parser.add_argument("--junction-depth", type=float, default=0.05, help="Source/drain junction depth from the silicon surface in µm.")
    parser.add_argument("--lateral-smoothing-length", type=float, default=0.015, help="Lateral source/drain junction smoothing length in µm.")
    parser.add_argument("--vertical-smoothing-length", type=float, default=0.01, help="Vertical source/drain implant smoothing length in µm.")
    parser.add_argument("--min-doping", type=float, default=1.0e11, help="Minimum clipped silicon doping concentration in cm^-3.")
    parser.add_argument("--max-doping", type=float, default=1.0e20, help="Maximum clipped silicon doping concentration in cm^-3.")

    parser.add_argument("--T-source", dest="temperature_source", type=float, default=300.0, help="Temperature at x = 0 in K.")
    parser.add_argument("--T-drain", dest="temperature_drain", type=float, default=300.0, help="Temperature at x = length in K.")

    parser.add_argument("--no-plot", action="store_true", help="Do not create the channel doping plot.")
    parser.add_argument("--show", action="store_true", help="Show the plot interactively.")

    return parser.parse_args()


def validate_args(args):
    positive_fields = {
        "--length": args.length,
        "--silicon-thickness": args.silicon_thickness,
        "--oxide-thickness": args.oxide_thickness,
        "--source-length": args.source_length,
        "--drain-length": args.drain_length,
        "--hmin": args.hmin,
        "--hmax": args.hmax,
        "--interface-refinement-length": args.interface_refinement_length,
        "--junction-refinement-length": args.junction_refinement_length,
        "--gate-edge-refinement-length": args.gate_edge_refinement_length,
        "--n-plus-level": args.n_plus_level,
        "--p-body-level": args.p_body_level,
        "--background-doping": args.background_doping,
        "--junction-depth": args.junction_depth,
        "--lateral-smoothing-length": args.lateral_smoothing_length,
        "--vertical-smoothing-length": args.vertical_smoothing_length,
        "--min-doping": args.min_doping,
        "--max-doping": args.max_doping,
        "--T-source": args.temperature_source,
        "--T-drain": args.temperature_drain,
    }

    for name, value in positive_fields.items():
        if value <= 0.0:
            raise ValueError(f"{name} must be positive.")

    if args.hmin > args.hmax:
        raise ValueError("--hmin cannot be larger than --hmax.")

    if args.source_length + args.drain_length >= args.length:
        raise ValueError("--source-length + --drain-length must be smaller than --length.")

    if args.junction_depth >= args.silicon_thickness:
        raise ValueError("--junction-depth must be smaller than --silicon-thickness.")

    if args.min_doping > args.max_doping:
        raise ValueError("--min-doping cannot be larger than --max-doping.")


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
            silicon_thickness=args.silicon_thickness,
            source_length=args.source_length,
            drain_length=args.drain_length,
            n_plus_level=args.n_plus_level,
            p_body_level=args.p_body_level,
            background_doping=args.background_doping,
            junction_depth=args.junction_depth,
            lateral_smoothing_length=args.lateral_smoothing_length,
            vertical_smoothing_length=args.vertical_smoothing_length,
            min_doping=args.min_doping,
            max_doping=args.max_doping,
            show_plot=args.show,
        )

    generate_mesh(
        mesh_file=mesh_file,
        length=args.length,
        silicon_thickness=args.silicon_thickness,
        oxide_thickness=args.oxide_thickness,
        source_length=args.source_length,
        drain_length=args.drain_length,
        h_min=args.hmin,
        h_max=args.hmax,
        interface_refinement_length=args.interface_refinement_length,
        junction_refinement_length=args.junction_refinement_length,
        gate_edge_refinement_length=args.gate_edge_refinement_length,
        n_plus_level=args.n_plus_level,
        p_body_level=args.p_body_level,
        background_doping=args.background_doping,
        junction_depth=args.junction_depth,
        lateral_smoothing_length=args.lateral_smoothing_length,
        vertical_smoothing_length=args.vertical_smoothing_length,
        min_doping=args.min_doping,
        max_doping=args.max_doping,
        temperature_source=args.temperature_source,
        temperature_drain=args.temperature_drain,
    )

    print(f"Wrote {mesh_file}")
    if not args.no_plot:
        print(f"Wrote {output_prefix.with_suffix('.channel_profile.csv')}")
        print(f"Wrote {output_prefix.with_suffix('.doping_profile.png')}")


if __name__ == "__main__":
    main()