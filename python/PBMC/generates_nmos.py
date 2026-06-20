#!/usr/bin/env python3

import argparse
from pathlib import Path

import gmsh
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d


SILICON_PHYSICAL_NAME = "Silicon_1"
OXIDE_PHYSICAL_NAME = "Oxide_1"


def add_physical_group(dim, tags, name):
    physical_tag = gmsh.model.addPhysicalGroup(dim, tags)
    gmsh.model.setPhysicalName(dim, physical_tag, name)

    return physical_tag


def smooth_step(value, edge, smoothing_length, direction):
    if smoothing_length <= 0.0:
        if direction > 0:
            return np.where(value >= edge, 1.0, 0.0)

        return np.where(value <= edge, 1.0, 0.0)

    argument = (value - edge) / smoothing_length

    if direction > 0:
        return 0.5 * (1.0 + np.tanh(argument))

    return 0.5 * (1.0 - np.tanh(argument))


def box_window(value, lower, upper, smoothing_length):
    left = smooth_step(value, lower, smoothing_length, 1)
    right = smooth_step(value, upper, smoothing_length, -1)

    return left * right


def smooth_log(x, values, smoothing_length):
    values = np.clip(values, 1.0e11, None)

    if smoothing_length <= 0.0:
        return values

    dx = x[1] - x[0]
    sigma = smoothing_length / dx

    log_values = np.log10(values)
    smoothed_log_values = gaussian_filter1d(log_values, sigma)

    return 10.0**smoothed_log_values


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

    acceptor = np.full_like(x, p_body_level, dtype=float)
    donor = np.full_like(x, background_doping, dtype=float)

    source_mask = smooth_step(x, source_length, lateral_smoothing_length, -1)
    drain_mask = smooth_step(x, length - drain_length, lateral_smoothing_length, 1)
    surface_mask = smooth_step(
        y,
        silicon_thickness - junction_depth,
        vertical_smoothing_length,
        1,
    )

    n_implant_mask = np.clip((source_mask + drain_mask) * surface_mask, 0.0, 1.0)

    donor = donor + n_plus_level * n_implant_mask

    donor = np.clip(donor, min_doping, max_doping)
    acceptor = np.clip(acceptor, min_doping, max_doping)

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


def get_boundary_entities(surface_tag):
    boundary = gmsh.model.getBoundary([(2, surface_tag)], oriented=False, recursive=False)

    return [tag for dim, tag in boundary if dim == 1]


def curve_center(curve_tag):
    return gmsh.model.occ.getCenterOfMass(1, curve_tag)


def curve_bbox(curve_tag):
    xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(1, curve_tag)

    return xmin, ymin, xmax, ymax


def collect_curves(curve_tags, predicate):
    return [tag for tag in curve_tags if predicate(tag)]


def point_bbox(point_tag):
    xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(0, point_tag)

    return xmin, ymin, xmax, ymax


def set_point_mesh_sizes(
    length,
    silicon_thickness,
    oxide_thickness,
    source_length,
    drain_length,
    gate_start,
    gate_end,
    h_min,
    h_max,
    gate_edge_refinement_length,
):
    oxide_top = silicon_thickness + oxide_thickness
    tolerance = 1.0e-8 * max(length, silicon_thickness + oxide_thickness, 1.0)

    points = [tag for _, tag in gmsh.model.getEntities(0)]

    if points:
        gmsh.model.mesh.setSize([(0, tag) for tag in points], h_max)

    refined_points = []

    for point_tag in points:
        xmin, ymin, xmax, ymax = point_bbox(point_tag)
        x = 0.5 * (xmin + xmax)
        y = 0.5 * (ymin + ymax)

        on_silicon_oxide_interface = abs(y - silicon_thickness) <= tolerance
        on_gate_top = abs(y - oxide_top) <= tolerance and gate_start - tolerance <= x <= gate_end + tolerance
        near_source_junction = abs(x - source_length) <= gate_edge_refinement_length and y >= silicon_thickness - gate_edge_refinement_length
        near_drain_junction = abs(x - (length - drain_length)) <= gate_edge_refinement_length and y >= silicon_thickness - gate_edge_refinement_length
        near_gate_edge = (
            abs(x - gate_start) <= gate_edge_refinement_length
            or abs(x - gate_end) <= gate_edge_refinement_length
        ) and y >= silicon_thickness - gate_edge_refinement_length

        if on_silicon_oxide_interface or on_gate_top or near_source_junction or near_drain_junction or near_gate_edge:
            refined_points.append(point_tag)

    if refined_points:
        gmsh.model.mesh.setSize([(0, tag) for tag in refined_points], h_min)


def create_geometry(
    length,
    silicon_thickness,
    oxide_thickness,
    source_length,
    drain_length,
    gate_length,
):
    gate_start = 0.5 * (length - gate_length)
    gate_end = gate_start + gate_length

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

    silicon_surfaces = []
    oxide_surfaces = []

    tolerance = 1.0e-6 * max(length, silicon_thickness + oxide_thickness, 1.0)

    for _, tag in gmsh.model.getEntities(2):
        _, ymin, _, _, ymax, _ = gmsh.model.getBoundingBox(2, tag)

        if ymax <= silicon_thickness + tolerance:
            silicon_surfaces.append(tag)
        elif ymin >= silicon_thickness - tolerance:
            oxide_surfaces.append(tag)

    if not oxide_surfaces:
        for _, tag in gmsh.model.getEntities(2):
            _, ymin, _, _, ymax, _ = gmsh.model.getBoundingBox(2, tag)
            center_y = 0.5 * (ymin + ymax)

            if center_y > silicon_thickness:
                oxide_surfaces.append(tag)
                if tag in silicon_surfaces:
                    silicon_surfaces.remove(tag)

    if not silicon_surfaces:
        raise RuntimeError("No silicon surface found.")

    if not oxide_surfaces:
        raise RuntimeError("No oxide surface found.")

    add_physical_group(2, silicon_surfaces, SILICON_PHYSICAL_NAME)
    add_physical_group(2, oxide_surfaces, OXIDE_PHYSICAL_NAME)

    all_curves = [tag for _, tag in gmsh.model.getEntities(1)]

    def has_horizontal_extent(curve_tag):
        xmin, ymin, xmax, ymax = curve_bbox(curve_tag)
        return xmax - xmin > tolerance and ymax - ymin <= 10.0 * tolerance

    def has_vertical_extent(curve_tag):
        xmin, ymin, xmax, ymax = curve_bbox(curve_tag)
        return ymax - ymin > tolerance and xmax - xmin <= 10.0 * tolerance

    def is_horizontal_at(y_target):
        def predicate(curve_tag):
            xmin, ymin, xmax, ymax = curve_bbox(curve_tag)
            return has_horizontal_extent(curve_tag) and ymin <= y_target + tolerance and ymax >= y_target - tolerance

        return predicate

    def is_vertical_at(x_target):
        def predicate(curve_tag):
            xmin, ymin, xmax, ymax = curve_bbox(curve_tag)
            return has_vertical_extent(curve_tag) and xmin <= x_target + tolerance and xmax >= x_target - tolerance

        return predicate

    def interval_overlaps(a_min, a_max, b_min, b_max):
        return min(a_max, b_max) >= max(a_min, b_min) - tolerance

    top_silicon_curves = collect_curves(all_curves, is_horizontal_at(silicon_thickness))
    bottom_silicon_curves = collect_curves(all_curves, is_horizontal_at(0.0))
    top_oxide_curves = collect_curves(all_curves, is_horizontal_at(silicon_thickness + oxide_thickness))
    left_silicon_curves = collect_curves(all_curves, is_vertical_at(0.0))
    right_silicon_curves = collect_curves(all_curves, is_vertical_at(length))

    source_top = []
    drain_top = []

    for curve_tag in top_silicon_curves:
        xmin, _, xmax, _ = curve_bbox(curve_tag)
        center_x = 0.5 * (xmin + xmax)

        if interval_overlaps(xmin, xmax, 0.0, source_length) and center_x <= source_length + tolerance:
            source_top.append(curve_tag)
        elif interval_overlaps(xmin, xmax, length - drain_length, length) and center_x >= length - drain_length - tolerance:
            drain_top.append(curve_tag)

    source_curves = sorted(set(source_top + left_silicon_curves))
    drain_curves = sorted(set(drain_top + right_silicon_curves))

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

    interface_curves = []
    for curve_tag in top_silicon_curves:
        xmin, _, xmax, _ = curve_bbox(curve_tag)
        if interval_overlaps(xmin, xmax, gate_start, gate_end):
            interface_curves.append(curve_tag)

    gate_edge_curves = []
    for curve_tag in all_curves:
        xmin, ymin, xmax, ymax = curve_bbox(curve_tag)
        center_x = 0.5 * (xmin + xmax)
        if (
            has_vertical_extent(curve_tag)
            and (abs(center_x - gate_start) <= tolerance or abs(center_x - gate_end) <= tolerance)
            and ymin <= silicon_thickness + tolerance
            and ymax >= silicon_thickness - tolerance
        ):
            gate_edge_curves.append(curve_tag)

    return {
        "gate_start": gate_start,
        "gate_end": gate_end,
        "silicon_surfaces": silicon_surfaces,
        "oxide_surfaces": oxide_surfaces,
        "interface_curves": sorted(set(interface_curves)),
        "gate_edge_curves": sorted(set(gate_edge_curves)),
    }


def add_box_refinement_field(x_min, x_max, y_min, y_max, h_min, h_max):
    field = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(field, "VIn", h_min)
    gmsh.model.mesh.field.setNumber(field, "VOut", h_max)
    gmsh.model.mesh.field.setNumber(field, "XMin", x_min)
    gmsh.model.mesh.field.setNumber(field, "XMax", x_max)
    gmsh.model.mesh.field.setNumber(field, "YMin", y_min)
    gmsh.model.mesh.field.setNumber(field, "YMax", y_max)

    return field


def add_curve_refinement_field(curve_tags, h_min, h_max, distance_min, distance_max):
    if not curve_tags:
        return None

    distance_field = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(distance_field, "CurvesList", curve_tags)
    gmsh.model.mesh.field.setNumber(distance_field, "Sampling", 100)

    threshold_field = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(threshold_field, "InField", distance_field)
    gmsh.model.mesh.field.setNumber(threshold_field, "SizeMin", h_min)
    gmsh.model.mesh.field.setNumber(threshold_field, "SizeMax", h_max)
    gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", distance_min)
    gmsh.model.mesh.field.setNumber(threshold_field, "DistMax", distance_max)

    return threshold_field


def set_mesh_fields(
    length,
    silicon_thickness,
    oxide_thickness,
    source_length,
    drain_length,
    gate_start,
    gate_end,
    interface_curves,
    gate_edge_curves,
    h_min,
    h_max,
    interface_refinement_length,
    junction_refinement_length,
    gate_edge_refinement_length,
):
    oxide_top = silicon_thickness + oxide_thickness
    source_channel_x = source_length
    drain_channel_x = length - drain_length

    fields = []

    interface_field = add_curve_refinement_field(
        curve_tags=interface_curves,
        h_min=h_min,
        h_max=h_max,
        distance_min=0.0,
        distance_max=interface_refinement_length,
    )
    if interface_field is not None:
        fields.append(interface_field)

    gate_edge_field = add_curve_refinement_field(
        curve_tags=gate_edge_curves,
        h_min=h_min,
        h_max=h_max,
        distance_min=0.0,
        distance_max=gate_edge_refinement_length,
    )
    if gate_edge_field is not None:
        fields.append(gate_edge_field)

    fields.append(
        add_box_refinement_field(
            x_min=gate_start,
            x_max=gate_end,
            y_min=silicon_thickness - interface_refinement_length,
            y_max=oxide_top,
            h_min=h_min,
            h_max=h_max,
        )
    )

    for x0 in [source_channel_x, drain_channel_x]:
        fields.append(
            add_box_refinement_field(
                x_min=x0 - junction_refinement_length,
                x_max=x0 + junction_refinement_length,
                y_min=silicon_thickness - 2.0 * junction_refinement_length,
                y_max=oxide_top,
                h_min=h_min,
                h_max=h_max,
            )
        )

    for x0 in [gate_start, gate_end]:
        fields.append(
            add_box_refinement_field(
                x_min=x0 - gate_edge_refinement_length,
                x_max=x0 + gate_edge_refinement_length,
                y_min=silicon_thickness - 2.0 * gate_edge_refinement_length,
                y_max=oxide_top,
                h_min=h_min,
                h_max=h_max,
            )
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
    gate_length,
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

        geometry = create_geometry(
            length=length,
            silicon_thickness=silicon_thickness,
            oxide_thickness=oxide_thickness,
            source_length=source_length,
            drain_length=drain_length,
            gate_length=gate_length,
        )

        set_point_mesh_sizes(
            length=length,
            silicon_thickness=silicon_thickness,
            oxide_thickness=oxide_thickness,
            source_length=source_length,
            drain_length=drain_length,
            gate_start=geometry["gate_start"],
            gate_end=geometry["gate_end"],
            h_min=h_min,
            h_max=h_max,
            gate_edge_refinement_length=gate_edge_refinement_length,
        )

        set_mesh_fields(
            length=length,
            silicon_thickness=silicon_thickness,
            oxide_thickness=oxide_thickness,
            source_length=source_length,
            drain_length=drain_length,
            gate_start=geometry["gate_start"],
            gate_end=geometry["gate_end"],
            interface_curves=geometry["interface_curves"],
            gate_edge_curves=geometry["gate_edge_curves"],
            h_min=h_min,
            h_max=h_max,
            interface_refinement_length=interface_refinement_length,
            junction_refinement_length=junction_refinement_length,
            gate_edge_refinement_length=gate_edge_refinement_length,
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
        description="Generate a simple realistic 2D nMOS mesh with silicon and oxide physical groups."
    )

    parser.add_argument(
        "mesh_name",
        help="Output mesh filename, for example nmos.msh.",
    )

    parser.add_argument(
        "--length",
        type=float,
        default=1.0,
        help="Total device length in µm.",
    )

    parser.add_argument(
        "--silicon-thickness",
        type=float,
        default=0.2,
        help="Silicon thickness in µm.",
    )

    parser.add_argument(
        "--oxide-thickness",
        type=float,
        default=0.01,
        help="Gate oxide thickness in µm.",
    )

    parser.add_argument(
        "--source-length",
        type=float,
        default=0.2,
        help="Source contact/implant length in µm.",
    )

    parser.add_argument(
        "--drain-length",
        type=float,
        default=0.2,
        help="Drain contact/implant length in µm.",
    )

    parser.add_argument(
        "--gate-length",
        type=float,
        default=0.6,
        help="Gate length in µm. The gate is centered between source and drain by default.",
    )

    parser.add_argument(
        "--hmin",
        type=float,
        default=0.002,
        help="Minimum mesh size in µm.",
    )

    parser.add_argument(
        "--hmax",
        type=float,
        default=0.02,
        help="Maximum mesh size in µm.",
    )

    parser.add_argument(
        "--interface-refinement-length",
        type=float,
        default=0.025,
        help="Refinement distance around the silicon/oxide interface in µm.",
    )

    parser.add_argument(
        "--junction-refinement-length",
        type=float,
        default=0.035,
        help="Refinement half-width around source/channel and drain/channel junctions in µm.",
    )

    parser.add_argument(
        "--gate-edge-refinement-length",
        type=float,
        default=0.025,
        help="Refinement half-width around gate edges in µm.",
    )

    parser.add_argument(
        "--n-plus-level",
        type=float,
        default=1.0e19,
        help="Peak N+ source/drain donor concentration in cm^-3.",
    )

    parser.add_argument(
        "--p-body-level",
        type=float,
        default=1.0e16,
        help="P-body acceptor concentration in cm^-3.",
    )

    parser.add_argument(
        "--background-doping",
        type=float,
        default=1.0e11,
        help="Minimum background donor concentration in cm^-3.",
    )

    parser.add_argument(
        "--junction-depth",
        type=float,
        default=0.05,
        help="Approximate source/drain junction depth from the silicon surface in µm.",
    )

    parser.add_argument(
        "--lateral-smoothing-length",
        type=float,
        default=0.015,
        help="Lateral smoothing length of source/drain junctions in µm.",
    )

    parser.add_argument(
        "--vertical-smoothing-length",
        type=float,
        default=0.01,
        help="Vertical smoothing length of source/drain implants in µm.",
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
        default=1.0e20,
        help="Maximum clipped doping concentration in cm^-3.",
    )

    parser.add_argument(
        "--T-source",
        dest="temperature_source",
        type=float,
        default=300.0,
        help="Temperature at x = 0 in K.",
    )

    parser.add_argument(
        "--T-drain",
        dest="temperature_drain",
        type=float,
        default=300.0,
        help="Temperature at x = length in K.",
    )

    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Do not create the channel doping plot.",
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

    if args.silicon_thickness <= 0.0:
        raise ValueError("--silicon-thickness must be positive.")

    if args.oxide_thickness <= 0.0:
        raise ValueError("--oxide-thickness must be positive.")

    if args.source_length <= 0.0:
        raise ValueError("--source-length must be positive.")

    if args.drain_length <= 0.0:
        raise ValueError("--drain-length must be positive.")

    if args.gate_length <= 0.0:
        raise ValueError("--gate-length must be positive.")

    if args.source_length + args.drain_length >= args.length:
        raise ValueError("--source-length + --drain-length must be smaller than --length.")

    if args.gate_length > args.length - args.source_length - args.drain_length:
        raise ValueError("--gate-length must fit between source and drain regions.")

    expected_gate_length = args.length - args.source_length - args.drain_length
    if abs(args.gate_length - expected_gate_length) > 1.0e-12:
        gate_start = 0.5 * (args.length - args.gate_length)
        gate_end = gate_start + args.gate_length

        if gate_start < args.source_length or gate_end > args.length - args.drain_length:
            raise ValueError("Centered gate overlaps source or drain. Adjust --gate-length or contact lengths.")

    if args.hmin <= 0.0:
        raise ValueError("--hmin must be positive.")

    if args.hmax <= 0.0:
        raise ValueError("--hmax must be positive.")

    if args.hmin > args.hmax:
        raise ValueError("--hmin cannot be larger than --hmax.")

    if args.interface_refinement_length <= 0.0:
        raise ValueError("--interface-refinement-length must be positive.")

    if args.junction_refinement_length <= 0.0:
        raise ValueError("--junction-refinement-length must be positive.")

    if args.gate_edge_refinement_length <= 0.0:
        raise ValueError("--gate-edge-refinement-length must be positive.")

    if args.n_plus_level <= 0.0:
        raise ValueError("--n-plus-level must be positive.")

    if args.p_body_level <= 0.0:
        raise ValueError("--p-body-level must be positive.")

    if args.background_doping <= 0.0:
        raise ValueError("--background-doping must be positive.")

    if args.junction_depth <= 0.0:
        raise ValueError("--junction-depth must be positive.")

    if args.junction_depth >= args.silicon_thickness:
        raise ValueError("--junction-depth must be smaller than --silicon-thickness.")

    if args.lateral_smoothing_length <= 0.0:
        raise ValueError("--lateral-smoothing-length must be positive.")

    if args.vertical_smoothing_length <= 0.0:
        raise ValueError("--vertical-smoothing-length must be positive.")

    if args.min_doping <= 0.0:
        raise ValueError("--min-doping must be positive.")

    if args.max_doping <= 0.0:
        raise ValueError("--max-doping must be positive.")

    if args.min_doping > args.max_doping:
        raise ValueError("--min-doping cannot be larger than --max-doping.")

    if args.temperature_source <= 0.0:
        raise ValueError("--T-source must be positive.")

    if args.temperature_drain <= 0.0:
        raise ValueError("--T-drain must be positive.")


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
        gate_length=args.gate_length,
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