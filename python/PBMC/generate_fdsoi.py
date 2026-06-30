
import argparse
from pathlib import Path

import gmsh
import matplotlib.pyplot as plt
import numpy as np


SILICON_FILM_PHYSICAL_NAME = "Silicon_1"
SUBSTRATE_PHYSICAL_NAME = "Silicon_2"
FRONT_OXIDE_PHYSICAL_NAME = "SiO2_1"
BOX_PHYSICAL_NAME = "SiO2_2"


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
    substrate_thickness,
    box_thickness,
    film_thickness,
    source_length,
    drain_length,
    n_plus_level,
    channel_acceptor_level,
    substrate_acceptor_level,
    background_donor_level,
    source_drain_depth,
    lateral_smoothing_length,
    vertical_smoothing_length,
    min_doping,
    max_doping,
):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    box_bottom = substrate_thickness
    box_top = substrate_thickness + box_thickness
    film_bottom = box_top
    film_top = film_bottom + film_thickness

    tolerance = 1.0e-12 * max(length, film_top, 1.0)

    film_mask = (y >= film_bottom - tolerance) & (y <= film_top + tolerance)
    substrate_mask = y <= substrate_thickness + tolerance

    donor = np.zeros_like(x, dtype=float)
    acceptor = np.zeros_like(x, dtype=float)

    donor[film_mask] = background_donor_level
    acceptor[film_mask] = channel_acceptor_level

    donor[substrate_mask] = background_donor_level
    acceptor[substrate_mask] = substrate_acceptor_level

    source_mask = smooth_step(x, source_length, lateral_smoothing_length, -1)
    drain_mask = smooth_step(x, length - drain_length, lateral_smoothing_length, 1)

    if source_drain_depth >= film_thickness:
        vertical_mask = np.ones_like(x, dtype=float)
    else:
        implant_bottom = film_top - source_drain_depth
        vertical_mask = smooth_step(y, implant_bottom, vertical_smoothing_length, 1)

    implant_mask = np.clip((source_mask + drain_mask) * vertical_mask, 0.0, 1.0) * film_mask
    donor += n_plus_level * implant_mask

    semiconductor_mask = film_mask | substrate_mask
    donor = np.where(semiconductor_mask, np.clip(donor, min_doping, max_doping), 0.0)
    acceptor = np.where(semiconductor_mask, np.clip(acceptor, min_doping, max_doping), 0.0)

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


def create_geometry(
    length,
    substrate_thickness,
    box_thickness,
    film_thickness,
    front_oxide_thickness,
    source_length,
    drain_length,
):
    gate_start = source_length
    gate_end = length - drain_length
    gate_length = gate_end - gate_start

    box_bottom = substrate_thickness
    box_top = substrate_thickness + box_thickness
    film_bottom = box_top
    film_top = film_bottom + film_thickness
    oxide_top = film_top + front_oxide_thickness

    substrate = gmsh.model.occ.addRectangle(
        0.0,
        0.0,
        0.0,
        length,
        substrate_thickness,
    )
    box = gmsh.model.occ.addRectangle(
        0.0,
        box_bottom,
        0.0,
        length,
        box_thickness,
    )
    film = gmsh.model.occ.addRectangle(
        0.0,
        film_bottom,
        0.0,
        length,
        film_thickness,
    )
    front_oxide = gmsh.model.occ.addRectangle(
        gate_start,
        film_top,
        0.0,
        gate_length,
        front_oxide_thickness,
    )

    gmsh.model.occ.fragment(
        [(2, substrate), (2, box), (2, film)],
        [(2, front_oxide)],
    )
    gmsh.model.occ.synchronize()

    tolerance = 1.0e-6 * max(length, oxide_top, 1.0)

    substrate_surfaces = []
    box_surfaces = []
    film_surfaces = []
    front_oxide_surfaces = []

    for _, tag in gmsh.model.getEntities(2):
        _, ymin, _, _, ymax, _ = gmsh.model.getBoundingBox(2, tag)
        center_y = 0.5 * (ymin + ymax)

        if center_y < substrate_thickness - tolerance:
            substrate_surfaces.append(tag)
        elif center_y < box_top - tolerance:
            box_surfaces.append(tag)
        elif center_y < film_top + tolerance:
            film_surfaces.append(tag)
        else:
            front_oxide_surfaces.append(tag)

    if not substrate_surfaces:
        raise RuntimeError("No substrate surface found.")
    if not box_surfaces:
        raise RuntimeError("No BOX surface found.")
    if not film_surfaces:
        raise RuntimeError("No silicon film surface found.")
    if not front_oxide_surfaces:
        raise RuntimeError("No front oxide surface found.")

    add_physical_group(2, film_surfaces, SILICON_FILM_PHYSICAL_NAME)
    add_physical_group(2, substrate_surfaces, SUBSTRATE_PHYSICAL_NAME)
    add_physical_group(2, front_oxide_surfaces, FRONT_OXIDE_PHYSICAL_NAME)
    add_physical_group(2, box_surfaces, BOX_PHYSICAL_NAME)

    def is_horizontal_at(curve_tag, y_target):
        xmin, ymin, xmax, ymax = curve_bbox(curve_tag)
        center_y = 0.5 * (ymin + ymax)
        return (xmax - xmin > tolerance) and (abs(center_y - y_target) <= tolerance)

    all_curves = [tag for _, tag in gmsh.model.getEntities(1)]

    film_top_curves = [tag for tag in all_curves if is_horizontal_at(tag, film_top)]
    oxide_top_curves = [tag for tag in all_curves if is_horizontal_at(tag, oxide_top)]
    substrate_bottom_curves = [tag for tag in all_curves if is_horizontal_at(tag, 0.0)]

    source_curves = []
    drain_curves = []

    for curve_tag in film_top_curves:
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
    if not oxide_top_curves:
        raise RuntimeError("No gate boundary found.")
    if not substrate_bottom_curves:
        raise RuntimeError("No back-gate boundary found.")

    add_physical_group(1, source_curves, "source")
    add_physical_group(1, drain_curves, "drain")
    add_physical_group(1, oxide_top_curves, "gate")
    add_physical_group(1, substrate_bottom_curves, "back_gate")

    return gate_start, gate_end


def set_mesh_fields(
    length,
    substrate_thickness,
    box_thickness,
    film_thickness,
    front_oxide_thickness,
    source_length,
    drain_length,
    h_min,
    h_max,
    interface_refinement_length,
    junction_refinement_length,
    gate_edge_refinement_length,
    source_drain_depth,
):
    gate_start = source_length
    gate_end = length - drain_length

    box_bottom = substrate_thickness
    box_top = substrate_thickness + box_thickness
    film_bottom = box_top
    film_top = film_bottom + film_thickness
    oxide_top = film_top + front_oxide_thickness

    fields = []

    def add_box(x_min, x_max, y_min, y_max, h_in=None):
        field = gmsh.model.mesh.field.add("Box")
        gmsh.model.mesh.field.setNumber(field, "VIn", h_min if h_in is None else h_in)
        gmsh.model.mesh.field.setNumber(field, "VOut", h_max)
        gmsh.model.mesh.field.setNumber(field, "XMin", max(0.0, x_min))
        gmsh.model.mesh.field.setNumber(field, "XMax", min(length, x_max))
        gmsh.model.mesh.field.setNumber(field, "YMin", max(0.0, y_min))
        gmsh.model.mesh.field.setNumber(field, "YMax", min(oxide_top, y_max))
        fields.append(field)

    add_box(
        0.0,
        length,
        film_bottom - interface_refinement_length,
        oxide_top,
    )

    add_box(
        0.0,
        length,
        box_top - interface_refinement_length,
        film_top + interface_refinement_length,
    )

    add_box(
        0.0,
        length,
        box_bottom - interface_refinement_length,
        box_top + interface_refinement_length,
        h_in=max(h_min, 0.5 * h_max),
    )

    for x0 in (gate_start, gate_end):
        add_box(
            x0 - gate_edge_refinement_length,
            x0 + gate_edge_refinement_length,
            film_bottom,
            oxide_top,
        )

    for x0 in (source_length, length - drain_length):
        add_box(
            x0 - junction_refinement_length,
            x0 + junction_refinement_length,
            film_top - min(source_drain_depth, film_thickness) - junction_refinement_length,
            film_top,
        )

    minimum_field = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(minimum_field, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(minimum_field)


def generate_mesh(
    mesh_file,
    length,
    substrate_thickness,
    box_thickness,
    film_thickness,
    front_oxide_thickness,
    source_length,
    drain_length,
    h_min,
    h_max,
    interface_refinement_length,
    junction_refinement_length,
    gate_edge_refinement_length,
    n_plus_level,
    channel_acceptor_level,
    substrate_acceptor_level,
    background_donor_level,
    source_drain_depth,
    lateral_smoothing_length,
    vertical_smoothing_length,
    min_doping,
    max_doping,
    temperature_source,
    temperature_drain,
):
    model_name = "FDSOI_NMOS_2D"

    gmsh.initialize()

    try:
        gmsh.model.add(model_name)

        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.MeshSizeMin", h_min)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h_max)

        create_geometry(
            length=length,
            substrate_thickness=substrate_thickness,
            box_thickness=box_thickness,
            film_thickness=film_thickness,
            front_oxide_thickness=front_oxide_thickness,
            source_length=source_length,
            drain_length=drain_length,
        )

        set_mesh_fields(
            length=length,
            substrate_thickness=substrate_thickness,
            box_thickness=box_thickness,
            film_thickness=film_thickness,
            front_oxide_thickness=front_oxide_thickness,
            source_length=source_length,
            drain_length=drain_length,
            h_min=h_min,
            h_max=h_max,
            interface_refinement_length=interface_refinement_length,
            junction_refinement_length=junction_refinement_length,
            gate_edge_refinement_length=gate_edge_refinement_length,
            source_drain_depth=source_drain_depth,
        )

        gmsh.model.mesh.generate(2)

        node_tags, node_coordinates, _ = gmsh.model.mesh.getNodes()
        x = node_coordinates[0::3]
        y = node_coordinates[1::3]

        acceptor, donor, net_doping = compute_doping_profile(
            x=x,
            y=y,
            length=length,
            substrate_thickness=substrate_thickness,
            box_thickness=box_thickness,
            film_thickness=film_thickness,
            source_length=source_length,
            drain_length=drain_length,
            n_plus_level=n_plus_level,
            channel_acceptor_level=channel_acceptor_level,
            substrate_acceptor_level=substrate_acceptor_level,
            background_donor_level=background_donor_level,
            source_drain_depth=source_drain_depth,
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
    substrate_thickness,
    box_thickness,
    film_thickness,
    source_length,
    drain_length,
    n_plus_level,
    channel_acceptor_level,
    substrate_acceptor_level,
    background_donor_level,
    source_drain_depth,
    lateral_smoothing_length,
    vertical_smoothing_length,
    min_doping,
    max_doping,
    show_plot,
):
    box_top = substrate_thickness + box_thickness
    film_bottom = box_top
    film_top = film_bottom + film_thickness

    x = np.linspace(0.0, length, 1200)
    y_front = np.full_like(x, film_top)
    y_mid = np.full_like(x, 0.5 * (film_bottom + film_top))
    y_back = np.full_like(x, film_bottom)

    rows = []
    labels = [
        ("front", y_front),
        ("mid_film", y_mid),
        ("back", y_back),
    ]

    profiles = {}

    for label, y_values in labels:
        acceptor, donor, net = compute_doping_profile(
            x=x,
            y=y_values,
            length=length,
            substrate_thickness=substrate_thickness,
            box_thickness=box_thickness,
            film_thickness=film_thickness,
            source_length=source_length,
            drain_length=drain_length,
            n_plus_level=n_plus_level,
            channel_acceptor_level=channel_acceptor_level,
            substrate_acceptor_level=substrate_acceptor_level,
            background_donor_level=background_donor_level,
            source_drain_depth=source_drain_depth,
            lateral_smoothing_length=lateral_smoothing_length,
            vertical_smoothing_length=vertical_smoothing_length,
            min_doping=min_doping,
            max_doping=max_doping,
        )
        profiles[label] = (acceptor, donor, net)
        rows.extend([acceptor, donor, net])

    profile_file = output_prefix.with_suffix(".channel_profile.csv")
    np.savetxt(
        profile_file,
        np.column_stack([x] + rows),
        delimiter=",",
        header=(
            "x,"
            "FrontAcceptorConcentration,FrontDonorConcentration,FrontDopingConcentration,"
            "MidFilmAcceptorConcentration,MidFilmDonorConcentration,MidFilmDopingConcentration,"
            "BackAcceptorConcentration,BackDonorConcentration,BackDopingConcentration"
        ),
        comments="",
    )

    fig, ax = plt.subplots()
    for label, linestyle in (("front", "-"), ("mid_film", "--"), ("back", ":")):
        acceptor, donor, net = profiles[label]
        ax.plot(x, np.clip(donor, min_doping, None), label=f"Donor {label}", linestyle=linestyle)
        ax.plot(x, np.clip(acceptor, min_doping, None), label=f"Acceptor {label}", linestyle=linestyle)
        ax.plot(x, np.clip(np.abs(net), min_doping, None), label=f"|Net| {label}", linestyle=linestyle)

    ax.set_xlabel("x (µm)")
    ax.set_ylabel("Concentration (cm⁻³)")
    ax.set_yscale("log")
    ax.set_ylim(min_doping, 2.0 * max(n_plus_level, channel_acceptor_level, substrate_acceptor_level, background_donor_level))
    ax.set_title(output_prefix.name)
    ax.grid(True, which="both")
    ax.legend(fontsize="small", ncol=2)

    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".doping_profile.png"), dpi=200)

    if show_plot:
        plt.show()

    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a 2D FD-SOI nMOS mesh with front oxide, thin silicon film, BOX, and back gate."
    )

    parser.add_argument("mesh_name", help="Output mesh filename, for example fdsoi_nmos.msh.")

    parser.add_argument("--length", type=float, default=0.20, help="Total device length in µm.")
    parser.add_argument("--substrate-thickness", type=float, default=0.05, help="Back-plane silicon thickness in µm.")
    parser.add_argument("--box-thickness", type=float, default=0.025, help="Buried oxide thickness in µm.")
    parser.add_argument("--film-thickness", type=float, default=0.010, help="Thin silicon film thickness in µm.")
    parser.add_argument("--front-oxide-thickness", type=float, default=0.002, help="Front gate oxide thickness in µm.")
    parser.add_argument("--source-length", type=float, default=0.050, help="Source access and implant length in µm.")
    parser.add_argument("--drain-length", type=float, default=0.050, help="Drain access and implant length in µm.")

    parser.add_argument("--hmin", type=float, default=0.0005, help="Minimum mesh size in µm.")
    parser.add_argument("--hmax", type=float, default=0.0050, help="Maximum mesh size in µm.")
    parser.add_argument(
        "--interface-refinement-length",
        type=float,
        default=0.005,
        help="Vertical refinement distance around oxide/silicon interfaces in µm.",
    )
    parser.add_argument(
        "--junction-refinement-length",
        type=float,
        default=0.010,
        help="Refinement half-width around source/channel and drain/channel junctions in µm.",
    )
    parser.add_argument(
        "--gate-edge-refinement-length",
        type=float,
        default=0.010,
        help="Refinement half-width around the two gate edges in µm.",
    )

    parser.add_argument("--n-plus-level", type=float, default=1.0e19, help="Peak N+ source/drain donor concentration in cm^-3.")
    parser.add_argument("--channel-acceptor-level", type=float, default=1.0e15, help="Thin-film channel acceptor concentration in cm^-3.")
    parser.add_argument("--substrate-acceptor-level", type=float, default=1.0e17, help="Back-plane substrate acceptor concentration in cm^-3.")
    parser.add_argument("--background-donor-level", type=float, default=1.0e11, help="Background donor concentration in silicon in cm^-3.")
    parser.add_argument("--source-drain-depth", type=float, default=0.010, help="Source/drain implant depth from the front film surface in µm.")
    parser.add_argument("--lateral-smoothing-length", type=float, default=0.004, help="Lateral source/drain smoothing length in µm.")
    parser.add_argument("--vertical-smoothing-length", type=float, default=0.002, help="Vertical source/drain smoothing length in µm.")
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
        "--substrate-thickness": args.substrate_thickness,
        "--box-thickness": args.box_thickness,
        "--film-thickness": args.film_thickness,
        "--front-oxide-thickness": args.front_oxide_thickness,
        "--source-length": args.source_length,
        "--drain-length": args.drain_length,
        "--hmin": args.hmin,
        "--hmax": args.hmax,
        "--interface-refinement-length": args.interface_refinement_length,
        "--junction-refinement-length": args.junction_refinement_length,
        "--gate-edge-refinement-length": args.gate_edge_refinement_length,
        "--n-plus-level": args.n_plus_level,
        "--channel-acceptor-level": args.channel_acceptor_level,
        "--substrate-acceptor-level": args.substrate_acceptor_level,
        "--background-donor-level": args.background_donor_level,
        "--source-drain-depth": args.source_drain_depth,
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

    if args.source_drain_depth > args.film_thickness:
        raise ValueError("--source-drain-depth cannot be larger than --film-thickness.")

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
            substrate_thickness=args.substrate_thickness,
            box_thickness=args.box_thickness,
            film_thickness=args.film_thickness,
            source_length=args.source_length,
            drain_length=args.drain_length,
            n_plus_level=args.n_plus_level,
            channel_acceptor_level=args.channel_acceptor_level,
            substrate_acceptor_level=args.substrate_acceptor_level,
            background_donor_level=args.background_donor_level,
            source_drain_depth=args.source_drain_depth,
            lateral_smoothing_length=args.lateral_smoothing_length,
            vertical_smoothing_length=args.vertical_smoothing_length,
            min_doping=args.min_doping,
            max_doping=args.max_doping,
            show_plot=args.show,
        )

    generate_mesh(
        mesh_file=mesh_file,
        length=args.length,
        substrate_thickness=args.substrate_thickness,
        box_thickness=args.box_thickness,
        film_thickness=args.film_thickness,
        front_oxide_thickness=args.front_oxide_thickness,
        source_length=args.source_length,
        drain_length=args.drain_length,
        h_min=args.hmin,
        h_max=args.hmax,
        interface_refinement_length=args.interface_refinement_length,
        junction_refinement_length=args.junction_refinement_length,
        gate_edge_refinement_length=args.gate_edge_refinement_length,
        n_plus_level=args.n_plus_level,
        channel_acceptor_level=args.channel_acceptor_level,
        substrate_acceptor_level=args.substrate_acceptor_level,
        background_donor_level=args.background_donor_level,
        source_drain_depth=args.source_drain_depth,
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

