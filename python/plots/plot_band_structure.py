import matplotlib as mpl
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from matplotlib.lines import Line2D
from scipy.interpolate import CubicSpline


try:
    plt.style.use(['science', 'high-vis'])
except Exception:
    None
plt.style.use(['seaborn-paper'])

mpl.rcParams['figure.figsize'] = [3.5, 2.8]


BZ_points = {
    "G":  np.array([0, 0, 0]),
    "X":  np.array([0, 1/2, 0]),
    "L":  np.array([1/4, 1/4, 1/4]),
    "W":  np.array([1/4, 1/2, 0]),
    "U":  np.array([1/8, 1/2, 1/8]),
    "K":  np.array([3/8, 3/8, 0]),
}

def get_epm_parameters_from_file(filename: str) -> dict:
    """Extract the EPM parameters from the band structure file.

    Args:
        filename (str): Filename of the band structure file.

    Returns:
        dic: The EPM p0arameters.
    """
    comment_char = '#'
    dict_params = {}
    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            if comment_char in line:
                line = line.split(comment_char)[1]
                split = line.split(' ')
                if len(split) == 2:
                    dict_params[split[0]] = split[1]
                elif len(split) > 3:
                    dict_params[split[0]] = split[1:]
            else:
                break
    return dict_params


def get_gap(bands: np.array) -> tuple:
    """Get the gap of the band structure.

    Args:
        bands (np.array): Band structure.

    Returns:
        tuple: The gap of the band structure.
    """
    EPSILON = 1.0e-6
    max_valence = -np.inf
    min_conduction = np.inf
    for band in bands:
        if np.all(band < EPSILON):
            max_valence = max(max_valence, max(band))
        else:
            min_conduction = min(min_conduction, min(band))
    gap = min_conduction - max_valence
    print(f"Min conduction: {min_conduction}")
    print(f"Max valence: {max_valence}")
    return gap


def get_path_from_filename(filename: str) -> tuple:
    """Extract the path from the filename produced by the band structure calculation.
    Also compute the distance between each point of the path.

    Args:
        filename(str): Filename of the band structure file.

    Returns:
        tuple: The result.
    """
    filename = Path(filename).stem
    print(f"Filename: {filename}")
    spliteed_filename = filename.split("_")
    index_material = spliteed_filename.index("EPM") + 1
    index_path = spliteed_filename.index("path") + 1
    material = spliteed_filename[index_material]
    path = spliteed_filename[index_path]
    list_points_string = list(path)
    list_points_plot = [point if point !=
                        "G" else "$\Gamma$" for point in list_points_string]
    print(f"Path: {list_points_string}")
    list_points = [BZ_points[point] for point in list_points_string]
    dist_btw_points = [np.linalg.norm(
        list_points[k+1] - list_points[k]) for k in range(len(list_points)-1)]

    # Treatment of UK path, which distance is set to 0.0
    if "UK" in path:
        index_UK = path.index("UK")
        dist_btw_points[index_UK] = 0.0
    return list_points_plot, dist_btw_points, material


def plot_band_structure(filename: str, ax, index_plot, nb_bands=10):
    """_summary_

    Args:
        filename (str): Filename of the band structure file.
        out_dir (str, optional): Path of the directory used to stor the plot figure. Defaults to ".".
        nb_bands (int, optional): Number of band to keep in the plot. Defaults to 10.
        plot (bool, optional): Wheter or not to show the plot in a pyplot windows. Defaults to True.

    Returns:
        _type_: _description_
    """
    dict_params = get_epm_parameters_from_file(filename)
    print(len(dict_params))
    print("Plotting band structure")
    list_points_string, dist_btw_points, material = get_path_from_filename(
        filename)
    point_sym_positions = [0]
    linestyles = ["-", "--", "-.", ":"]
    plot_ls = linestyles[index_plot % len(linestyles)]
    for k in range(1, len(list_points_string)):
        point_sym_positions.append(
            point_sym_positions[k-1] + dist_btw_points[k-1])

    cnv = {1: lambda s: np.float(s.strip() or 'Nan')}
    band_energies = np.loadtxt(
        filename, delimiter=",", usecols=tuple(i for i in range(nb_bands)), skiprows=5)
    band_energies = band_energies.T
    print(dist_btw_points)
    

    for band in band_energies[::]:
        ax.plot(band, ls=plot_ls, color='k', lw=0.5)

    band_gap = get_gap(band_energies)
    print(f"Band gap: {band_gap}")

    return band_gap


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-f", "--file", dest="band_files",
                        help="The file to parse.", required=True, nargs="+")
    parser.add_argument("-b", "--nbbands", dest="nb_bands", type=int,
                        help="The nb of bands to plot.", required=False, default=10)
    parser.add_argument("-o", "--outputdir", dest="output_path_dir",
                        help="The directory to save the results.", default="./")
    parser.add_argument("-p", "--plot", dest="show_plot",
                        help="If true, the plot is displayed.", default=False, type=bool)
    args = parser.parse_args()
    OUT_DIR = args.output_path_dir

    list_files = args.band_files
    FILE_PATH_0 = list_files[0]
    list_points_string, dist_btw_points, material = get_path_from_filename(
        FILE_PATH_0)
    nb_points = num_lines = sum(1 for line in open(FILE_PATH_0)) - 1

    fig, ax = plt.subplots()
    band_gap = 0.0
    for index, file in enumerate(list_files):
        get_epm_parameters_from_file(file)
        band_gap = plot_band_structure(file, ax, index, nb_bands=args.nb_bands)

    # ax.set_ylim(bottom=-20, top=12)
    ax.grid(True, which='both', axis='both', ls="-",
            lw=0.25, alpha=0.5, color='grey')
    ax.set_title(f"Band structure of {material}")
    ax.set_xlabel("$\mathbf{k}$")
    ax.set_ylabel("Energy (eV)")
    point_sym_positions = [0]
    for k in range(1, len(list_points_string)):
        point_sym_positions.append(
            point_sym_positions[k-1] + dist_btw_points[k-1])
    point_sym_positions = np.array(point_sym_positions)
    point_sym_positions /= np.max(point_sym_positions)
    point_sym_positions *= nb_points

    ax.set_xticks(point_sym_positions, labels=list_points_string)
    ax.set_xticks([], minor=True)

    lines = [Line2D([0], [0], color='k', lw=0.5, ls=lstyle)
             for lstyle in ["-", "--"]]
    labels = ["non-local", "local"]

    ax.legend(lines, labels, loc='lower right', fancybox=True,  frameon=True,
              edgecolor='k', facecolor='w', fontsize=6, framealpha=0.75)
    fig.tight_layout()
    fig.savefig(f"{OUT_DIR}/band_structure_{material}.png", dpi=300)

    if args.show_plot:
        plt.show()

    ax.set_ylim(-2.0, band_gap + 1.0)
    fig.savefig(f"{OUT_DIR}/band_structure_{material}_zoom.png", dpi=300)

    # plot_band_structure(FILE_PATH, OUT_DIR, args.nb_bands)
