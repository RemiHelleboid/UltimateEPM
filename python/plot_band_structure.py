import numpy as np
import scipy.stats as st
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm, title
import scipy.stats as st
from scipy.stats import skewnorm
import glob
import os
from argparse import ArgumentParser

try:
    plt.style.use(['science', 'high-vis'])
except Exception:
    plt.style.use(['seaborn'])

BZ_points = {
    "G":  np.array([0, 0, 0]),
    "X":  np.array([0, 1/2, 0]),
    "L":  np.array([1/4, 1/4, 1/4]),
    "W":  np.array([1/4, 1/2, 0]),
    "U":  np.array([1/8, 1/2, 1/8]),
    "K":  np.array([3/8, 3/8, 0]),
}

def get_path_from_filename(filename):
    filename = Path(filename).stem
    print(f"Filename: {filename}")
    spliteed_filename = filename.split("_")
    index_material = spliteed_filename.index("EEP") + 1 
    index_path = spliteed_filename.index("path") + 1 
    material = spliteed_filename[index_material]
    path = spliteed_filename[index_path]
    list_points_string = list(path)
    list_points_plot = [point if point != "G" else "$\Gamma$" for point in list_points_string ]
    print(f"Path: {list_points_string}")
    list_points = [BZ_points[point] for point in list_points_string]
    dist_btw_points = [np.linalg.norm(list_points[k+1] - list_points[k]) for k in range(len(list_points)-1)]
    return list_points_plot, dist_btw_points, material

def plot_band_structure(filename, OUT_DIR, nb_bands=10):
    list_points_string, dist_btw_points, material = get_path_from_filename(filename)
    point_sym_positions = [0]
    for k in range(1, len(list_points_string)):
        point_sym_positions.append(point_sym_positions[k-1] + dist_btw_points[k-1])
        
    cnv = {1: lambda s: np.float(s.strip() or 'Nan')}
    band_energies = np.loadtxt(
        filename, delimiter=" ", usecols=tuple(i for i in range(nb_bands)))
    band_energies = band_energies.T
    fig, ax = plt.subplots()
    ax.set_title(f"Band structure of {material}")
    ax.set_xlabel("$\mathbf{k}$")
    ax.set_ylabel("Energy (eV)")

    for band in band_energies[::]:
        ax.plot(band, ls="-", color='k', lw=0.5)
        print(f'[{band.min()},{band.max()}],')

    ax.set_ylim(bottom=-14, top=10)
    ax.grid(True, which='both', axis='both', ls="-",
            lw=0.25, alpha=0.5, color='grey')
    point_sym_positions = np.array(point_sym_positions)
    point_sym_positions /= np.max(point_sym_positions)
    point_sym_positions *= len(band_energies[0])
    ax.set_xticks(point_sym_positions, labels=list_points_string)
    ax.set_xticks([], minor=True)
    fig.tight_layout()
    filename = Path(filename).stem
    fig.savefig(f"{OUT_DIR}/{filename[:-4:]}.png", dpi=600)
    # fig.savefig(f"{filename[:-4:]}.pdf", dpi=600)
    # plt.show()


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-f", "--file", dest="path_file",
                        help="The file to parse.", required=True)
    parser.add_argument("-b", "--nbbands", dest="nb_bands", type=int,
                        help="The nb of bands to plot.", required=False)
    parser.add_argument("-o", "--outputdir", dest="output_path_dir",
                        help="The directory to save the results.", default="./")
    args = parser.parse_args()
    FILE_PATH = args.path_file
    OUT_DIR = args.output_path_dir

    plot_band_structure(FILE_PATH, OUT_DIR, args.nb_bands)
