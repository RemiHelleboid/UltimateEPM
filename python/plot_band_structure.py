import numpy as np
import scipy.stats as st
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm, title
import scipy.stats as st
from scipy.stats import skewnorm
import glob
import os
from argparse import ArgumentParser

plt.style.use(['science', 'high-vis'])

BZ_points = {
    "G":  np.array([0, 0, 0]),
    "X":  np.array([0, 1/2, 0]),
    "L":  np.array([1/4, 1/4, 1/4]),
    "W":  np.array([1/4, 1/2, 0]),
    "U":  np.array([1/8, 1/2, 1/8]),
    "K":  np.array([3/8, 3/8, 0]),
}

def get_path_from_filename(filename):
    spliteed_filename = filename.split("_")
    index_path = spliteed_filename.index("path") + 1 
    path = spliteed_filename[index_path]
    list_points_string = list(path)
    list_points_plot = [point if point != "G" else "$\Gamma$" for point in list_points_string ]
    print(f"Path: {list_points_string}")
    list_points = [BZ_points[point] for point in list_points_string]
    print(f"Path: {list_points}")
    dist_btw_points = [np.linalg.norm(list_points[k+1] - list_points[k]) for k in range(len(list_points)-1)]
    print(f"Distance between points: {dist_btw_points}")
    return list_points_plot, dist_btw_points

def plot_band_structure(filename):
    list_points_string, dist_btw_points = get_path_from_filename(filename)
    point_sym_positions = [0]
    for k in range(1, len(list_points_string)):
        point_sym_positions.append(point_sym_positions[k-1] + dist_btw_points[k-1])
    
    print(f"Point sym positions: {list_points_string}")
    
    cnv = {1: lambda s: np.float(s.strip() or 'Nan')}
    nb_bands = 10
    band_energies = np.loadtxt(
        filename, delimiter=" ", usecols=tuple(i for i in range(nb_bands)))
    band_energies = band_energies.T
    fig, ax = plt.subplots()
    # ax.set_title(f"Band structure from {filename}")
    ax.set_xlabel("$\mathbf{k}$")
    ax.set_ylabel("Energy (eV)")

    for band in band_energies[::]:
        ax.plot(band, ls="-", color='k', lw=0.5)

    ax.set_ylim(bottom=-14, top=10)
    ax.grid(True, which='both', axis='both', ls="-",
            lw=0.25, alpha=0.5, color='grey')
    point_sym_positions = np.array(point_sym_positions)
    point_sym_positions /= np.max(point_sym_positions)
    point_sym_positions *= len(band_energies[0])
    ax.set_xticks(point_sym_positions, labels=list_points_string)
    ax.set_xticks([], minor=True)
    
    print(point_sym_positions)
    fig.tight_layout()
    
    # fig.savefig(f"{filename[:-4:]}.png", dpi=600)
    fig.savefig(f"{filename[:-4:]}.pdf", dpi=600)
    plt.show()


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-f", "--file", dest="path_file",
                        help="The file to parse.", required=True)
    parser.add_argument("-o", "--outputdir", dest="output_path_dir",
                        help="The directory to save the results.", default="./")
    args = parser.parse_args()
    FILE_PATH = args.path_file

    plot_band_structure(FILE_PATH)