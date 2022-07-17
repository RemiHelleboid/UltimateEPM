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

plt.style.use(['science', "high-vis"])


def plot_dos_per_band(filename, ax_plot=None, band_type="all"):
    if ax_plot is not None:
        axs = ax_plot
    else:
        fig, axs = plt.subplots()
    BANDS = np.loadtxt(filename, delimiter=',', skiprows=1)
    print(f"Shape bands data: {BANDS.shape}")
    number_bands = int((BANDS.shape[1] - 1) / 2)
    BANDS = BANDS.T
    band_counter = 0
    for band_index in range(0, number_bands):
        energies, dos = zip(*sorted(zip(BANDS[2*band_index], BANDS[2*band_index + 1])))
        if np.mean(energies) <= 0.0 and (band_type != "valence" and band_type != "all"):
            continue
        if np.mean(energies) >= 0.5 and (band_type!="conduction" and band_type!="all"):
            continue
        band_counter += 1
        axs.plot(energies, dos, lw=0.5 , label=f"{band_counter}")
    axs.legend(fontsize='xx-small', title_fontsize='x-small', title="Band index", fancybox=True)
    axs.set_xlabel("Energy (eV)")
    axs.set_ylabel("Density of state (a.u.)")
    if band_type == "conduction":
        axs.set_xlim(-1, )


def plot_dos_sum_bands(filename, ax_plot=None, band_type="all"):
    if ax_plot is not None:
        axs = ax_plot
    else:
        fig, axs = plt.subplots()
    BANDS = np.loadtxt(filename, delimiter=',', skiprows=1)
    print(f"Shape bands data: {BANDS.shape}")
    number_bands = int((BANDS.shape[1] - 1) / 2)
    BANDS = BANDS.T
    
    nb_points = 1000
    min_linspace = 0 if band_type in ["all", "conduction"] else -10.0
    max_linspace = 0 if band_type in ["all", "valence"] else 10.0
    energies_plot = np.linspace(min_linspace, max_linspace, nb_points)
    dos_total = np.zeros_like(energies_plot)

    count_band = 0
    for band_index in range(0, number_bands):
        energies, dos = zip(*sorted(zip(BANDS[2*band_index], BANDS[2*band_index + 1])))
        if np.mean(energies) <= 0.0 and (band_type != "valence" and band_type != "all"):
            continue
        if np.mean(energies) >= 0.5 and (band_type!="conduction" and band_type!="all"):
            continue
        dos_interp = np.interp(energies_plot, energies, dos)
        dos_total += dos_interp
        count_band += 1
        axs.plot(energies_plot, dos_total, label=f"{count_band}", lw=0.5)
    axs.legend(fontsize='x-small', title_fontsize='x-small', title="Number\n of bands", fancybox=True)
    axs.set_xlabel("Energy (eV)")
    axs.set_ylabel("Density of state (a.u.)")







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
    
    band_type = "valence"

    fig, axs = plt.subplots()
    plot_dos_per_band(FILE_PATH, axs, band_type)
    fig.tight_layout()
    fig.savefig("DOS_PER_BAND.pdf", dpi=600)
    fig.savefig("DOS_PER_BAND.png", dpi=600)
    plt.show()
    
    fig, axs = plt.subplots()
    plot_dos_sum_bands(FILE_PATH, axs, band_type)
    fig.tight_layout()
    fig.savefig("DOS_TOTAL.pdf", dpi=600)
    fig.savefig("DOS_TOTAL.png", dpi=600)
    plt.show()
