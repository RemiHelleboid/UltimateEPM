import numpy as np
import scipy.stats as st
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm, title
import scipy.stats as st
from scipy.stats import skewnorm
import glob, os
from argparse import ArgumentParser

plt.style.use(['science', 'high-vis'])


def plot_band_structure(filename):
    cnv = {1: lambda s: np.float(s.strip() or 'Nan')}
    band_energies = np.loadtxt(filename, delimiter=" ", usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8))
    band_energies = band_energies.T   
    fig, ax = plt.subplots()
    ax.set_title(f"Band structure from {filename}")
    ax.set_xlabel("$\mathbf{k}$")
    ax.set_ylabel("Energy (eV)")
    
    for band in band_energies[::]:
        ax.plot(band)
    plt.show()
    

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-f", "--file", dest="path_file", help="The file to parse.", required=True)
    parser.add_argument("-o", "--outputdir", dest="output_path_dir", help="The directory to save the results.", default="./")
    args = parser.parse_args()
    FILE_PATH = args.path_file
    
    plot_band_structure(FILE_PATH)

    