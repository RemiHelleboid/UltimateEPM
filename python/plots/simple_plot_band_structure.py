import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from matplotlib.lines import Line2D
import matplotlib as mpl
from scipy.interpolate import CubicSpline
import os, sys


try:
    import scienceplots
    plt.style.use(['science', 'high-vis', "grid"])
except Exception:
    print("Could not use science style")
    plt.style.use(['seaborn-paper'])

mpl.rcParams['figure.figsize'] = [3.5, 2.8]



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
    
    band_energies = np.loadtxt(
        filename, delimiter=",", usecols=tuple(i for i in range(nb_bands)), skiprows=5)
    band_energies = band_energies.T
    for i in range(nb_bands):
        ax.plot(band_energies[i], color="black", lw=0.5, ls='-')
        
    


if __name__ == "__main__":
    file = sys.argv[1]
    fig, ax = plt.subplots()
    plot_band_structure(file, ax, 0)
    plt.show()
   
