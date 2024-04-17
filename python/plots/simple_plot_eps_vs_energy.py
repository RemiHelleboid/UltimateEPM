import matplotlib as mpl
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from matplotlib.lines import Line2D
from scipy.interpolate import CubicSpline
from pathlib import Path

import sys, glob

try:
    import scienceplots
    plt.style.use(['science', 'muted'])
except Exception:
    print('Could not load style sheet')
    None
    plt.style.use(['seaborn-paper'])

mpl.rcParams['figure.figsize'] = [3.5, 2.8]


def plot_dielectric_function_vs_energy(filename):
    energy, eps_r, eps_i = np.loadtxt(filename, delimiter=",", skiprows=1, unpack=True)
    fig, axs = plt.subplots()
    axs.plot(energy, eps_r, label="$\eps_r", c="r")
    axs.plot(energy, eps_i, label="$\eps_u", c='b')
    
    axs.set_ylabel("Relative dielectric function")
    axs.set_xlabel("Energy (eV)")
    fig.tight_layout(   )
    
    plt.show()
    
    
if __name__ == "__main__":
    file = sys.argv[1]
    plot_dielectric_function_vs_energy(file)

    