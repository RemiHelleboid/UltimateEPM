import matplotlib as mpl
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from matplotlib.lines import Line2D
from scipy.interpolate import CubicSpline
from pathlib import Path


import sys

try:
    import scienceplots
    plt.style.use(['science', 'high-vis'])
except Exception:
    print('Could not load style sheet')
    None
    plt.style.use(['seaborn-paper'])

mpl.rcParams['figure.figsize'] = [3.5, 2.8]


def plot_dielectric_function_vs_energy(filename: str, ax):
    # energies, eps_r, eps_i = np.loadtxt(filename, unpack=True, skiprows=1, delimiter=',')
    data = np.loadtxt(filename, unpack=True, skiprows=1, delimiter=',')
    energies = data[0]
    eps_r = data[1]
    if len(data) == 3:
        eps_i = data[2]
    else:
        eps_i = np.zeros_like(eps_r)
    label = Path(filename).stem.split('_')[2::]
    label = ["{:.0f}".format(float(x)) for x in label]
    s = "$\mathbf{k} = (" + " \ ".join(label) + ")$"
    ax.plot(energies, eps_r, label=s)
    # ax.plot(energies, eps_i, "-", label="Imaginary", c='r')
    abs_eps = np.sqrt(eps_r**2 + eps_i**2)
    # ax.plot(energies, abs_eps, label=s)
    
    
if __name__ == '__main__':
    list_files = sys.argv[1:]
    fig, ax = plt.subplots()
    for filename in list_files:
        plot_dielectric_function_vs_energy(filename, ax)
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('$\epsilon$')
    ax.legend()
    # ax.set_title('Eps vs Energy')
    fig.tight_layout()
    filename =f"{Path(list_files[0]).with_suffix('.png')}"
    fig.savefig(filename, dpi=600)
    plt.show()