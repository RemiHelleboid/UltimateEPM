import matplotlib as mpl
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from matplotlib.lines import Line2D
from scipy.interpolate import CubicSpline
from pathlib import Path


try:
    plt.style.use(['science', 'high-vis'])
except Exception:
    print('Could not load style sheet')
    None
    plt.style.use(['seaborn-paper'])

mpl.rcParams['figure.figsize'] = [3.5, 2.8]


def plot_dielectric_function_vs_energy(filename: str):
    # energies, eps_r, eps_i = np.loadtxt(filename, unpack=True, skiprows=1, delimiter=',')
    data = np.loadtxt(filename, unpack=True, skiprows=1, delimiter=',')
    energies = data[0]
    eps_r = data[1]
    if len(data) == 3:
        eps_i = data[2]
    else:
        eps_i = np.zeros_like(eps_r)
    fig, ax = plt.subplots()
    ax.plot(energies, eps_r, label="Real", c='b')
    ax.plot(energies, eps_i, "-", label="Imaginary", c='r')
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('$\epsilon$')
    ax.legend()
    # ax.set_title('Eps vs Energy')
    fig.tight_layout()
    filename =f"{Path(filename).with_suffix('.pdf')}"
    fig.savefig(filename)
    
    
if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-f', '--filename', type=str, help='Filename of the dielectric function.')
    args = parser.parse_args()
    plot_dielectric_function_vs_energy(args.filename)
    plt.show()