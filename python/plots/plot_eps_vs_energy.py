import matplotlib as mpl
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from matplotlib.lines import Line2D
from scipy.interpolate import CubicSpline
from pathlib import Path

filename = "mysequence.fasta"
new_filename = Path(filename).stem + ".aln"


try:
    plt.style.use(['science', 'high-vis'])
except Exception:
    print('Could not load style sheet')
    None
    plt.style.use(['seaborn-paper'])

mpl.rcParams['figure.figsize'] = [3.5, 2.8]

def plot_dielectric_function_vs_energy(filename: str):
    energies, eps = np.loadtxt(filename, unpack=True, skiprows=1, delimiter=',')
    fig, ax = plt.subplots()
    ax.plot(energies, eps, label="Real part")
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
    # plt.show()