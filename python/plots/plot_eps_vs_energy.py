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


def plot_dielectric_function_vs_energy(dirname):
    list_files = glob.glob(f"{dirname}/*.csv")
    list_q = []
    for filename in list_files:
        filestem = Path(filename).stem
        qz = float(filestem.split('_')[-1])
        qy = float(filestem.split('_')[-2])
        qx = float(filestem.split('_')[-3])
        print(qx, qy, qz)
        list_q.append(np.sqrt(qx**2 + qy**2 + qz**2))
    print(list_q)
        
    list_files = [x for _, x in sorted(zip(list_q, list_files))]
    cmap = plt.get_cmap('jet')
    
    norm = mpl.colors.Normalize(vmin=min(list_q), vmax=max(list_q))
    
    fig, ax = plt.subplots()
    for idx, filename in enumerate(list_files):
        data = np.loadtxt(filename, unpack=True, skiprows=1, delimiter=',')
        energies = data[0]
        eps_r = data[1]
        
        ax.plot(energies, eps_r, color=cmap(norm(list_q[idx])))
        
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label='$|\mathbf{k}|$')
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('$\epsilon$')
    # ax.legend(title="$|\mathbf{k}|$", fontsize=8)
    # ax.set_title('Eps vs Energy')
    fig.tight_layout()
    filename =f"{Path(list_files[0]).with_suffix('.png')}"
    fig.savefig(filename, dpi=600)
    plt.show()

    
if __name__ == '__main__':
    dirname = Path(sys.argv[1])
    plot_dielectric_function_vs_energy(dirname)
    
    