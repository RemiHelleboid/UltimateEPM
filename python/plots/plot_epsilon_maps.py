import matplotlib as mpl
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from matplotlib.lines import Line2D
from scipy.interpolate import CubicSpline
from pathlib import Path
import os, glob, sys
import re

def parse_epsilon_files(dirname):
    files = glob.glob(dirname + "/epsilon*.csv")
    list_qx = []
    data = []
    for filename in files:
        qx = re.findall(r"Qx[-+]?(?:\d*\.\d+|\d+)", filename)[0].replace("Qx", "")
        print(qx)
        list_qx.append(qx)
        energies, eps_r = np.loadtxt(filename, unpack=True, skiprows=1, delimiter=',')
        data.append([energies, eps_r])
    return list_qx, data

def plot_q_epsilon_map(dirname):
    fig, ax = plt.subplots()
    list_q, data = parse_epsilon_files(dirname)
    for i, (energies, eps_r) in enumerate(data):
        ax.plot(energies, eps_r, label=f"Qx={list_q[i]}")
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('$\epsilon$')
    ax.legend()
    # ax.set_title('Eps vs Energy')
    fig.tight_layout()
    filename =f"plot_epsilon_map.pdf"
    fig.savefig(filename + ".pdf")
    plt.show()
    
    
if __name__ == '__main__':
    dirname = sys.argv[1]
    qx, data = parse_epsilon_files(dirname)
    plot_q_epsilon_map(dirname)