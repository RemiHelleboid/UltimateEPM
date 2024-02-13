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


import matplotlib.style
import matplotlib as mpl

import scienceplots

plt.style.use(['science', 'muted', 'scatter', 'grid'])


def scatter_plot_rates(filename):
    fig, ax = plt.subplots()
    data = np.loadtxt(filename)
    energy = data[:,0]
    energy -= np.min(energy)
    for i in range(1, data.shape[1]):
        ax.scatter(energy, data[:,i], label=f"Mode {i}", s=1)   
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("Rate (s$^-1$)")
    
    fig.tight_layout()
    
    ax.legend()
    
def plot_rates(filename):
    fig, ax = plt.subplots()
    data = np.loadtxt(filename)
    energy = data[:,0]
    energy -= np.min(energy)
    dos = data[:,1]
    ax.plot(energy, dos, label="DOS")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("DOS")
    fig.tight_layout()
    ax.legend()
    
    
    
    
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-f", "--filename", type=str, required=True, help="Filename to plot")
    args = parser.parse_args()
    # scatter_plot_rates(args.filename)
    plot_rates(args.filename)
    plt.show()