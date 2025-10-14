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

plt.style.use(['science', 'muted', 'grid', 'no-latex'])


def scatter_plot_rates(filename):
    fig, ax = plt.subplots(figsize=(8, 6))
    data = np.loadtxt(filename, delimiter=",")
    energy = data[:,1]
    nb_modes = data.shape[1] - 2
    colors = cm.viridis(np.linspace(0, 1, nb_modes))
    # ax.scatter(energy, np.sum(data[:,2:], axis=1), label="Total Rate", s=1, color="black")
    # energy -= np.min(energy)
    for i in range(2, data.shape[1]):
        ax.scatter(energy, data[:,i], label=f"Mode {i}", s=1, color=colors[i-2])  
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("Rate (s$^-1$)")
    
    fig.tight_layout()
    
    ax.legend()
    
def plot_rates(filename, comparison_csv=None):
    data = np.loadtxt(filename, delimiter=",", skiprows=1)
    energy = data[:,0]

    # Sum of rates
    fig2, ax2 = plt.subplots(figsize=(6, 4))
    ax2.plot(energy, np.sum(data[:,2:], axis=1), label="Total Rate", color="black")
    ax2.set_xlabel("Energy (eV)")
    ax2.set_ylabel("Rate (s$^-1$)")
    ax2.set_title("Total Electron-Phonon Rate")
    # ax2.plot(energy, data[:,1]/np.max(data[:,2:]), label="Wavefunction Overlap", color="blue", linestyle="--")

    if comparison_csv is not None:
        # Load comparison data
        energy, rates = np.loadtxt(comparison_csv, delimiter=",", unpack=True)
        ax2.plot(energy, rates, label="Comparison", color="red")
    ax2.legend()


    fig2.tight_layout()    
    
    
def plo_dos(filename):
    data = np.loadtxt(filename)
    energy = data[:,0]
    energy -= np.min(energy)
    dos = data[:,1]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(energy, dos, label="DOS", color="black", linestyle="--")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("DOS")
    fig.tight_layout()
    ax.legend()
    
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-f", "--filename", type=str, required=True, help="Filename to plot")
    args = parser.parse_args()
    try:
        scatter_plot_rates("rates_all.csv")
    except:
        pass
    try:
        plot_rates("rates_vs_energy.csv", "../examples/RatesSiFischetti1988.csv")
        plot_rates("rates_vs_energy.csv", "../examples/RatesSiKunikiyo1994.csv")
    except:
        pass
    plt.show()