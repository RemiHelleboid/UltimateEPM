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

plt.style.use(['science', 'muted', 'scatter', 'grid', 'no-latex'])


def scatter_plot_rates(filename):
    fig, ax = plt.subplots(figsize=(8, 6))
    data = np.loadtxt(filename)
    energy = data[:,0]
    # energy -= np.min(energy)
    for i in range(2, data.shape[1]):
        ax.scatter(energy, data[:,i], label=f"Mode {i}", s=1)   
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("Rate (s$^-1$)")
    
    fig.tight_layout()
    
    ax.legend()
    
def plot_rates(filename, comparison_csv=None):
    data = np.loadtxt(filename)
    energy = data[:,0]
    # energy -= np.min(energy)
    # Smallest strictly positive energy : 
    dos = data[:,1]
    min_energy = np.min(energy[dos > 0])
    energy -= min_energy
    nb_modes = data.shape[1] - 2
    fig, ax = plt.subplots(nb_modes, figsize=(6, 4*nb_modes))
    for i in range(1, data.shape[1]-3):
        ax[i-2].plot(energy, data[:,i], label=f"Mode {i}", color="black")
        ax[i-2].set_xlabel("Energy (eV)")
        ax[i-2].set_ylabel("Rate (s$^-1$)")
        ax[i-2].set_title(f"Mode {i}")

    # Sum of rates
    fig2, ax2 = plt.subplots(figsize=(6, 4))
    ax2.plot(energy, np.sum(data[:,2:], axis=1), label="Total Rate", color="black")
    ax2.set_xlabel("Energy (eV)")
    ax2.set_ylabel("Rate (s$^-1$)")
    ax2.set_title("Total Electron-Phonon Rate")

    if comparison_csv is not None:
        # Load comparison data
        energy, rates = np.loadtxt(comparison_csv, delimiter=",", unpack=True)
        ax2.plot(energy, rates, label="Comparison", color="red")
    ax2.legend()

    # max_rate = np.max(data[:,2:])
    # dos /= np.max(dos)
    # ax.plot(energy, dos, label="DOS", color="black", linestyle="--")
    
    # ax.legend()
    fig.tight_layout()    
    
    
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
    scatter_plot_rates("rates_all.csv", "../examples/RatesSiKunikiyo1994.csv")
    # plo_dos()
    plot_rates("rates_vs_energy.csv")
    plt.show()