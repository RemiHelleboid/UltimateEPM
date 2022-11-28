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
import matplotlib.cm as cm

def parse_epsilon_files(dirname):
    files = glob.glob(dirname + "/epsilon*.csv")
    list_qx = []
    data = []
    list_energies = []
    for filename in files:
        qx = re.findall(r"Qx[-+]?(?:\d*\.\d+|\d+)", filename)[0].replace("Qx", "")
        print(qx)
        list_qx.append(qx)
        energies, eps_r = np.loadtxt(filename, unpack=True, skiprows=1, delimiter=',', dtype=np.float64)
        data.append(eps_r)
        list_energies = energies
    data_sorted = [x for _, x in sorted(zip(list_qx, data))]
    list_qx.sort()

    # data = np.array(data_sorted)
    # list_qx = np.array(list_qx)
    # print(data.shape)
    return list_qx, list_energies, data_sorted

def plot_q_epsilon_map(dirname):
    STEP = 20
    fig, ax = plt.subplots()
    list_q, list_energies, data = parse_epsilon_files(dirname)
    MyMap = cm.get_cmap('autumn')
    for i, eps_r in enumerate(data[::STEP]):
        color = MyMap(float(list_q[STEP*i]))
        ax.plot(list_energies, eps_r, label=f"{float(list_q[STEP*i]):.2f}", alpha=0.85, c=color)
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('$\epsilon$')
    ax.legend(title="$q_x$")
    # ax.set_title('Eps vs Energy')
    fig.tight_layout()
    filename =f"plot_epsilon_map.pdf"
    fig.savefig(filename + ".pdf")
    plt.show()
    
def plot_heatmap(dirname):
    fig, ax = plt.subplots()
    list_q, list_energies, data = parse_epsilon_files(dirname)
    list_q = np.array(list_q, dtype=np.float64)
    data = np.array(data)
    list_energies = np.array(list_energies, dtype=np.float64)
    # print(list_q.shape)
    X, Y = np.meshgrid(list_q,list_energies)
    map = 'jet'
    pc = ax.pcolormesh(Y, X, data.T, cmap=map, shading='auto', rasterized=True)
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("q$_x$")
    # ax.set_title("Silicon dielectric function in the (100) direction")
    fig.tight_layout()
    # cs = ax.contourf(Y, X, data.T, levels=25, cmap='seismic')
    fig.colorbar(pc, ax=ax, label="$\epsilon_r$")
    # ax.contour(cs, colors='k')
    # ax.imshow(data)
    # list_energies = data[0][0]
    # heat_map_e_q = np.zeros((list_q.shape[0], list_energies.shape[0]))
    # for idx_q in range(list_q.shape[0]):
    #     for idx_e in range(list_energies.shape[0]):
            
    #         heat_map_e_q[idx_q, idx_e] = data[idx_q, idx_e]
    # ax.imshow(heat_map_e_q)
    fig.savefig("silicon_epsilon_100.png", dpi=300)
    fig.savefig("silicon_epsilon_100.pdf")
    plt.show()
    
def plot_3d_surface(dirname):
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    list_q, list_energies, data = parse_epsilon_files(dirname)
    list_q = np.array(list_q, dtype=np.float64)
    data = np.array(data)
    list_energies = np.array(list_energies, dtype=np.float64)
    # print(list_q.shape)
    X, Y = np.meshgrid(list_q,list_energies)
    surf = ax.plot_surface(X, Y, data.T, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    plt.show()

   
   
if __name__ == '__main__':
    dirname = sys.argv[1]
    # qx, data = parse_epsilon_files(dirname)
    # plot_q_epsilon_map(dirname)
    plot_heatmap(dirname)
    # plot_3d_surface(dirname)