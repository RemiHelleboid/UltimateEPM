import matplotlib as mpl
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from matplotlib.lines import Line2D
from scipy.interpolate import CubicSpline
from pathlib import Path
import os
import glob
import sys
import re
import matplotlib.cm as cm

def get_crystal_dir(list_q):
    dir = 111
    for q in list_q:
        if q[2] != q[1] or q[2] != q[0]:
            dir = 110
            break
    for q in list_q:
        print(q)
        if q[1] != q[0]:
            dir = 100
            break
    return dir


def parse_epsilon_files(dirname):
    files = glob.glob(dirname + "/*.csv")
    list_qxyz = []
    list_norm_q = []
    data = []
    list_energies = []
    dir = 0
    for idx, filename in enumerate(files):
        print(f"\r File n° {idx} / {len(files)}", end="", flush=True)
        stem_name = Path(filename).stem
        qx, qy, qz = stem_name.split(
            "_")[-3], stem_name.split("_")[-2], stem_name.split("_")[-1]
        norm_q = np.sqrt(float(qx)**2 + float(qy)**2 + float(qz)**2)
        list_norm_q.append(norm_q)
        # print(f"q = (qx, qy, qz) = ({qx}, {qy}, {qz})")
        list_qxyz.append((qx, qy, qz))
        energies, eps_r = np.loadtxt(
            filename, unpack=True, skiprows=1, delimiter=',')
        # print("Type and shape eps_r: ", type(eps_r), eps_r.shape)
        list_energies = energies
        data.append(list(eps_r))
    dir = get_crystal_dir(list_qxyz)
    data_sorted = [x for _, x in sorted(zip(list_norm_q, data))]
    array_data = np.array(data_sorted)
    list_norm_q.sort()
    return dir, list_norm_q, list_energies, array_data

def plot_q_epsilon_map(dirname):
    fig, ax = plt.subplots()
    dir, list_q, list_energies, data = parse_epsilon_files(dirname)
    STEP = len(list_q) // 8 
    MyMap = cm.get_cmap('autumn')
    for i, eps_r in enumerate(data[::STEP]):
        color = MyMap(float(list_q[STEP*i]))
        ax.plot(list_energies, eps_r,
                label=f"{float(list_q[STEP*i]):.2f}", alpha=0.85, c=color)
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel('$\epsilon$')
    ax.legend(title="$\left| \mathbf{q} \\right|$")
    ax.set_title(f"Silicon dielectric function in the ({dir}) direction.")
    fig.tight_layout()
    filename = f"plot_epsilon_at_diff_q_{dir}"
    fig.savefig(filename + ".pdf")
    fig.savefig(filename + ".png")
    # plt.show()


def plot_heatmap(dirname):
    fig, ax = plt.subplots()
    dir, list_q, list_energies, data = parse_epsilon_files(dirname)
    list_q = np.array(list_q, dtype=np.float64)
    print("Shape of data: ", data.shape)
    # datac = np.array(data)
    list_energies = np.array(list_energies, dtype=np.float64)
    # print(list_q.shape)
    X, Y = np.meshgrid(list_q, list_energies)
    map = 'jet'
    pc = ax.pcolormesh(Y, X, data.T, cmap=map, shading='auto', rasterized=True)
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("$\left| \mathbf{q} \\right|$")
    ax.set_title(f"Silicon dielectric function in the {dir} direction.")
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
    fig.tight_layout()
    filename = f"plot_epsilon_maps_q_{dir}"
    fig.savefig(filename + ".pdf")
    fig.savefig(filename + ".png", dpi=300)

    # plt.show()


def plot_3d_surface(dirname):
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    dir, list_q, list_energies, data = parse_epsilon_files(dirname)
    list_q = np.array(list_q, dtype=np.float64)
    # data = np.array(data)
    list_energies = np.array(list_energies, dtype=np.float64)
    # print(list_q.shape)
    X, Y = np.meshgrid(list_q, list_energies)
    surf = ax.plot_surface(X, Y, data.T, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    plt.show()


if __name__ == '__main__':
    dirname = sys.argv[1]
    # qx, data = parse_epsilon_files(dirname)
    plot_q_epsilon_map(dirname)
    plot_heatmap(dirname)
    # plot_3d_surface(dirname)
