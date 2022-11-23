import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def IsInIrreducibleWedge(k):
    return (k[2] >= 0.0 and k[2] <= k[1] and k[1] <= k[0] and k[0] <= 1.0) and \
        (np.sum(k) <= 3.0/2.0)


def Monkhorst_Pack(Nk1: int, Nk2: int, Nk3: int, irreducible_wedge=True) -> np.ndarray:
    b_1 = np.array([-1, 1, 1])
    b_2 = np.array([1, -1, 1])
    b_3 = np.array([1, 1, -1])
    list_k = []
    for i in range(0, Nk1):
        u_1 = i/Nk1
        for j in range(0, Nk2):
            u_2 = j/Nk2
            for k in range(0, Nk3):
                u_3 = k/Nk3
                k = u_1*b_1 + u_2*b_2 + u_3*b_3
                list_k.append(k)
    if irreducible_wedge:
        list_k = [x for x in list_k if IsInIrreducibleWedge(x)]
    # print(f"Number of kpoints in the irreducible wedge MP: {len(list_k)}")
    return np.array(list_k)

# def Monkhorst_Pack(Nk1: int, Nk2: int, Nk3: int, irreducible_wedge=True) -> np.ndarray:
#     b_1 = np.array([-1, 1, 1])
#     b_2 = np.array([1, -1, 1])
#     b_3 = np.array([1, 1, -1])
#     list_k = []
#     for i in range(0, Nk1):
#         u_1 = (2 * i - Nk1 + 1) / (2 * Nk1)
#         for j in range(0, Nk2):
#             u_2 = (2 * j - Nk2 + 1) / (2 * Nk2)
#             for k in range(0, Nk3):
#                 u_3 = (2 * k - Nk3 + 1) / (2 * Nk3)
#                 k = u_1 * b_1 + u_2 * b_2 + u_3 * b_3
#                 list_k.append(k)
#     if irreducible_wedge:
#         list_k = [x for x in list_k if IsInIrreducibleWedge(x)]
#     print(f"Number of kpoints in the irreducible wedge MP: {len(list_k)}")
#     return np.array(list_k)


def plot_segment(dot_1, dot_2, ax):
    ax.plot3D([dot_1[0], dot_2[0]], [dot_1[1], dot_2[1]],
              [dot_1[2], dot_2[2]], c='b')


def plot_bz_outline(ax):
    y_max_xp = 2.0 * np.array([1/4,  1/2,  0])
    y_max_xm = 2.0 * np.array([-1/4, 1/2,  0])
    y_max_zp = 2.0 * np.array([0,    1/2,  1/4])
    y_max_zm = 2.0 * np.array([0,    1/2,  -1/4])
    y_min_xp = 2.0 * np.array([1/4,  -1/2,  0])
    y_min_xm = 2.0 * np.array([-1/4, -1/2,  0])
    y_min_zp = 2.0 * np.array([0,    -1/2,  1/4])
    y_min_zm = 2.0 * np.array([0,    -1/2,  -1/4])

    x_max_yp = 2.0 * np.array([1/2, 1/4,   0])
    x_max_ym = 2.0 * np.array([1/2, -1/4,  0])
    x_max_zp = 2.0 * np.array([1/2, 0,     1/4])
    x_max_zm = 2.0 * np.array([1/2, 0,     -1/4])

    x_min_yp = 2.0 * np.array([-1/2, 1/4,   0])
    x_min_ym = 2.0 * np.array([-1/2, -1/4,  0])
    x_min_zp = 2.0 * np.array([-1/2, 0,     1/4])
    x_min_zm = 2.0 * np.array([-1/2, 0,     -1/4])

    z_max_yp = 2.0 * np.array([0,    1/4,   1/2])
    z_max_ym = 2.0 * np.array([0,    -1/4,  1/2])
    z_max_xp = 2.0 * np.array([1/4,  0,     1/2])
    z_max_xm = 2.0 * np.array([-1/4, 0,     1/2])
    z_min_yp = 2.0 * np.array([0,    1/4,   -1/2])
    z_min_ym = 2.0 * np.array([0,    -1/4,  -1/2])
    z_min_xp = 2.0 * np.array([1/4,  0,     -1/2])
    z_min_xm = 2.0 * np.array([-1/4, 0,     -1/2])
    all_lines = [
        [y_max_zm, y_max_xm],  # ymax square
        [y_max_xm, y_max_zp],
        [y_max_zp, y_max_xp],
        [y_max_xp, y_max_zm],
        [z_max_ym, z_max_xm],  # zmax square
        [z_max_xm, z_max_yp],
        [z_max_yp, z_max_xp],
        [z_max_xp, z_max_ym],
        [x_max_ym, x_max_zm],  # xmax square
        [x_max_zm, x_max_yp],
        [x_max_yp, x_max_zp],
        [x_max_zp, x_max_ym],
        [y_min_zm, y_min_xm],  # ymin square
        [y_min_xm, y_min_zp],
        [y_min_zp, y_min_xp],
        [y_min_xp, y_min_zm],
        [z_min_ym, z_min_xm],  # zmin square
        [z_min_xm, z_min_yp],
        [z_min_yp, z_min_xp],
        [z_min_xp, z_min_ym],
        [x_min_ym, x_min_zm],  # xmin square
        [x_min_zm, x_min_yp],
        [x_min_yp, x_min_zp],
        [x_min_zp, x_min_ym],
        [x_max_zp, z_max_xp],  # Here start the link between square, from xmax
        [x_max_zm, z_min_xp],
        [x_max_yp, y_max_xp],
        [x_max_ym, y_min_xp],
        [x_min_zp, z_max_xm],  # From xmin
        [x_min_zm, z_min_xm],
        [x_min_yp, y_max_xm],
        [x_min_ym, y_min_xm],
        [z_max_yp, y_max_zp],  # The others
        [z_max_ym, y_min_zp],
        [z_min_yp, y_max_zm],
        [z_min_ym, y_min_zm]
    ]
    for lines in all_lines:
        plot_segment(lines[0], lines[1], ax)



def RandomKPointsIrreducibleWedge(Npoints, plot=False):
    kpointsBZ = []
    kpoints_notBZ = []
    while len(kpointsBZ) < Npoints:
        k = np.random.uniform(-1.0, 1.0, 3)
        if IsInIrreducibleWedge(k):
            kpointsBZ.append(k)
        else:
            kpoints_notBZ.append(k)
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        plot_bz_outline(ax)
        ax.scatter([k[0] for k in kpointsBZ], [k[1]
                for k in kpointsBZ], [k[2] for k in kpointsBZ], c='r', marker='.')
        # ax.scatter([k[0] for k in kpoints_notBZ], [k[1]
        #            for k in kpoints_notBZ], [k[2] for k in kpoints_notBZ], c='b', alpha=0.25)
        plt.show()
    return kpointsBZ

def comparison(Nk):
    k_random = RandomKPointsIrreducibleWedge(Nk, plot=False)
    k_Monkhorst = Monkhorst_Pack(10, 10, 10)
    k_Monkhorst = [x for x in k_Monkhorst if IsInIrreducibleWedge(x)]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plot_bz_outline(ax)
    # ax.scatter([k[0] for k in k_random], [k[1]
    #         for k in k_random], [k[2] for k in k_random], c='r', marker='.', alpha=0.85)
    ax.scatter([k[0] for k in k_Monkhorst], [k[1]
            for k in k_Monkhorst], [k[2] for k in k_Monkhorst], c='b', marker='+')
    print(f"Number of kpoints in the irreducible wedge Mon: {len(k_Monkhorst)}")
    plt.show()


if __name__ == "__main__":
    Nk1 = 8
    Nk2 = 8
    Nk3 = 8
    list_k = Monkhorst_Pack(Nk1, Nk2, Nk3)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plot_bz_outline(ax)
    ax.scatter(list_k[:, 0], list_k[:, 1], list_k[:, 2], c='r', s=4)
    plt.show()
    # RandomKPointsIrreducibleWedge(5000, plot=True)
    # comparison(150)
