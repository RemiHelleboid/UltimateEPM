import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm


def scalar_triple_product(v1, v2, v3):
    A = v1.get_coords_as_array()
    B = v2.get_coords_as_array()
    C = v3.get_coords_as_array()
    return np.dot(A, np.cross(B, C))


class vertex:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.energy = 0

    def __str__(self):
        return "(" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ")"

    def __repr__(self):
        return "(" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ', ' + str(self.energy) +")"

    def set_energy(self, energy):
        self.energy = energy

    def randomize_vtx(self, mean, std, randomize_energy, min=None, max=None):
        self.x = np.random.normal(mean, std)
        self.y = np.random.normal(mean, std)
        self.z = np.random.normal(mean, std)
        if min is not None and max is not None:
            self.x = np.clip(self.x, min, max)
            self.y = np.clip(self.y, min, max)
            self.z = np.clip(self.z, min, max)
        if randomize_energy:
            self.energy = np.random.normal(mean, std) ** 2
            self.energy = np.clip(self.energy, 0, 1e100)

    def get_coords_as_array(self):
        return np.array([self.x, self.y, self.z])

    def __sub__(self, other):
        return vertex(self.x - other.x, self.y - other.y, self.z - other.z)
    
    def __add__(self, other):
        return vertex(self.x + other.x, self.y + other.y, self.z + other.z)


class Triangle:
    def __init__(self, v1, v2, v3):
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.energy = 0
        self.list_vtx = [self.v1, self.v2, self.v3]

    def __str__(self):
        return "(" + str(self.v1) + ", " + str(self.v2) + ", " + str(self.v3) + ")"

    def __repr__(self):
        return "(" + str(self.v1) + ", " + str(self.v2) + ", " + str(self.v3) + ")"

    def set_energy(self, energy):
        self.energy = energy
        
    def get_list_vtx(self):
        return [self.v1, self.v2, self.v3]

    def randomize_triangle(self, mean, std, randomize_energy, min=None, max=None):
        self.v1.randomize_vtx(mean, std, randomize_energy, min, max)
        self.v2.randomize_vtx(mean, std, randomize_energy, min, max)
        self.v3.randomize_vtx(mean, std, randomize_energy, min, max)
        if randomize_energy:
            self.energy = np.random.normal(mean, std) ** 2
            self.energy = np.clip(self.energy, 0, 1e100)

    def get_coords_as_array(self):
        return np.array([self.v1.get_coords_as_array(), self.v2.get_coords_as_array(), self.v3.get_coords_as_array()])

    def get_energy(self):
        return self.energy

    def get_vertices(self):
        return [self.v1, self.v2, self.v3]

    def get_vertex_by_index(self, index):
        if index == 0:
            return self.v1
        elif index == 1:
            return self.v2
        elif index == 2:
            return self.v3
        else:
            return None

    def get_edges(self):
        res = list(combinations(self.list_vtx, 2))
        return res

    def commpute_surface_triangle(self):
        return np.dot(self.list_vtx[0] - self.list_vtx[1], self.list_vtx[0] - self.list_vtx[2]) / 2.0

    def compute_barycentric_coordinates(self, point):
        v0 = self.v1 - point
        v1 = self.v2 - point
        v2 = self.v3 - point
        d00 = np.dot(v0, v0)
        d01 = np.dot(v0, v1)
        d11 = np.dot(v1, v1)
        d20 = np.dot(v2, v0)
        d21 = np.dot(v2, v1)
        denom = d00 * d11 - d01 * d01
        v = (d11 * d20 - d01 * d21) / denom
        w = (d00 * d21 - d01 * d20) / denom
        u = 1.0 - v - w
        return u, v, w

    def compute_barycenter(self):
        return (self.v1 + self.v2 + self.v3) / 3.0


    def plot_triangle(self, ax):
        list_edges = self.get_edges()
        for edge in list_edges:
            ax.plot([edge[0].x, edge[1].x], [edge[0].y, edge[1].y], [edge[0].z, edge[1].z], color='blue',
                    lw=1.0e-1)


class Tetra:
    def __init__(self, v1, v2, v3, v4):
        self.list_vtx = [v1, v2, v3, v4]

    def unitary_tetra(self):
        self.list_vtx[0] = vertex(0, 0, 0)
        self.list_vtx[1] = vertex(1, 0, 0)
        self.list_vtx[2] = vertex(0, 1, 0)
        self.list_vtx[3] = vertex(0, 0, 1)

    def randomize_tetra(self, mean, std, randomize_energies=False, min=None, max=None):
        for vtx in self.list_vtx:
            vtx.randomize_vtx(mean, std, randomize_energies, min, max)

    def __str__(self):
        return "(" + str(self.vA) + ", " + str(self.vB) + ", " + str(self.vC) + ", " + str(self.vD) + ")"

    def __repr__(self):
        return "(" + str(self.vA) + ", " + str(self.vB) + ", " + str(self.vC) + ", " + str(self.vD) + ")"

    def set_values_energy_from_list(self, values):
        self.values_energy = values

    def get_values_energy(self):
        return [vtx.energy for vtx in self.list_vtx]

    def set_values_energy(self, e0, e1, e2, e3):
        self.list_vtx[0].set_energy(e0)
        self.list_vtx[1].set_energy(e1)
        self.list_vtx[2].set_energy(e2)
        self.list_vtx[3].set_energy(e3)

    def sort_vertices_order_from_energy(self):
        print(self.list_vtx)
        self.list_vtx.sort(key=lambda vtx: vtx.energy)
        print(f"New order: {self.list_vtx}")

    def get_edges(self):
        res = list(combinations(self.list_vtx, 2))
        return res

    def compute_tetra_volume(self):
        volume = scalar_triple_product(self.list_vtx[1] - self.list_vtx[0], self.list_vtx[2] - self.list_vtx[0],
                                       self.list_vtx[3] - self.list_vtx[0])
        return volume / 6.0

    def plot_tetra(self, ax, show=False):
        list_edges = self.get_edges()
        for edge in list_edges:
            ax.plot([edge[0].x, edge[1].x], [edge[0].y, edge[1].y], [edge[0].z, edge[1].z], color='black',
                    lw=1.0e-1)
        if show:
            plt.show()

    def plot_tetra_with_energy(self, ax, show=False):
        list_edges = self.get_edges()
        list_energies = self.get_values_energy()
        print(f"List energies: {list_energies}")
        print(list_energies)
        for edge in list_edges:
            ax.plot([edge[0].x, edge[1].x], [edge[0].y, edge[1].y], [edge[0].z, edge[1].z], color='black',
                    lw=1.0)
        for idx, vtx in enumerate(self.list_vtx):
            ax.scatter(vtx.x, vtx.y, vtx.z, c=list_energies[idx], s=200.0, cmap=cm.get_cmap(
                "jet"), vmin=0, vmax=np.max(list_energies))
        if show:
            plt.show()

    def compute_iso_surface(self, iso_value):
        self.sort_vertices_order_from_energy()
        list_vertices_iso = []
        eps_0 = iso_value

        eps_01 = eps_0 - self.list_vtx[0].energy
        eps_20 = self.list_vtx[1].energy - eps_0
        eps_21 = self.list_vtx[1].energy - self.list_vtx[0].energy
        eps_30 = self.list_vtx[2].energy - eps_0
        eps_31 = self.list_vtx[2].energy - self.list_vtx[0].energy
        eps_32 = self.list_vtx[2].energy - self.list_vtx[1].energy
        eps_40 = self.list_vtx[3].energy - eps_0
        eps_41 = self.list_vtx[3].energy - self.list_vtx[0].energy
        eps_42 = self.list_vtx[3].energy - self.list_vtx[1].energy
        eps_43 = self.list_vtx[3].energy - self.list_vtx[2].energy

        # The vertices are order in such a way that the first vertex is the one with the lowest energy and the energies are increasing.

        # CASE 1: If the energy of the iso_value is lower than the energy of the first vertex, the iso_value is not in the tetrahedron.
        if (iso_value > self.list_vtx[3].energy):
            print("CASE 1")
            return []

        # CASE 2:  If the energy of the iso_value is higher than the energy of the last vertex, the iso_value is not in the tetrahedron.
        if (iso_value < self.list_vtx[0].energy):
            print("CASE 2")
            return []

        # CASE 3: If the energy of the iso_value is between the energy of the first and the last vertex, the iso_value is in the tetrahedron.
        # In this case, the iso surface is a tringle made of vertices T1 T2 T3.
        if (iso_value < self.list_vtx[1].energy and iso_value > self.list_vtx[0].energy):
            print("CASE 3")
            T1 = vertex(eps_01 / eps_21, 0, 0)
            T2 = vertex(0, eps_01 / eps_31, 0)
            T3 = vertex(0, 0, eps_01 / eps_41)
            return [T1, T2, T3]

        # CASE 4: If the energy of the iso_value is between the energy of the second and the third vertex, the iso_value is in made of 4 vertices:
        # T1 T2 T3 T4.
        if iso_value < self.list_vtx[2].energy and iso_value > self.list_vtx[1].energy:
            print("CASE 4")
            T1 = vertex(eps_40/eps_42, 0, - eps_20/eps_42)
            T2 = vertex(eps_30 / eps_32, -eps_20 / eps_32, 0)
            T3 = vertex(0.0, eps_01/eps_31, 0.0)
            T4 = vertex(0.0, 0.0, eps_01/eps_41)
            return [T1, T2, T3, T4]

        # CASE 5: If the energy of the iso_value is between the energy of the third and the fourth vertex, the iso_value is in made of 3 vertices:
        # T1 T2 T3.
        if iso_value < self.list_vtx[3].energy and iso_value > self.list_vtx[2].energy:
            print("CASE 5")
            T1 = vertex(eps_40/eps_42, 0, - eps_20/eps_42)
            T2 = vertex(0.0, eps_40 / eps_43, -eps_30/eps_43)
            T3 = vertex(0.0, 0.0, eps_01/eps_41)
            return [T1, T2, T3]

        else:
            print("Error: The iso_value is not in the tetrahedron.")
            return []


def get_unit_tetra():
    v1 = vertex(0, 0, 0)
    v2 = vertex(0, 0, 0)
    v3 = vertex(0, 0, 0)
    v4 = vertex(0, 0, 0)
    new_tetra = Tetra(v1, v2, v3, v4)
    new_tetra.unitary_tetra()
    return new_tetra


def generates_random_tetra(nb_tetra):
    list_tetras = []
    for index_tetra in range(nb_tetra):
        v1 = vertex(0, 0, 0)
        v2 = vertex(0, 0, 0)
        v3 = vertex(0, 0, 0)
        v4 = vertex(0, 0, 0)
        new_tetra = Tetra(v1, v2, v3, v4)
        new_tetra.randomize_tetra(0, 1, True)
        # new_tetra.unitary_tetra()
        list_tetras.append(new_tetra)

    # 3D plot of tetras
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, projection='3d')

    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.set_axis_off()
    ax.set_rasterization_zorder(4)

    bound = 2
    ax.set_xlim([-bound, bound])
    ax.set_ylim([-bound, bound])
    ax.set_zlim([-bound, bound])
    ax.set_aspect('auto')

    for tetra in list_tetras:
        tetra.plot_tetra_with_energy(ax)

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    UNIT_TETRA = get_unit_tetra()

    UNIT_TETRA.set_values_energy(0.499, 0.05, 0.10, 0.9900)

    print(
        f"Volume of the unitary tetrahedron: {UNIT_TETRA.compute_tetra_volume()}")
    energy_iso = 0.5
    
    ISO_TRIANGLE = UNIT_TETRA.compute_iso_surface(energy_iso)
    ISO_SURF = Triangle(ISO_TRIANGLE[0], ISO_TRIANGLE[1], ISO_TRIANGLE[2])
    
    print(f"Iso surface of the tetrahedron: {ISO_SURF}")
    

    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, projection='3d')
    
    UNIT_TETRA.plot_tetra_with_energy(ax)
    ISO_SURF.plot_triangle(ax)
    
    plt.show()