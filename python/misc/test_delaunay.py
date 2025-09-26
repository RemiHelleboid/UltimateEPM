import numpy as np
from scipy.spatial import Delaunay
import gmsh
import sys, os
import itertools


# Generate random 3D points
# points = np.random.rand(30, 3)
            # Gamma = gmsh.model.occ.addPoint(0.0, 0.0, 0.0, h * gamma_mesh) # Γ
            # L = gmsh.model.occ.addPoint(0.5, 0.5, 0.5, h)     # L
            # X = gmsh.model.occ.addPoint(1.0, 0.0, 0.0, h)     # X
            # K = gmsh.model.occ.addPoint(0.75, 0.75, 0.0, h)   # K
            # W = gmsh.model.occ.addPoint(1.0, 0.5, 0.0, h)     # W
            # U = gmsh.model.occ.addPoint(1.0, 0.25, 0.25, h)   # U


def write_msh(filename, points, tetrahedra):
    """
    Write a Gmsh 2.2 ASCII .msh file with 3D tetrahedral elements.
    
    Parameters
    ----------
    filename : str
        Output filename (e.g. "mesh.msh").
    points : ndarray (N, 3)
        Array of 3D points (x,y,z).
    tetrahedra : ndarray (M, 4)
        Array of tetrahedral connectivity (node indices, 0-based).
    """
    with open(filename, "w") as f:
        # Header
        f.write("$MeshFormat\n")
        f.write("2.2 0 8\n")  # version, file-type (0=ASCII), data-size
        f.write("$EndMeshFormat\n")

        # Nodes
        f.write("$Nodes\n")
        f.write(f"{len(points)}\n")
        for i, p in enumerate(points, start=1):
            f.write(f"{i} {p[0]} {p[1]} {p[2]}\n")
        f.write("$EndNodes\n")

        # Elements
        f.write("$Elements\n")
        f.write(f"{len(tetrahedra)}\n")
        for i, tet in enumerate(tetrahedra, start=1):
            # Element type 4 = Tetrahedron (4-node)
            # No tags (we just put 0)
            f.write(f"{i} 4 0 {tet[0]+1} {tet[1]+1} {tet[2]+1} {tet[3]+1}\n")
        f.write("$EndElements\n")

        # 8 symetry and 6 permutations
def _permutation_matrices():
    mats = []
    for perm in itertools.permutations(range(3)):
        P = np.zeros((3, 3), dtype=int)
        for i, j in enumerate(perm):
            P[i, j] = 1
        mats.append(P)
    return mats  # 6

def _reflection_matrices():
    mats = []
    for sx, sy, sz in itertools.product([-1, 1], repeat=3):
        mats.append(np.diag([sx, sy, sz]))
    return mats  # 8

def _symmetry_ops_full():
    ops, seen = [], set()
    for P in _permutation_matrices():
        for R in _reflection_matrices():    
            M = R @ P
            key = tuple(M.flatten())
            if key not in seen:
                seen.add(key)
                ops.append(M)
    assert len(ops) == 48
    for op in ops:
        print(op)
        print()
    return ops


# Example usage:
if __name__ == "__main__":
    import numpy as np
    from scipy.spatial import Delaunay

    # Random 3D points
    points = np.array([
    [0.0, 0.0, 0.0],   # Gamma
    [0.5, 0.5, 0.5],   # L
    [1.0, 0.0, 0.0],   # X
    [0.75, 0.75, 0.0], # K
    [1.0, 0.5, 0.0],   # W
    [1.0, 0.25, 0.25]  # U
])
    # Delaunay tetrahedralization
    tri = Delaunay(points)

    # Write to .msh
    write_msh("delaunay_mesh.msh", points, tri.simplices)

    # Load and visualize in Gmsh
    gmsh.initialize()
    gmsh.open("delaunay_mesh.msh")

    nb_levels = sys.argv[1] if len(sys.argv) > 1 else 0
    nb_levels = int(nb_levels)
    for i in range(nb_levels):
        gmsh.model.mesh.refine()

    # Get nodes and elements
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements()
    print(f"Nodes: {len(node_tags)}, Elements: {len(elem_tags[0])}")
    gmsh.finalize()

    # Save the coords in a csv file
    np.savetxt("delaunay_mesh_nodes.csv", node_coords.reshape(-1, 3), delimiter=",")

    points = node_coords.reshape(-1, 3)

    List_new_nodes = []
    all_ops = _symmetry_ops_full()
    print(f"Number of symmetry operations: {len(all_ops)}")
    for op in all_ops:
        for p in points:
            p_new = op @ p
            # if np.all(p_new >= -1e-8) and np.all(p_new <= 1 + 1e-8):
            List_new_nodes.append(tuple(np.round(p_new, 8)))
    List_new_nodes = list(set(List_new_nodes))
    print(f"Number of unique new nodes: {len(List_new_nodes)}")
    np.savetxt("delaunay_mesh_nodes_full_symmetry.csv", np.array(List_new_nodes), delimiter=",")

    print("List of new nodes:")
    print(List_new_nodes)

    list_coords = []
    for node in List_new_nodes:
        list_coords.append(node[0])
        list_coords.append(node[1])
        list_coords.append(node[2])
    
    gmsh.initialize()
    gmsh.model.add("delaunay_mesh_refined")

    gmsh.option.setNumber('Mesh.Algorithm3D', 10) # new algo
    N = len(List_new_nodes)
    points = np.array(list_coords)
    tets = gmsh.algorithm.tetrahedralize(points)
    vol = gmsh.model.addDiscreteEntity(3)
    gmsh.model.mesh.addNodes(3, vol, range(1, N+1), points)
    gmsh.model.mesh.addElementsByType(vol, 4, [], tets)

    # save the mesh
    gmsh.write("delaunay_mesh_refined.msh")

    gmsh.fltk.run()

    gmsh.finalize()


