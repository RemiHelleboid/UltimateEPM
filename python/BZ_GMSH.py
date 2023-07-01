from argparse import ArgumentParser
import gmsh
import sys
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.spatial as spatial


gmsh.initialize()



# Define 8 reflection matrices
R1 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
R2 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
R3 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
R4 = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
R7 = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
R5 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])
R6 = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
R8 = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])

# Define 6 Permutation matrices
P1 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
P4 = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
P3 = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
P5 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
P2 = np.array([[0, 0, 1], [0, 1, 0], [0, 0, 1]])
P6 = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])

ListReflectionMatrices = [R1, R2, R3, R4, R5, R6, R7, R8]
ListPermutationMatrices = [P1, P2, P3, P4, P5, P6]

ListAllTransformations = []
for R in ListReflectionMatrices:
    for P in ListPermutationMatrices:
        ListAllTransformations.append(np.matmul(R, P))
TotalNumberTransformations = len(ListAllTransformations)

BZ_points = {
    "Ga": np.array([0, 0, 0]),
    "X":  np.array([1, 0, 0]),
    "L":  np.array([1/2, 1/2, 1/2]),
    "W":  np.array([1, 1/2, 0]),
    "U":  np.array([1, 1/4, 1/4]),
    "K":  np.array([3/4, 3/4, 0]),
}

def generate_all_points(initialPoint: np.array):
    """
    Generate all the points of the cube from the initial point
    :param initialPoint: np.array
    :return: list of np.array
    """
    listOfPoints = []
    k = 0
    for R in ListAllTransformations[:]:
        k += 1
        listOfPoints.append(np.matmul(R, initialPoint))
    # print(k)
    return listOfPoints



def createGeometryAndMesh():
    # Clear all models and create a new one
    gmsh.clear()
    gmsh.model.add("BZ")

    ListInitialPoints = BZ_points.values()
    PointIndex = {"Ga": 1, "X": 2, "L": 3, "W": 4, "U": 5, "K": 6}
    
    lc = 1e-2
    for p in ListInitialPoints:
        gmsh.model.geo.addPoint(p[0], p[1], p[2], lc)
    gmsh.model.geo.addLine(PointIndex["Ga"], PointIndex["X"], 1)
    gmsh.model.geo.addLine(PointIndex["X"], PointIndex["W"], 2)
    gmsh.model.geo.addLine(PointIndex["W"], PointIndex["K"], 3)
    gmsh.model.geo.addLine(PointIndex["K"], PointIndex["Ga"], 4)
    gmsh.model.geo.addLine(PointIndex["Ga"], PointIndex["L"], 5)
    gmsh.model.geo.addLine(PointIndex["L"], PointIndex["K"], 6)
    gmsh.model.geo.addLine(PointIndex["L"], PointIndex["W"], 7)
    gmsh.model.geo.addLine(PointIndex["L"], PointIndex["U"], 8)
    gmsh.model.geo.addLine(PointIndex["U"], PointIndex["X"], 9)
    gmsh.model.geo.addLine(PointIndex["U"], PointIndex["W"], 10)

    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
    gmsh.model.geo.addCurveLoop([5, 6, 4], 2)
    gmsh.model.geo.addCurveLoop([7, 3, -6], 3)
    gmsh.model.geo.addCurveLoop([8, 10, -7], 4)
    gmsh.model.geo.addCurveLoop([9, 2, -10], 5)
    gmsh.model.geo.addCurveLoop([1, -9, -8, -5], 6)

    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.addPlaneSurface([2], 2)
    gmsh.model.geo.addPlaneSurface([3], 3)
    gmsh.model.geo.addPlaneSurface([4], 4)
    gmsh.model.geo.addPlaneSurface([5], 5)
    gmsh.model.geo.addPlaneSurface([6], 6)

    # gmsh.model.geo.addSurfaceLoop([1, 2, 3, 4, 5, 6], 1)
    # gmsh.model.geo.addVolume([1], 1)

    # Physical groups
    gmsh.model.addPhysicalGroup(1, [1, 2, 3, 4, 5, 6], 1)
    # Mesh

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    # Get the vertices of the mesh
    nodes = gmsh.model.mesh.getNodes()
    nodes = np.array(nodes[1])
    print(nodes.shape)
    ListNodes = []
    for k in range(len(nodes)//3):
        ListNodes.append(np.array([nodes[3*k], nodes[3*k+1], nodes[3*k+2]]))
    return ListNodes

    
nodes = createGeometryAndMesh()
# exit(0)
AllNodes = []
for p in nodes:
    AllNodes += generate_all_points(p)
AllNodes = np.array(AllNodes)
X = AllNodes[:, 0]
Y = AllNodes[:, 1]
Z = AllNodes[:, 2]
np.savetxt("BZ_points2.csv", AllNodes, delimiter=",", header="X,Y,Z", comments="")



def checkForEvent():
    action = gmsh.onelab.getString("ONELAB/Action")
    if len(action) and action[0] == "check":
        gmsh.onelab.setString("ONELAB/Action", [""])
        createGeometryAndMesh()
        gmsh.graphics.draw()
    return True


if "-nopopup" not in sys.argv:
    gmsh.fltk.initialize()
    while gmsh.fltk.isAvailable() and checkForEvent():
        gmsh.fltk.wait()

