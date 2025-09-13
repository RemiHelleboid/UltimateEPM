import gmsh
import sys
import numpy as np

gmsh.initialize(sys.argv)

gmsh.model.add("boolean")

# from http://en.wikipedia.org/wiki/Constructive_solid_geometry

gmsh.option.setNumber("Mesh.Algorithm", 6)

minSize = 0.2
maxSize = 0.25
gmsh.option.setNumber("Mesh.MeshSizeMin", minSize)
gmsh.option.setNumber("Mesh.MeshSizeMax", maxSize)

Rt = 1.0

R = 3.0
Rs = R * .7
Rt = R 

gmsh.model.occ.addBox(0, 0, 0, 2 * R, 2 * R, 2 * R, 1)
gmsh.model.occ.addSphere(0, 0, 0, Rt, 2)
gmsh.model.occ.intersect([(3, 1)], [(3, 2)], 3)

# gmsh.model.occ.addBox(-R, -R, -R, 2 * R,  R, 2*R, 4)
# gmsh.model.occ.intersect([(3, 3)], [(3, 4)], 5)

# gmsh.model.occ.addBox(-R, -R, -R,  R, 2* R, 2*R, 6)
# gmsh.model.occ.intersect([(3, 5)], [(3, 6)], 7)


gmsh.model.occ.synchronize()
gmsh.model.occ.removeAllDuplicates()

gmsh.model.mesh.generate(3)



# Get list of nodes
nodes = gmsh.model.mesh.getNodes()
Coordinates = nodes[1].reshape((nodes[0].shape[0], 3))

filename = f"nodes_eight_sphere_{R}.csv"
np.savetxt(filename, Coordinates, delimiter=" ")

gmsh.write(f"EightSphere_{R}.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()


    dict_alg = {"1": "Delaunay", 
                "3": "Initial mesh only",
                "4": "Frontal",
                "7": "MMG3D",
                "9": "R-tree",
                "10": "HXT"}
    # 1 ok; 3 bad; 4 weird; 7 bad; 9 okish; 10 ok.
