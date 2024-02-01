import gmsh
import sys
import numpy as np

gmsh.initialize(sys.argv)

gmsh.model.add("boolean")

# from http://en.wikipedia.org/wiki/Constructive_solid_geometry

gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.05)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.10)

Rt = 1.0

R = 3.0
Rs = R * .7
Rt = R 

gmsh.model.occ.addBox(-R, -R, -R, 2 * R, 2 * R, R, 1)
gmsh.model.occ.addSphere(0, 0, 0, Rt, 2)
gmsh.model.occ.intersect([(3, 1)], [(3, 2)], 3)

gmsh.model.occ.addBox(-R, -R, -R, 2 * R,  R, 2*R, 4)
gmsh.model.occ.intersect([(3, 3)], [(3, 4)], 5)

gmsh.model.occ.addBox(-R, -R, -R,  R, 2* R, 2*R, 6)
gmsh.model.occ.intersect([(3, 5)], [(3, 6)], 7)


gmsh.model.occ.synchronize()

gmsh.model.mesh.generate(3)



# Get list of nodes
nodes = gmsh.model.mesh.getNodes()
Coordinates = nodes[1].reshape((nodes[0].shape[0], 3))

filename = f"nodes_eight_sphere_{Rt}.csv"
np.savetxt(filename, Coordinates, delimiter=" ")

gmsh.write(f"EightSphere_{R}.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
