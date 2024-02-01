import gmsh
import sys
import numpy as np
import os, sys
import pathlib, shutil, glob
from pathlib import Path

gmsh.initialize(sys.argv)

filename = "EightSphere_3.0.msh"

gmsh.open(filename)
model_name = gmsh.model.getCurrent()

print(f"model_name = {model_name}")


def extract_dielectric_file(filename):
    energy, eps_r, eps_i = np.loadtxt(filename, delimiter=",", unpack=True, skiprows=1)
    return energy, eps_r, eps_i



DIR_OUT = "OUTDIR/"
list_dielectric_files = glob.glob(DIR_OUT + "Si*.csv")
print(f"list_dielectric_files = {list_dielectric_files[1]}")
list_idx_dielectric = [int(Path(file).stem.split("_")[2]) for file in list_dielectric_files]

NbFiles = len(list_dielectric_files)

# Get list of nodes
nodes = gmsh.model.mesh.getNodes()
TagNodes = nodes[0]
NbNodes = TagNodes.shape[0]

if NbFiles != NbNodes:
    raise ValueError(f"NbFiles = {NbFiles} != NbNodes = {NbNodes}")
else:
    print(f"NbFiles = {NbFiles} == NbNodes = {NbNodes}")

# Open the first file to get the energy
energy, eps_r, eps_i = extract_dielectric_file(list_dielectric_files[0])
print(f"energy = {energy}")
LIST_ENERGY = []
LIST_EPS_R = []
LIST_EPS_I = []

for file in list_dielectric_files:
    energy, eps_r, eps_i = extract_dielectric_file(file)
    LIST_ENERGY.append(energy)
    LIST_EPS_R.append(eps_r)
    LIST_EPS_I.append(eps_i)
    
LIST_ENERGY = np.array(LIST_ENERGY)
LIST_EPS_R = np.array(LIST_EPS_R)
LIST_EPS_I = np.array(LIST_EPS_I)

print(f"LIST_ENERGY = {LIST_ENERGY.shape}")
print(f"LIST_EPS_R = {LIST_EPS_R.shape}")
print(f"LIST_EPS_I = {LIST_EPS_I.shape}")

energies = LIST_ENERGY[0]
print(f"energies = {energies}")
NbEnergies = energies.shape[0]

FileName = "EightSphere_3.0_with_dielectric.msh"

for ix_e in range(NbEnergies//10):
    print(f"idx_e = {ix_e}")
    e = energies[ix_e]
    eps_r = LIST_EPS_R[:, ix_e]
    eps_i = LIST_EPS_I[:, ix_e]
    # Add view
    NameView = f"eps_r_{e}"
    data_tag = gmsh.view.add(NameView)
    gmsh.view.addHomogeneousModelData(data_tag, 0, model_name, "NodeData", TagNodes, eps_r)
    gmsh.view.write(data_tag, FileName, True)
    



# test_data = np.linalg.norm(Coordinates, axis=1)

# # Add view
# data_tag = gmsh.view.add("data_test")
# # gmsh::view::addHomogeneousModelData(data_tag, 0, "mesh_discrete", "NodeData", index_vtx,
# #                                                     vertex_index_values_of_function.second);
# gmsh.view.addHomogeneousModelData(data_tag, 0, model_name, "NodeData", TagNodes, test_data)

# # Export 
# # gmsh::view::write(data_tag, output_path, true);
# gmsh.view.write(data_tag, "test_data3.msh", True)

# test_data2 = np.sin(np.linalg.norm(Coordinates, axis=1))

# # Add view
# data_tag2 = gmsh.view.add("data_test_2")
# # gmsh::view::addHomogeneousModelData(data_tag, 0, "mesh_discrete", "NodeData", index_vtx,
# #                                                     vertex_index_values_of_function.second);
# gmsh.view.addHomogeneousModelData(data_tag2, 0, model_name, "NodeData", TagNodes, test_data2)

# # Export 
# # gmsh::view::write(data_tag, output_path, true);
# gmsh.view.write(data_tag2, "test_data3.msh", True)
