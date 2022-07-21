# BzMeshBandsDos

Band structure and DOS computation using empirical pseudopotentials on the full Brillouin-Zone.

[![Build & Unit Test](https://github.com/RemiHelleboid/EmpiricalPseudopotential/actions/workflows/build_code.yaml/badge.svg)](https://github.com/RemiHelleboid/EmpiricalPseudopotential/actions/workflows/build_code.yaml)
[![Codacy Code Quality](https://app.codacy.com/project/badge/Grade/da70f725be754c928f4506a2bf86caea)](https://www.codacy.com/gh/RemiHelleboid/BzMeshBandsDos/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=RemiHelleboid/BzMeshBandsDos&amp;utm_campaign=Badge_Grade)

---

This repository was initially a fork from : [EmpiricalPseudopotential](https://github.com/aromanro/EmpiricalPseudopotential) of [Adrian Roman](https://compphys.go.ro/empirical-pseudopotential/).

---

## Features

You can do two types of calculations:

### Compute the electronic band structure over a path of high-symmetry points (e.g. $L\Gamma X U K \Gamma$) for a given material, and plot the results 
  ![Silicon Band Structure over LGXUG path](doc/EEP_Si_nb_bands_16_path_LGXUG_size_basis_181.png "SiliconBandStructure")

---

### Compute the electronic band structure over all k-points of an input mesh of the Brillouin Zone (or a fraction of it). The result can then be visualized, for example, through iso-energy surface
  ![Animation of iso-energy surfaces of the 1st conduction band of Silicon](doc/rotation_animation_4th_band_iso.gif "Silicon1stCB_isoenergy")

---

### Compute the density of states over the all Brillouin Zone
  ![Silicon Density of States For 10 Bands](doc/DOS_PER_BAND_DOS_bz_mesh_1_mediumEEP_Si_nb_bands_10_path__size_basis_137_all_bands.png "SiliconODS10Bands")
  ![Silicon Totoal Density of States](doc/DOS_TOTAL_DOS_bz_mesh_1_mediumEEP_Si_nb_bands_10_path__size_basis_137_all_bands.png "SiliconDOSTOTALBands")

## Build and Compilation

### Dependencies

This project relies on these libraries:

-   [Eigen](https://eigen.tuxfamily.org)
-   [GMSH](https://gmsh.info/)
-   [OpenMP](https://www.openmp.org/)
-   [tclap](http://tclap.sourceforge.net/)

You will also need standard packages such as CMake, Make, a C++ compiler, etc.
A **minimal installation** command would look like:  
`sudo apt-get update && sudo apt-get install -y apt-utils cmake g++ libopenmpi-dev`

If you don't have **Eigen** install on your system, the sources will be automatically fetched by CMake when you'll compile the project. You don't need to do anything. [Eigen Website](https://eigen.tuxfamily.org).

The **tclap** library is headers-only and embedded in the sources of the project.
You do not need to install it.

You can install **GMSH** from sources straightforwardly :  
Go wherever you want to install the library and run:  
`git clone https://gitlab.onelab.info/gmsh/gmsh.git && cd gmsh && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_BUILD_SHARED=1 -DENABLE_PRIVATE_API=1 .. && make -j 8 shared && make install/fast && cd .. && rm -rf gmsh `  
Find more informations on [GMSH Website](https://gmsh.info/).

### Compilation

*    `git clone https://github.com/RemiHelleboid/EmpiricalPseudopotential.git`  
*    `cd EmpiricalPseudopotential`  
*    `mkdir build && cd build`  
*    `cmake ..`  
*    `make`  

---

## Usage

## Available Material

## Brillouin Zone Meshing
To get the required mesh of the Brillouin Zone, you can use the BZ.py script from the great J. Grebot, there: [fcc-bz-mesh](https://github.com/JGrebot/fcc-bz-mesh).
  ![Diamond Brillouin Zone Mesh](doc/bz_mesh_jg_8.png "DiamondBZMesh")