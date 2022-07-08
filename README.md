# EmpiricalPseudopotential
Band structure computation using empirical pseudopotentials


This repository was initially a fork from : [EmpiricalPseudopotential](https://github.com/aromanro/EmpiricalPseudopotential) of [Adrian Roman](https://compphys.go.ro/empirical-pseudopotential/).

---

## Features
You can do two types of calculations:
* Compute the electronic band structure over a path of high-symmetry points (e.g. $L\Gamma X U K \Gamma$) for a given material, and plot the results :  
![Silicon Band Structure over LGXUG path](doc/EEP_Si_nb_bands_16_path_LGXUG_size_basis_181.png "SiliconBandStructure")

---

* Compute the electronic band structure over all k-points of an input mesh of the Brillouin Zone (or a fraction of it). The result can then be visualized, for example, through iso-energy surface:
![Animation of iso-energy surfaces of the 1st conduction band of Silicon](doc/rotation_animation_4th_band_iso.gif "Silicon1stCB_isoenergy")

---

## Compilation
### Dependencies  
This project relies on these libraries:
* Eigen
* GMSH
* OpenMP (optional)
* tclap (embedded with the project)


You will also need standard packages such as CMake, Make, a C++ compiler, etc.
A minimal installation command would look like:  
`sudo apt-get update && sudo apt-get install -y apt-utils git curl cmake g++ 
libopenmpi-dev`  

If you don't have Eigen install on your system, the sources will be automatically 
fetched by CMake when you'll compile the project. You don't need to do anything. 

You can install GMSH from sources straightforwardly :   
Go wherever you want to install the library and run:   
`git clone https://gitlab.onelab.info/gmsh/gmsh.git && cd gmsh && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_BUILD_SHARED=1 -DENABLE_PRIVATE_API=1 .. && make -j 8 shared && make install/fast && cd .. && rm -rf gmsh `  
Find more informations on [GMSH Website](https://gmsh.info/)


### Compilation  
1. `git clone https://github.com/RemiHelleboid/EmpiricalPseudopotential.git`
2. `cd EmpiricalPseudopotential`
3. `mkdir build && cd build`
4. `cmake ..`
5. `make`
   
---
## Usage


## Available Material




