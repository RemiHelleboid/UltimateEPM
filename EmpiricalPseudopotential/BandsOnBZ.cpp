/**
 * @file compute_bands_on_bz_mesh.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2022-07-07
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "BandStructure.h"
#include "Options.h"
#include "Material.h"
#include "bz_mesh.hpp"

int main(int argc, char* argv[]) {
    EmpiricalPseudopotential::Materials materials;

    Options my_options;
    my_options.print_options();
    EmpiricalPseudopotential::Material mat = materials.materials.at(my_options.materialName);

    // bz_mesh my_mesh("mesh.msh");
    bz_mesh my_mesh("mymesh_1.msh");
    my_mesh.read_mesh();
    std::vector<Vector3D<double>>& mesh_kpoints = my_mesh.get_kpoints();
    
    EmpiricalPseudopotential::BandStructure my_bandstructure;
    my_bandstructure.Initialize(mat, my_options.nrLevels, mesh_kpoints, my_options.nearestNeighbors);
    my_bandstructure.Compute();
    my_bandstructure.AdjustValues();

    unsigned int band_index = 0;
    std::vector<double> band_0 = my_bandstructure.get_band(band_index);
    // export band in file
    my_bandstructure.export_result_in_file_with_kpoints("BZ_BANDS_SI.csv");

    // my_mesh.add_band_on_mesh("band_0", band_0);

    return 0;
}