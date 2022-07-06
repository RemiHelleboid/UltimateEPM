#include "Options.h"
#include "BandStructure.h"
#include "BandStructure.h"


int main(int argc, char* argv[])
{
    EmpiricalPseudopotential::Materials materials;

    Options my_options;
    my_options.print_options();

    std::size_t path_index = 6;
    EmpiricalPseudopotential::BandStructure my_bandstructure;
    my_bandstructure.Initialize(my_options.paths[path_index], my_options.nrPoints, my_options.nearestNeighbors);

    const unsigned int nrPoints = my_bandstructure.GetPointsNumber();
    
    auto res = my_bandstructure.Compute(materials.materials[my_options.materialName], 0, nrPoints, my_options.nrLevels);
    // Print res
    for (auto& row : res) {
        for (auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    my_bandstructure.AdjustValues();
    my_bandstructure.export_result_in_file("result_si.txt");

    
    return 0;
}