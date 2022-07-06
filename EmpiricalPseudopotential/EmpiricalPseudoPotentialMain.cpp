#include "BandStructure.h"
#include "Options.h"

int main(int argc, char* argv[]) {
    EmpiricalPseudopotential::Materials materials;

    Options my_options;
    my_options.print_options();

    std::size_t                             path_index = 1;
    EmpiricalPseudopotential::BandStructure my_bandstructure;
    my_bandstructure.Initialize(materials.materials[my_options.materialName],
                                my_options.nrLevels,
                                my_options.paths[path_index],
                                my_options.nrPoints,
                                my_options.nearestNeighbors);

    const unsigned int nrPoints = my_bandstructure.GetPointsNumber();

    auto res = my_bandstructure.Compute();
    // Print res
    for (auto& row : res) {
        for (auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    my_bandstructure.AdjustValues();
    my_bandstructure.export_result_in_file(my_bandstructure.path_band_filename() + ".txt");

    return 0;
}