#include "BandStructure.h"
#include "Options.h"

int main(int argc, char* argv[]) {
    EmpiricalPseudopotential::Materials materials;

    Options my_options;
    my_options.print_options();

    std::size_t path_index = 1;

    for (auto const& [name, mat] : materials.materials) {
        std::cout << "----------------------------------------------------" << std::endl;
        std::cout << "Material: " << name << std::endl;
        for (int path_index = 0; path_index < my_options.paths.size(); path_index++) {
            std::vector<std::string> path = my_options.paths[path_index];
            std::cout << "path: ";
            for (auto& point : path) {
                std::cout << point << " ";
            }
            std::cout << std::endl;
            EmpiricalPseudopotential::BandStructure my_bandstructure;
            my_bandstructure.Initialize(mat,
                                        my_options.nrLevels,
                                        my_options.paths[path_index],
                                        my_options.nrPoints,
                                        my_options.nearestNeighbors);

            const unsigned int nrPoints = my_bandstructure.GetPointsNumber();

            // auto res = my_bandstructure.Compute_parralel(8);
            auto res = my_bandstructure.Compute();
            my_bandstructure.AdjustValues();
            my_bandstructure.export_result_in_file(my_bandstructure.path_band_filename() + ".txt");
        }
    }

    return 0;
}