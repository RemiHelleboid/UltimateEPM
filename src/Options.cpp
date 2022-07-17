#include "Options.h"

#include <iostream>

Options::Options()
    : nrThreads(4),
      materialName("Si"),
      nrPoints(60),
      nearestNeighbors(12),
      nrLevels(30),
      pathNo(1),
      paths{{
          {"K", "W", "X", "G", "L", "W"},
          {"W", "G", "X", "W", "L", "G"},
          {"W", "L", "G", "X", "W", "K"},
          {"L", "G", "X", "K", "G"},
          {"L", "G", "X", "W", "K", "G"},
          {"L", "G", "X", "U", "K", "G"},
          {"L", "G", "X", "U", "G"},
          {"L", "K", "W", "G", "X", "W", "L", "G", "K"},
          {"G", "X", "K", "G", "L", "K", "W", "X"},
          {"G", "X", "W", "L", "G", "K", "W", "U"},
          {"G", "X", "W", "L", "G", "K"},
          {"G", "X", "W", "L", "G", "K", "X"},
          {"G", "X", "W", "L", "G", "X"},
          {"G", "X", "W", "G", "U", "X"},
          {"G", "X", "W", "K", "L", "G"},
          {"G", "X", "W", "K", "G", "L", "U", "W", "L", "K"},
          {"G", "X", "U", "K", "G", "L", "W", "X"},
      }},
      m_fileconfig(nullptr) {}

void Options::print_options() {
    std::cout << "Options:" << std::endl;
    std::cout << "nrThreads: " << nrThreads << std::endl;
    std::cout << "materialName: " << materialName << std::endl;
    std::cout << "nrPoints: " << nrPoints << std::endl;
    std::cout << "nearestNeighbors: " << nearestNeighbors << std::endl;
    std::cout << "nrLevels: " << nrLevels << std::endl;
    std::cout << "pathNo: " << pathNo << std::endl;
    // std::cout << "paths: " << std::endl;
    // for (auto& path : paths) {
    //     for (auto& point : path) {
    //         std::cout << point << " ";
    //     }
    //     std::cout << std::endl;
    // }
}


