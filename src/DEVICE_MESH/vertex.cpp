/**
 * @file vertex.cpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief
 * @version 0.1
 * @date 2021-07-10
 *
 * @copyright Copyright (c) 2021
 *
 */
#include "vertex.hpp"

namespace uepm {

namespace mesh {

void vertex::print_info() const {
    std::cout << "Index    : " << m_index << std::endl;
    std::cout << "X        : " << x() << std::endl;
    std::cout << "Y        : " << y() << std::endl;
    std::cout << "Z        : " << z() << std::endl;
}

}  // namespace mesh

}  // namespace uepm