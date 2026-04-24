/**
 * @file dataset.cpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief
 * @version 0.1
 * @date 2021-08-17
 *
 * @copyright Copyright (c) 2021
 *
 */
#include "dataset.hpp"
#include "export_vector_to_csv.hpp"


#include <cassert>

namespace uepm {

namespace mesh {

std::ostream &operator<<(std::ostream &os, const DataType &type) {
    if (type == DataType::scalar) {
        os << "scalar";
    } else if (type == DataType::vector) {
        os << "vector";
    } else if (type == DataType::tensor) {
        os << "tensor";
    } else {
        assert("Error: operator<< overloaded for DataType does not support the type.");
    }
    return os;
}
template <>
void dataset<double>::print_info() const {
    std::cout << "\nDataset name                  : " << m_name << std::endl;
    std::cout << "Dataset index                 : " << m_index << std::endl;
    std::cout << "Dataset index region validity : " << m_index_region_validity << std::endl;
    std::cout << "Dataset data type             : " << m_data_type << std::endl;
    std::cout << "Dataset data dimension        : " << m_dimension << std::endl;
    std::cout << "Dataset number values         : " << m_values->size() << std::endl;
}

template <>
void dataset<double>::dump_values_in_file(const std::string &filename) const{
    utils::export_vector_to_csv(filename, m_name, *m_values);
}


}  //  namespace mesh

}  // namespace uepm