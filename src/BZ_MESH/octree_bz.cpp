/**
 * @file Octree_mesh.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-08-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "octree_bz.hpp"


#include <algorithm>
#include <array>
#include <memory>
#include <vector>

namespace bz_mesh {

Octree_mesh::Octree_mesh(const std::vector<Tetra *> &list_tetras, const bbox_mesh &bounding_box) {
    if (list_tetras.size() <= max_number_of_elements_per_node || bounding_box.get_diagonal_size() < min_size_of_a_node) {
        m_is_leaf     = true;
        m_list_tetras = list_tetras;
        m_node_box    = bounding_box;
        return;
    }

    m_is_leaf                        = false;
    m_node_box                       = bounding_box;
    std::array<bbox_mesh, 8> sub_boxes = m_node_box.split_3d_box_in_octants();

    for (int i = 0; i < 8; i++) {
        m_list_sub_nodes.push_back(
            std::make_unique<Octree_mesh>(Octree_mesh::find_overlapping_tetras(list_tetras, sub_boxes[i]), sub_boxes[i]));
    }
}

/**
 * @brief Find the tetra from the list of tetras that are overlapping with the bounding box.
 *
 * @note The list of pointers to tetras is passed by copy because it is modified during the search.
 * Perf should be compared with a list of pointers to tetras that is passed by reference and copied within the function.
 *
 * @param list_p_tetras
 * @param bounding_box
 * @return std::vector<Tetra *>
 */
std::vector<Tetra *> Octree_mesh::find_overlapping_tetras(const std::vector<Tetra *> &list_p_tetras, const bbox_mesh &bounding_box) {
    std::vector<Tetra *> list_overlapping_tetras;
    for (auto &p_tetra : list_p_tetras) {
        if (bounding_box.is_overlapping(p_tetra->get_bounding_box())) {
            list_overlapping_tetras.push_back(p_tetra);
        }
    }
    return list_overlapping_tetras;
}

/**
 * @brief Find the tetra that contains the location.
 * If none is found, return nullptr, else return a pointer to the tetra.
 *
 * @param location
 * @return Tetra*
 */
Tetra *Octree_mesh::find_tetra_at_location(const vector3 &location) const {
    if (!m_is_leaf) {
        for (auto &&p_sub_node : m_list_sub_nodes) {
            if (p_sub_node->is_inside(location)) {
                return p_sub_node->find_tetra_at_location(location);
            }
        }
    }
    const auto it_tetra =
        std::find_if(m_list_tetras.begin(), m_list_tetras.end(), [&](const Tetra *p_tetra) { return p_tetra->is_location_inside(location); });
    return it_tetra != m_list_tetras.end() ? *it_tetra : nullptr;
}

}  // namespace bz_mesh
