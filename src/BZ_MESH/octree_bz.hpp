/**
 * @file Octree_mesh.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-08-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "bbox_mesh.hpp"
#include "mesh_tetra.hpp"
#include "vector.hpp"

namespace uepm::mesh_bz {

class Octree_mesh {
 protected:
    /**
     * @brief Maximum number of elements in a leaf node.
     * If the number of elements in a node is greater than this value, the node is split into 8 children.
     */
    static constexpr std::size_t max_number_of_elements_per_node = 32;

    /**
     * @brief Minimum size of a node.
     * If the size of a node is smaller than this value, the node is not split anymore.
     */
    static constexpr double min_size_of_a_node = 1e8;

    /**
     * @brief Bounding box of the node.
     */
    bbox_mesh m_node_box;

    /**
     * @brief A leaf node is a node that has no children.
     *
     */
    bool m_is_leaf = false;

    /**
     * @brief List of pointers to elements overlapping in the node.
     */
    std::vector<Tetra *> m_list_tetras;

    /**
     * @brief List of children of the node.
     */
    std::vector<std::unique_ptr<Octree_mesh>> m_list_sub_nodes;

    std::vector<Tetra *> find_overlapping_tetras(const std::vector<Tetra *> &list_p_tetras, const bbox_mesh &bounding_box);
    bool                 is_inside(const vector3 &location) const { return m_node_box.is_inside(location); }

 public:
    Octree_mesh(){};
    Octree_mesh(const std::vector<Tetra *> &list_tetras, const bbox_mesh &bounding_box);

    bool is_leaf() const { return m_is_leaf; }
    const bbox_mesh &get_bounding_box() const { return m_node_box; }
    const std::vector<Tetra *> &get_list_tetras() const { return m_list_tetras; }
    const std::vector<std::unique_ptr<Octree_mesh>> &get_list_sub_nodes() const { return m_list_sub_nodes; }

    Tetra *find_tetra_at_location(const vector3 &location) const;
};

}  // namespace uepm::mesh_bz
