/**
 * @file octree_node.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-31
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once

#include "tree_node.hpp"

namespace uepm {

namespace mesh {

class octree_node : public tree_node {
 private:
    std::vector<std::unique_ptr<octree_node>> m_sub_nodes;
    bbox                                      m_node_box;
    bbox                                      m_bottom_left_subbox_front;
    bbox                                      m_bottom_right_subbox_front;
    bbox                                      m_top_right_subbox_front;
    bbox                                      m_top_left_subbox_front;
    bbox                                      m_bottom_left_subbox_back;
    bbox                                      m_bottom_right_subbox_back;
    bbox                                      m_top_right_subbox_back;
    bbox                                      m_top_left_subbox_back;

    std::vector<element *> find_overlapping_elements(const std::vector<element *> &list_p_elements, const bbox &bounding_box) override;

 public:
    octree_node(){};
    octree_node(const std::vector<element *> &list_p_elements, const bbox &bounding_box, bool root);

    element *              find_element_at_location(const vector3 &position) const override;
    std::vector<element *> find_elements_overlapping_box(const bbox &box) const override;
    //  std::vector<element *> find_elements_overlapping_sphere(const vector3 sphere_center, double radius) const override;
};

}  //  namespace mesh

}  // namespace uepm