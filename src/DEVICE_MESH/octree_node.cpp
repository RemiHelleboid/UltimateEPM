/**
 * @file octree_node.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-31
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "octree_node.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

#include "bbox.hpp"
#include "element.hpp"
#include "element3d.hpp"

namespace uepm {

namespace mesh {

constexpr int    Max_Element_Per_Quadtree_Node = 32;
constexpr double Min_Box_Size                  = 1e-3;

static const std::string out_file = "OctreeSave.txt";

octree_node::octree_node(const std::vector<element *> &list_p_elements, const bbox &bounding_box, bool is_root) : m_node_box(bounding_box) {
    if (list_p_elements.size() <= Max_Element_Per_Quadtree_Node || m_node_box.get_diagonal_size() <= Min_Box_Size) {
        // std::ofstream my_file;
        // my_file.open(out_file, std::ios::app);
        // my_file << bounding_box << "," << list_p_elements.size() << std::endl;
        // my_file.close();
        this->m_list_p_elements = list_p_elements;
        m_is_leaf               = true;
        return;
    }
    const std::vector<bbox> list_sub_boxes = bounding_box.split_3d_box_in_octants();
    m_bottom_left_subbox_front             = list_sub_boxes[0];
    m_bottom_right_subbox_front            = list_sub_boxes[1];
    m_top_right_subbox_front               = list_sub_boxes[2];
    m_top_left_subbox_front                = list_sub_boxes[3];
    m_bottom_left_subbox_back              = list_sub_boxes[4];
    m_bottom_right_subbox_back             = list_sub_boxes[5];
    m_top_right_subbox_back                = list_sub_boxes[6];
    m_top_left_subbox_back                 = list_sub_boxes[7];

    constexpr uint number_octant = 8;
    m_sub_nodes.resize(number_octant);
    std::array<std::vector<element *>, number_octant> list_overlap_element_quadrants;
    const bool                                        child_is_root = false;
#pragma omp parallel for if (is_root)
    for (uint index_sub_box = 0; index_sub_box < number_octant; ++index_sub_box) {
        m_sub_nodes[index_sub_box] =
            std::make_unique<octree_node>(octree_node::find_overlapping_elements(list_p_elements, list_sub_boxes[index_sub_box]),
                                          list_sub_boxes[index_sub_box],
                                          child_is_root);
    }
}

std::vector<element *> octree_node::find_overlapping_elements(const std::vector<element *> &list_p_elements, const bbox &bounding_box) {
    constexpr uint         number_octant = 8;
    std::vector<element *> list_overlaping_elements;
    list_overlaping_elements.reserve(list_p_elements.size() / number_octant);
    std::copy_if(list_p_elements.begin(), list_p_elements.end(), std::back_inserter(list_overlaping_elements), [&](const auto &p_element) {
        return bounding_box.is_overlapping_tetra(*p_element);
    });
    list_overlaping_elements.shrink_to_fit();
    return list_overlaping_elements;
}

element *octree_node::find_element_at_location(const vector3 &position) const {
    if (m_is_leaf) {
        auto it_element = std::find_if(m_list_p_elements.begin(), m_list_p_elements.end(), [&](auto &&p_element) {
            return p_element->is_location_inside_element(position);
        });
        if (it_element != m_list_p_elements.end()) {
            return *it_element;
        }
        // LOG_ERROR << "WARNING : NO ELEMENT AT LOCATION OCTREE FAILED." << std::endl;
        return nullptr;

    }  // Descend down in the quadtree
    if (m_bottom_left_subbox_front.is_inside(position)) {
        return m_sub_nodes[0]->find_element_at_location(position);
    }
    if (m_bottom_right_subbox_front.is_inside(position)) {
        return m_sub_nodes[1]->find_element_at_location(position);
    }
    if (m_top_right_subbox_front.is_inside(position)) {
        return m_sub_nodes[2]->find_element_at_location(position);
    }
    if (m_top_left_subbox_front.is_inside(position)) {
        return m_sub_nodes[3]->find_element_at_location(position);
    }
    if (m_bottom_left_subbox_back.is_inside(position)) {
        return m_sub_nodes[4]->find_element_at_location(position);
    }
    if (m_bottom_right_subbox_back.is_inside(position)) {
        return m_sub_nodes[5]->find_element_at_location(position);
    }
    if (m_top_right_subbox_back.is_inside(position)) {
        return m_sub_nodes[6]->find_element_at_location(position);
    }
    if (m_top_left_subbox_back.is_inside(position)) {
        return m_sub_nodes[7]->find_element_at_location(position);
    }

    return nullptr;
}

std::vector<element *> octree_node::find_elements_overlapping_box(const bbox &box) const {
    std::set<element *> set_overlapping_elements;
    if (m_is_leaf) {
        std::copy_if(m_list_p_elements.begin(),
                     m_list_p_elements.end(),
                     std::inserter(set_overlapping_elements, set_overlapping_elements.end()),
                     [&](auto &&p_element) { return box.is_overlapping_tetra(*p_element); });
        return std::vector<element *>(set_overlapping_elements.begin(), set_overlapping_elements.end());
    }
    if (m_bottom_left_subbox_front.is_overlapping(box)) {
        std::vector<element *> list_overlapping_elements = m_sub_nodes[0]->find_elements_overlapping_box(box);
        set_overlapping_elements.insert(list_overlapping_elements.begin(), list_overlapping_elements.end());
    }
    if (m_bottom_right_subbox_front.is_overlapping(box)) {
        std::vector<element *> list_overlapping_elements = m_sub_nodes[1]->find_elements_overlapping_box(box);
        set_overlapping_elements.insert(list_overlapping_elements.begin(), list_overlapping_elements.end());
    }
    if (m_top_right_subbox_front.is_overlapping(box)) {
        std::vector<element *> list_overlapping_elements = m_sub_nodes[2]->find_elements_overlapping_box(box);
        set_overlapping_elements.insert(list_overlapping_elements.begin(), list_overlapping_elements.end());
    }
    if (m_top_left_subbox_front.is_overlapping(box)) {
        std::vector<element *> list_overlapping_elements = m_sub_nodes[3]->find_elements_overlapping_box(box);
        set_overlapping_elements.insert(list_overlapping_elements.begin(), list_overlapping_elements.end());
    }
    if (m_bottom_left_subbox_back.is_overlapping(box)) {
        std::vector<element *> list_overlapping_elements = m_sub_nodes[4]->find_elements_overlapping_box(box);
        set_overlapping_elements.insert(list_overlapping_elements.begin(), list_overlapping_elements.end());
    }
    if (m_bottom_right_subbox_back.is_overlapping(box)) {
        std::vector<element *> list_overlapping_elements = m_sub_nodes[5]->find_elements_overlapping_box(box);
        set_overlapping_elements.insert(list_overlapping_elements.begin(), list_overlapping_elements.end());
    }
    if (m_top_right_subbox_back.is_overlapping(box)) {
        std::vector<element *> list_overlapping_elements = m_sub_nodes[6]->find_elements_overlapping_box(box);
        set_overlapping_elements.insert(list_overlapping_elements.begin(), list_overlapping_elements.end());
    }
    if (m_top_left_subbox_back.is_overlapping(box)) {
        std::vector<element *> list_overlapping_elements = m_sub_nodes[7]->find_elements_overlapping_box(box);
        set_overlapping_elements.insert(list_overlapping_elements.begin(), list_overlapping_elements.end());
    }
    return std::vector<element *>(set_overlapping_elements.begin(), set_overlapping_elements.end());
}

}  // namespace mesh

}  // namespace uepm