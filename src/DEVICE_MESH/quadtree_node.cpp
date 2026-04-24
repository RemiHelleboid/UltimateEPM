/**
 * @file quadtree_node.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-29
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "quadtree_node.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

#include "bbox.hpp"
#include "element.hpp"
#include "element2d.hpp"
#include "tree_node.hpp"

namespace uepm {

namespace mesh {

constexpr int    Max_Element_Per_Quadtree_Node = 8;
constexpr double Min_Box_Size                  = 1e-3;

static const std::string out_file = "QuadTreeSave.txt";

quadtree_node::quadtree_node(const std::vector<element *> &list_p_elements, const bbox &bounding_box, int generation)
    : m_node_box(bounding_box),
      m_generation_depth(generation) {
    if (list_p_elements.size() <= Max_Element_Per_Quadtree_Node || m_node_box.get_diagonal_size() <= Min_Box_Size) {
        // std::ofstream my_file;
        // my_file.open(out_file, std::ios::app);
        // my_file << bounding_box << "," << list_p_elements.size() << std::endl;
        // my_file.close();
        this->m_list_p_elements = list_p_elements;
        m_is_leaf               = true;
        return;
    }
    const std::vector<bbox> list_sub_boxes = bounding_box.split_2d_box_in_quadrants();
    m_bottom_left_subbox                   = list_sub_boxes[0];
    m_bottom_right_subbox                  = list_sub_boxes[1];
    m_top_right_subbox                     = list_sub_boxes[2];
    m_top_left_subbox                      = list_sub_boxes[3];

    constexpr uint number_quadrant = 4;
    m_sub_nodes.resize(number_quadrant);
    std::array<std::vector<element *>, number_quadrant> list_overlap_element_quadrants;
    int                                                 new_generation_depth = m_generation_depth + 1;
    // #pragma omp parallel for
    for (uint index_sub_box = 0; index_sub_box < number_quadrant; ++index_sub_box) {
        m_sub_nodes[index_sub_box] =
            std::make_unique<quadtree_node>(quadtree_node::find_overlapping_elements(list_p_elements, list_sub_boxes[index_sub_box]),
                                            list_sub_boxes[index_sub_box],
                                            new_generation_depth);
    }
}

std::vector<element *> quadtree_node::find_overlapping_elements(const std::vector<element *> &list_p_elements, const bbox &bounding_box) {
    std::vector<element *> list_overlaping_elements;
    list_overlaping_elements.reserve(list_p_elements.size() / 3);
    std::copy_if(list_p_elements.begin(), list_p_elements.end(), std::back_inserter(list_overlaping_elements), [&](const auto &p_element) {
        return bounding_box.is_overlapping_triangle(*p_element);
    });
    list_overlaping_elements.shrink_to_fit();
    return list_overlaping_elements;
}

element *quadtree_node::find_element_at_location(const vector3 &position) const {
    if (m_is_leaf) {
        auto it_element = std::find_if(m_list_p_elements.begin(), m_list_p_elements.end(), [&](auto &&p_element) {
            return p_element->is_location_inside_element(position);
        });
        if (it_element != m_list_p_elements.end()) {
            return *it_element;
        }
        // LOG_ERROR << "MAJOR FAILURE : NO ELEMENT AT LOCATION QUADTREE FAILED. " << position;
        return nullptr;
    }
    // Descend down in the quadtree
    if (m_bottom_left_subbox.is_inside(position)) {
        return m_sub_nodes[0]->find_element_at_location(position);
    }
    if (m_bottom_right_subbox.is_inside(position)) {
        return m_sub_nodes[1]->find_element_at_location(position);
    }
    if (m_top_right_subbox.is_inside(position)) {
        return m_sub_nodes[2]->find_element_at_location(position);
    }
    if (m_top_left_subbox.is_inside(position)) {
        return m_sub_nodes[3]->find_element_at_location(position);
    }
    return nullptr;
}

std::vector<element *> quadtree_node::find_elements_overlapping_box(const bbox &box) const {
    std::set<element *> set_overlaping_elements;
    if (m_is_leaf) {
        std::copy_if(m_list_p_elements.begin(),
                     m_list_p_elements.end(),
                     std::inserter(set_overlaping_elements, set_overlaping_elements.begin()),
                     [&](const auto &p_element) { return box.is_overlapping_triangle(*p_element); });
        return std::vector<element *>(set_overlaping_elements.begin(), set_overlaping_elements.end());
    }
    // Descend down in the quadtree
    if (m_bottom_left_subbox.is_overlapping(box)) {
        std::vector<element *> list_overlaping_elements = m_sub_nodes[0]->find_elements_overlapping_box(box);
        set_overlaping_elements.insert(list_overlaping_elements.begin(), list_overlaping_elements.end());
    }
    if (m_bottom_right_subbox.is_overlapping(box)) {
        std::vector<element *> list_overlaping_elements = m_sub_nodes[1]->find_elements_overlapping_box(box);
        set_overlaping_elements.insert(list_overlaping_elements.begin(), list_overlaping_elements.end());
    }
    if (m_top_right_subbox.is_overlapping(box)) {
        std::vector<element *> list_overlaping_elements = m_sub_nodes[2]->find_elements_overlapping_box(box);
        set_overlaping_elements.insert(list_overlaping_elements.begin(), list_overlaping_elements.end());
    }
    if (m_top_left_subbox.is_overlapping(box)) {
        std::vector<element *> list_overlaping_elements = m_sub_nodes[3]->find_elements_overlapping_box(box);
        set_overlaping_elements.insert(list_overlaping_elements.begin(), list_overlaping_elements.end());
    }
    return std::vector<element *>(set_overlaping_elements.begin(), set_overlaping_elements.end());
}

}  // namespace mesh

}  // namespace uepm