/**
 * @file tree_node.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2021-10-27
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <memory>
#include <vector>

#include "bbox.hpp"

namespace uepm {

namespace mesh {

class tree_node {
 protected:
    bool                   m_is_leaf = false;
    std::vector<element *> m_list_p_elements;

    virtual std::vector<element *> find_overlapping_elements(const std::vector<element *> &list_p_elements, const bbox &bounding_box) = 0;

 public:
    tree_node(){};
    virtual ~tree_node(){};
    virtual element *              find_element_at_location(const vector3 &location) const = 0;
    virtual std::vector<element *> find_elements_overlapping_box(const bbox &box) const    = 0;
    //  virtual std::vector<element *> find_elements_overlapping_sphere(const vector3 sphere_center, double radius) const = 0;
};

}  //  namespace mesh

}  // namespace uepm