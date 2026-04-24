/**
 * @file writer.hpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief
 * @version 0.1
 * @date 2021-07-23
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <string>
#include <vector>

#include "dataset.hpp"
#include "mesh.hpp"
#include "region.hpp"

namespace uepm {

namespace mesh {

static const std::vector<std::string> empty_vector_string{};

class writer {
 private:
 public:
    writer()                              = delete;
    writer(const writer &)                = delete;
    writer     &operator=(const writer &) = delete;
    static void export_as_msh(const mesh                     &myMesh,
                              const std::string              &filepath,
                              const std::vector<std::string> &datasets_to_export = {},
                              bool                            export_all_dataset = false);
};

}  //  namespace mesh

}  // namespace uepm