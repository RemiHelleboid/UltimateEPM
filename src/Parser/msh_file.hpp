/**
 * @file msh_file.hpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief header file for file::msh_file class.
 * @version 0.1
 * @date 2021-07-10
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <cassert>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "dataset.hpp"
#include "file.hpp"

static const std::map<std::string, uepm::mesh::DataLocationType> msh_to_armin_data_location_type{{"NodeData", uepm::mesh::DataLocationType::vertex},
                                                                                                   {"ElementData", uepm::mesh::DataLocationType::cell}};

#ifdef ST_VERSION
#include "STF_writer.hpp"
#endif

namespace uepm {

namespace file {

/**
 * @class msh_file
 * @brief Class representing .msh file format (gmsh file format).
 *
 */
class msh_file : public file {
 protected:
 public:
    ~msh_file() = default;

    msh_file() = default;
    explicit msh_file(const std::string &filepath);

    void read_mesh() override;
    void read_states(const std::vector<std::string> &list_dataset_to_import = {}) override;
    void open_file();

    void export_as_msh(const std::string &      filename,
                       std::vector<std::string> datasets_to_export = {},
                       bool                     export_all_dataset = 0) override {
        uepm::mesh::writer::export_as_msh(this->m_Mesh, filename, datasets_to_export, export_all_dataset);
    };


};

}  //  namespace file

}  // namespace uepm
