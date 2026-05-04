/**
 * @file file.hpp
 * @author Rémi Helleboid (remi@helleboid.net)
 * @brief header file for base class file::file.
 * @version 0.1
 * @date 2021-07-16
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <plog/Log.h>

#include <cassert>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "mesh.hpp"
#include "writer.hpp"

namespace uepm {
namespace file {

/**
 * @class file
 * @brief file class : abstract class for file mesh I/O, which is derived for each specific file format.
 *
 */
class file {
 protected:
    std::string      m_file_path;
    uepm::mesh::mesh m_Mesh;

 public:
    file() = default;
    explicit file(const std::string& filename);

    void set_file_path(const std::string& filename) { m_file_path = filename; }
    virtual ~file();

    uepm::mesh::mesh*  get_p_mesh() { return &m_Mesh; }
    uepm::mesh::mesh&& get_mesh() { return std::move(m_Mesh); }

    void apply_mundfabisation(const std::string& mundfab_data_name, const double value_silicon, const std::string& silicon_region_name) {
        m_Mesh.mundfabisation(mundfab_data_name, value_silicon, silicon_region_name);
    }
    virtual void read_mesh()                                                              = 0;
    virtual void read_states(const std::vector<std::string>& list_dataset_to_import = {}) = 0;

    void         export_vertices_to_csv(std::string filename) const { m_Mesh.export_vertices_to_csv(filename); }
    virtual void export_as_msh(const std::string&       filename,
                               std::vector<std::string> datasets_to_export = {},
                               bool                     export_all_dataset = 0) = 0;

#ifdef ST_VERSION
    virtual void export_as_STF_from_template(const std::string& filename, const std::string& template_file) = 0;
#endif
};

}  //  namespace file.

}  // namespace uepm