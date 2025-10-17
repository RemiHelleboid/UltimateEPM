/**
 * @file gmsh_guard.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2025-10-17
 * 
 * 
 */

#include "gmsh.h"

namespace uepm::mesh_bz {
struct GmshSession {
    explicit GmshSession(int verbosity = 0) {
        gmsh::initialize();
        gmsh::option::setNumber("General.Verbosity", verbosity);
    }
    ~GmshSession() noexcept {
        try {
            gmsh::finalize();
        } catch (...) {
        }
    }
};
}  // namespace uepm::mesh_bz
