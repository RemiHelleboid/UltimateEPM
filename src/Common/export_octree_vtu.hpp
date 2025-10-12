/**
 * @file export_octree_vtu.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2025-10-12
 * 
 * 
 */

#pragma once
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

// Adjust these includes/namespaces to your project
#include "bbox_mesh.hpp"
#include "octree_bz.hpp"

namespace uepm::mesh_bz {

// Adapt these accessors to your bbox_mesh API
struct Vec3 {
    double x, y, z;
};

inline Vec3 bb_min(const bbox_mesh& b) {
    const auto& v = b.min_corner();  // <-- change to your getter
    return {v.x(), v.y(), v.z()};
}
inline Vec3 bb_max(const bbox_mesh& b) {
    const auto& v = b.max_corner();  // <-- change to your getter
    return {v.x(), v.y(), v.z()};
}

struct VTUBuffer {
    std::vector<double>        points;   // xyz triples
    std::vector<std::uint32_t> conn;     // 8 indices per hex
    std::vector<std::uint32_t> offsets;  // 8,16,24,...
    std::vector<std::uint8_t>  types;    // all 12 (VTK_HEXAHEDRON)

    std::vector<int> depth;
    std::vector<int> is_leaf;
    std::vector<int> leaf_count;
};

inline void append_box(VTUBuffer& buf, const bbox_mesh& b, int d, bool leaf, int nleaf) {
    const auto mn = bb_min(b);
    const auto mx = bb_max(b);

    const std::uint32_t base = static_cast<std::uint32_t>(buf.points.size() / 3);

    // Push 8 points in VTK hex order (0..7)
    const double xs[2] = {mn.x, mx.x};
    const double ys[2] = {mn.y, mx.y};
    const double zs[2] = {mn.z, mx.z};

    auto add = [&](int xi, int yi, int zi) {
        buf.points.push_back(xs[xi]);
        buf.points.push_back(ys[yi]);
        buf.points.push_back(zs[zi]);
    };
    // zmin layer
    add(0, 0, 0);
    add(1, 0, 0);
    add(1, 1, 0);
    add(0, 1, 0);
    // zmax layer
    add(0, 0, 1);
    add(1, 0, 1);
    add(1, 1, 1);
    add(0, 1, 1);

    // connectivity: 8 consecutive indices
    for (std::uint32_t i = 0; i < 8; ++i) {
        buf.conn.push_back(base + i);
    }
    // offsets: cumulative node count
    buf.offsets.push_back(static_cast<std::uint32_t>(buf.conn.size()));
    // type: 12
    buf.types.push_back(12);

    buf.depth.push_back(d);
    buf.is_leaf.push_back(leaf ? 1 : 0);
    buf.leaf_count.push_back(leaf ? nleaf : 0);
}

inline void collect_octree_nodes(const Octree_mesh& node, VTUBuffer& buf, int depth, bool leaves_only) {
    if (node.is_leaf()) {
        append_box(buf, node.get_bounding_box(), depth, true, static_cast<int>(node.get_list_tetras().size()));
        return;
    }
    if (!leaves_only) {
        append_box(buf, node.get_bounding_box(), depth, false, 0);
    }
    for (const auto& child : node.get_list_sub_nodes()) {
        if (child) {
            collect_octree_nodes(*child, buf, depth + 1, leaves_only);
        }
    }
}

inline bool write_octree_as_vtu(const Octree_mesh& root, const std::string& filename, bool leaves_only = true) {
    VTUBuffer buf;
    buf.points.reserve(8 * 3 * 1024);
    buf.conn.reserve(8 * 1024);
    buf.offsets.reserve(1024);
    buf.types.reserve(1024);

    collect_octree_nodes(root, buf, /*depth=*/0, leaves_only);

    std::ofstream out(filename);
    if (!out) {
        return false;
    }

    // Minimal ASCII VTU
    out <<
        R"(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints=")"
        << (buf.points.size() / 3) << R"(" NumberOfCells=")" << buf.offsets.size() << R"(">
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
)";

    // Points
    for (size_t i = 0; i < buf.points.size(); i += 3) {
        out << buf.points[i] << " " << buf.points[i + 1] << " " << buf.points[i + 2] << "\n";
    }
    out <<
        R"(        </DataArray>
      </Points>
      <Cells>
        <DataArray type="UInt32" Name="connectivity" format="ascii">
)";
    // Connectivity
    for (size_t i = 0; i < buf.conn.size(); i += 8) {
        out << buf.conn[i] << " " << buf.conn[i + 1] << " " << buf.conn[i + 2] << " " << buf.conn[i + 3] << " " << buf.conn[i + 4] << " "
            << buf.conn[i + 5] << " " << buf.conn[i + 6] << " " << buf.conn[i + 7] << "\n";
    }
    out <<
        R"(        </DataArray>
        <DataArray type="UInt32" Name="offsets" format="ascii">
)";
    for (auto v : buf.offsets) {
        out << v << "\n";
    }
    out <<
        R"(        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">
)";
    for (auto v : buf.types) {
        out << static_cast<int>(v) << "\n";
    }
    out <<
        R"(        </DataArray>
      </Cells>
      <CellData Scalars="depth">
        <DataArray type="Int32" Name="depth" format="ascii">
)";
    for (auto v : buf.depth) {
        out << v << "\n";
    }
    out <<
        R"(        </DataArray>
        <DataArray type="Int32" Name="is_leaf" format="ascii">
)";
    for (auto v : buf.is_leaf) {
        out << v << "\n";
    }
    out <<
        R"(        </DataArray>
        <DataArray type="Int32" Name="leaf_count" format="ascii">
)";
    for (auto v : buf.leaf_count) {
        out << v << "\n";
    }
    out <<
        R"(        </DataArray>
      </CellData>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)";
    return true;
}

}  // namespace uepm::mesh_bz
