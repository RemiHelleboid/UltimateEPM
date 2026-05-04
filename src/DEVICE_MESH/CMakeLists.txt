set(HEADER_FILES_LIBMESH
    "geometry_entity.hpp"
    "vertex.hpp"
    "vector.hpp"
    "bbox.hpp"
    "element.hpp"
    "element1d.hpp"
    "element2d.hpp"
    "element3d.hpp"
    "region.hpp"
    "mesh.hpp"
    "dataset.hpp"
    "function.hpp"
    "writer.hpp"
    "tree_node.hpp"
    "quadtree_node.hpp"
    "octree_node.hpp"
    "utils_mesh.hpp"
    "quadrangular_interpolators.hpp"
    "grid_data.hpp"
    "device.hpp"
    "contact.hpp"
    )

set(SOURCE_FILES_LIBMESH
    "geometry_entity.cpp"
    "vertex.cpp"
    "bbox.cpp"
    "element.cpp"
    "element1d.cpp"
    "element2d.cpp"
    "element3d.cpp"
    "region.cpp"
    "mesh.cpp"
    "dataset.cpp"
    "writer.cpp"
    "quadtree_node.cpp"
    "device.cpp"
    "octree_node.cpp")

# Add the STF library for reading STF files
add_library(libmesh STATIC ${SOURCE_FILES_LIBMESH} ${HEADER_FILES_LIBMESH})

target_include_directories(libmesh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(libmesh PUBLIC ${GMSH_INC})

# Add the STF library for reading STF files
target_include_directories(libmesh PUBLIC ${GMSH_INC})
target_link_libraries(libmesh PUBLIC ${GMSH_LIB})
target_link_libraries(libmesh PUBLIC uepm_common plog rapidcsv)
target_link_libraries(libmesh PUBLIC uepm::common)

if(USE_OPENMP_ACCELERATION)
  target_link_libraries(libmesh PUBLIC OpenMP::OpenMP_CXX)
endif(USE_OPENMP_ACCELERATION)

# All users of this library will need at least C++20
target_compile_features(libmesh PRIVATE cxx_std_20)
