
set(HEADER_FILES_LIBFILE "file.hpp"
    "msh_file.hpp"
    "vtkWriter.hpp"
)

set(SOURCE_FILES_LIBFILE "file.cpp"
    "msh_file.cpp"
    "vtkWriter.cpp"
)

add_library(libfile ${HEADER_FILES_LIBFILE} ${SOURCE_FILES_LIBFILE})


target_include_directories(libfile PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(libfile PRIVATE ${GMSH_LIB})
target_include_directories(libfile PRIVATE ${GMSH_INC})

if (ENABLE_ST_VERSION)
    message(STATUS "Adding STF file support")
    target_include_directories(libfile PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(libfile PRIVATE ${HDF5_HL_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
endif()

target_link_libraries(libfile PRIVATE libmesh)

