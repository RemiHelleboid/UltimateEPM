set(HEADER_FILES_PHYSICS "materials.hpp"
)

set(SOURCE_FILES_PHYSICS materials.cpp
)

# Add the STF library for reading STF files
add_library(libphysics STATIC ${SOURCE_FILES_PHYSICS} ${HEADER_FILES_PHYSICS})
add_library(uepm::physics ALIAS libphysics)

target_include_directories(libphysics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(libphysics PUBLIC yaml-cpp)
target_link_libraries(libphysics PUBLIC uepm::common)
target_compile_definitions(libphysics PRIVATE UEPM_MATERIAL_DATA_DIR="${UEPM_MATERIAL_DATA_DIR}")
