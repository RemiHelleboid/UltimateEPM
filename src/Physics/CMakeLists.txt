set(HEADER_FILES_PHYSICS "materials.hpp"
    "physical_constants.hpp"
    "physical_statistics.hpp"
    # "impact_ionization.hpp"
    "physical_constants.hpp"
    "physical_statistics.hpp"
)

set(SOURCE_FILES_PHYSICS materials.cpp
)

# Add the STF library for reading STF files
add_library(libphysics STATIC ${SOURCE_FILES_PHYSICS} ${HEADER_FILES_PHYSICS})

target_include_directories(libphysics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(libphysics PUBLIC yaml-cpp)
target_link_libraries(libphysics PUBLIC uepm::common)

