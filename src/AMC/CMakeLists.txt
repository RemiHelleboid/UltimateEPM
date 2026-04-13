set(HEADER_FILES_LIBAMC
    valley_model.hpp
    particle_amc.hpp
    bulk_amc_simulation.hpp
    scattering_events.hpp
  )


set(SOURCE_FILES_LIBAMC
    particle_amc.cpp
    bulk_amc_simulation.cpp
    valley_model.cpp
)

add_library(lib_amc STATIC ${SOURCE_FILES_LIBAMC} ${HEADER_FILES_LIBAMC})

target_include_directories(lib_amc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link with reuired libraries
target_link_libraries(lib_amc PUBLIC libepp lib_bzmesh)
target_link_libraries(lib_amc PUBLIC Eigen3::Eigen rapidcsv fmt::fmt)
target_include_directories(lib_amc PUBLIC ${GMSH_INC})
target_link_libraries(lib_amc PUBLIC ${GMSH_LIB})


  target_link_libraries(lib_amc PUBLIC OpenMP::OpenMP_CXX)

# All users of this library will need at least C++20
target_compile_features(lib_amc PUBLIC cxx_std_20)