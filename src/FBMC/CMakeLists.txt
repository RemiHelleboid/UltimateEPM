set(HEADER_FILES_LIBFBMC single_part_fbmc.hpp particle.hpp

  )


set(SOURCE_FILES_LIBFBMC
  single_part_fbmc.cpp particle.cpp
)

add_library(lib_fbmc STATIC ${SOURCE_FILES_LIBFBMC} ${HEADER_FILES_LIBFBMC})

target_include_directories(lib_fbmc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link with reuired libraries
target_link_libraries(lib_fbmc PUBLIC libepp lib_bzmesh)
target_link_libraries(lib_fbmc PUBLIC Eigen3::Eigen rapidcsv fmt::fmt)
target_include_directories(lib_fbmc PUBLIC ${GMSH_INC})
target_link_libraries(lib_fbmc PUBLIC ${GMSH_LIB})

if(ENABLE_OPEN_MP)
  target_link_libraries(lib_fbmc PUBLIC OpenMP::OpenMP_CXX)
endif(ENABLE_OPEN_MP)

# All users of this library will need at least C++20
target_compile_features(lib_fbmc PRIVATE cxx_std_20)