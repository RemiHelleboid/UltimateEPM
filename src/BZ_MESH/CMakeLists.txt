set(HEADER_FILES_LIBBZMESH
  mesh_tetra.hpp
  mesh_vertex.hpp
  bz_mesh.hpp
  vector.hpp
  iso_triangle.hpp
  bz_states.hpp
  electron_phonon.hpp
  impact_ionization.hpp
  dielectric_mesh.hpp
  )


set(SOURCE_FILES_LIBBZMESH
  mesh_tetra.cpp
  bz_mesh.cpp
  bz_states.cpp
  electron_phonon.cpp
  impact_ionization.cpp
  dielectric_mesh.cpp
)

add_library(lib_bzmesh STATIC ${SOURCE_FILES_LIBBZMESH} ${HEADER_FILES_LIBBZMESH})

target_include_directories(lib_bzmesh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link with reuired libraries
target_link_libraries(lib_bzmesh PUBLIC libepp)
target_link_libraries(lib_bzmesh PUBLIC Eigen3::Eigen rapidcsv fmt::fmt)
target_include_directories(lib_bzmesh PUBLIC ${GMSH_INC})
target_link_libraries(lib_bzmesh PUBLIC ${GMSH_LIB})

if(USE_OPENMP_ACCELERATION)
  target_link_libraries(lib_bzmesh PUBLIC OpenMP::OpenMP_CXX)
endif(USE_OPENMP_ACCELERATION)

# All users of this library will need at least C++20
target_compile_features(lib_bzmesh PRIVATE cxx_std_20)