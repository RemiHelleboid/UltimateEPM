set(HEADER_FILES_LIBEPP
  SymmetryPoints.h
  bz_meshfile.hpp
  NonLocalParameters.hpp
  NonLocalFunctional.hpp
  SpinOrbitParameters.hpp
  SpinOrbitFunctional.hpp
  Options.h
  Material.h
  bessel_func.hpp
  Hamiltonian.h
  Vector3D.h
  Vector3D.inl
  BandStructure.h
  DielectricFunction.hpp
  Pseudopotential.h)

set(SOURCE_FILES_LIBEPP
  Material.cpp
  SymmetryPoints.cpp
  Options.cpp
  Hamiltonian.cpp
  NonLocalFunctional.cpp
  BandStructure.cpp
  Pseudopotential.cpp
  SpinOrbitFunctional.cpp
  DielectricFunction.cpp
  bz_meshfile.cpp)

add_library(libepp STATIC ${SOURCE_FILES_LIBEPP} ${HEADER_FILES_LIBEPP})

target_include_directories(libepp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link with reuired libraries
target_link_libraries(libepp PUBLIC Eigen3::Eigen yaml-cpp rapidcsv fmt::fmt)
target_include_directories(libepp PUBLIC ${GMSH_INC})
target_link_libraries(libepp PUBLIC ${GMSH_LIB})

if(USE_OPENMP_ACCELERATION)
  target_link_libraries(libepp PUBLIC OpenMP::OpenMP_CXX)
endif(USE_OPENMP_ACCELERATION)

# All users of this library will need at least C++20
target_compile_features(libepp PRIVATE cxx_std_20)
