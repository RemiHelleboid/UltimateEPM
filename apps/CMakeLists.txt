# Create the main executables
add_executable(EmpiricalPseudoPotentialMain EmpiricalPseudoPotentialMain.cpp)
target_link_libraries(EmpiricalPseudoPotentialMain PRIVATE libepp)
target_compile_features(EmpiricalPseudoPotentialMain PRIVATE cxx_std_20)

add_executable(BandsOnBZ BandsOnBZ.cpp)
target_link_libraries(BandsOnBZ PRIVATE libepp lib_bzmesh Eigen3::Eigen)
target_compile_features(BandsOnBZ PRIVATE cxx_std_20)

add_executable(DOS_MeshBZ DOS_MeshBZ.cpp)
target_link_libraries(DOS_MeshBZ PRIVATE libepp lib_bzmesh)
target_compile_features(DOS_MeshBZ PRIVATE cxx_std_20)

add_executable(fullstates.epm fullstates.cpp)
target_link_libraries(fullstates.epm PRIVATE libepp lib_bzmesh)
target_compile_features(fullstates.epm PRIVATE cxx_std_20)

add_executable(elph.epm elph.cpp)
target_link_libraries(elph.epm PRIVATE libepp lib_bzmesh fmt::fmt)
target_compile_features(elph.epm PRIVATE cxx_std_20)

add_executable(fbmc.epm fbmc.cpp)
target_link_libraries(fbmc.epm PRIVATE libepp lib_bzmesh lib_fbmc)
target_compile_features(fbmc.epm PRIVATE cxx_std_20)

add_executable(impact_io.epm impact_io.cpp)
target_link_libraries(impact_io.epm PRIVATE libepp lib_bzmesh)
target_compile_features(impact_io.epm PRIVATE cxx_std_20)

add_executable(bz_meshing.epm bz_meshing.cpp)
target_link_libraries(bz_meshing.epm PRIVATE libepp lib_bzmesh)
target_compile_features(bz_meshing.epm PRIVATE cxx_std_20)
target_include_directories(bz_meshing.epm PRIVATE ${GMSH_INC})
target_link_libraries(bz_meshing.epm PRIVATE ${GMSH_LIB})

if(USE_MPI_ACCELERATION)
  add_executable(mpiBandsOnBZ mpi_BandsOnBZ.cpp)
  target_link_libraries(mpiBandsOnBZ PRIVATE libepp Eigen3::Eigen)
  target_link_libraries(mpiBandsOnBZ PRIVATE MPI::MPI_CXX)
  target_compile_features(mpiBandsOnBZ PRIVATE cxx_std_20)

  add_executable(mpiDOS_MeshBZ mpi_DOS_MeshBZ.cpp)
  target_link_libraries(mpiDOS_MeshBZ PRIVATE libepp lib_bzmesh)
  target_link_libraries(mpiDOS_MeshBZ PRIVATE MPI::MPI_CXX)
  target_compile_features(mpiDOS_MeshBZ PRIVATE cxx_std_20)

  add_executable(epsilon.epm epsilon.cpp)
  target_link_libraries(epsilon.epm PRIVATE libepp lib_bzmesh)
  target_compile_features(epsilon.epm PRIVATE cxx_std_20)
  target_link_libraries(epsilon.epm PRIVATE MPI::MPI_CXX)
endif(USE_MPI_ACCELERATION)
