# Create the main executables
add_executable(EmpiricalPseudoPotentialMain EmpiricalPseudoPotentialMain.cpp)
target_link_libraries(EmpiricalPseudoPotentialMain PUBLIC libepp)
target_compile_features(EmpiricalPseudoPotentialMain PRIVATE cxx_std_20)

add_executable(BandsOnBZ BandsOnBZ.cpp)
target_link_libraries(BandsOnBZ PUBLIC libepp Eigen3::Eigen)
target_compile_features(BandsOnBZ PRIVATE cxx_std_20)

add_executable(DOS_MeshBZ DOS_MeshBZ.cpp)
target_link_libraries(DOS_MeshBZ PUBLIC libepp lib_bzmesh)
target_compile_features(DOS_MeshBZ PRIVATE cxx_std_20)

add_executable(epsilon.epm epsilon.cpp)
target_link_libraries(epsilon.epm PUBLIC libepp lib_bzmesh)
target_compile_features(epsilon.epm PRIVATE cxx_std_20)
target_link_libraries(epsilon.epm PUBLIC MPI::MPI_CXX)

add_executable(fullstates.epm fullstates.cpp)
target_link_libraries(fullstates.epm PUBLIC libepp lib_bzmesh)
target_compile_features(fullstates.epm PRIVATE cxx_std_20)
target_link_libraries(fullstates.epm PUBLIC MPI::MPI_CXX)

add_executable(elph.epm elph.cpp)
target_link_libraries(elph.epm PUBLIC libepp lib_bzmesh)
target_compile_features(elph.epm PRIVATE cxx_std_20)
target_link_libraries(elph.epm PUBLIC MPI::MPI_CXX)

add_executable(fbmc.epm fbmc.cpp)
target_link_libraries(fbmc.epm PUBLIC libepp lib_bzmesh lib_fbmc)
target_compile_features(fbmc.epm PRIVATE cxx_std_20)
target_link_libraries(fbmc.epm PUBLIC MPI::MPI_CXX)

add_executable(impact_io.epm impact_io.cpp)
target_link_libraries(impact_io.epm PUBLIC libepp lib_bzmesh)
target_compile_features(impact_io.epm PRIVATE cxx_std_20)
target_link_libraries(impact_io.epm PUBLIC MPI::MPI_CXX)

add_executable(bz_meshing.epm bz_meshing.cpp)
target_link_libraries(bz_meshing.epm PUBLIC libepp lib_bzmesh)
target_compile_features(bz_meshing.epm PRIVATE cxx_std_20)
target_link_libraries(bz_meshing.epm PUBLIC MPI::MPI_CXX)
target_include_directories(bz_meshing.epm  PUBLIC ${GMSH_INC})
target_link_libraries(bz_meshing.epm  PUBLIC ${GMSH_LIB})


if(USE_MPI_ACCELERATION)
  add_executable(mpiBandsOnBZ mpi_BandsOnBZ.cpp)
  target_link_libraries(mpiBandsOnBZ PUBLIC libepp Eigen3::Eigen)
  target_link_libraries(mpiBandsOnBZ PUBLIC MPI::MPI_CXX)
  target_compile_features(mpiBandsOnBZ PRIVATE cxx_std_20)

  add_executable(mpiDOS_MeshBZ mpi_DOS_MeshBZ.cpp)
  target_link_libraries(mpiDOS_MeshBZ PUBLIC libepp lib_bzmesh)
  target_link_libraries(mpiDOS_MeshBZ PUBLIC MPI::MPI_CXX)
  target_compile_features(mpiDOS_MeshBZ PRIVATE cxx_std_20)
endif(USE_MPI_ACCELERATION)

