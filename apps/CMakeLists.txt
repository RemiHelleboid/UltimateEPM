# Create the main executables
add_executable(EmpiricalPseudoPotentialMain EmpiricalPseudoPotentialMain.cpp)
target_link_libraries(EmpiricalPseudoPotentialMain PUBLIC libepp pthread)
target_compile_features(EmpiricalPseudoPotentialMain PRIVATE cxx_std_20)

add_executable(BandsOnBZ BandsOnBZ.cpp)
target_link_libraries(BandsOnBZ PUBLIC libepp pthread Eigen3::Eigen )
target_compile_features(BandsOnBZ PRIVATE cxx_std_20)

add_executable(DOS_MeshBZ DOS_MeshBZ.cpp)
target_link_libraries(DOS_MeshBZ PUBLIC libepp lib_bzmesh pthread)
target_compile_features(DOS_MeshBZ PRIVATE cxx_std_20)

if(USE_MPI_ACCELERATION)
  add_executable(mpiBandsOnBZ mpi_BandsOnBZ.cpp)
  target_link_libraries(mpiBandsOnBZ PUBLIC libepp pthread Eigen3::Eigen)
  target_link_libraries(mpiBandsOnBZ PUBLIC MPI::MPI_CXX)
  target_compile_features(mpiBandsOnBZ PRIVATE cxx_std_20)
endif(USE_MPI_ACCELERATION)
