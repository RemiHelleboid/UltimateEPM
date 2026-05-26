set(HEADER_FILES_LIBAMC
    valley_model.hpp
    particle_amc.hpp
    bulk_amc_simulation.hpp
    scattering_events.hpp
    intervalley_phonon.hpp
    scattering_channels.hpp
    device_amc_simulation.hpp
    amc_transport_kernel.hpp
    amc_material_model.hpp
    amc_scattering_model.hpp
    amc_device_history.hpp
    amc_self_consistent_device_simulation_3d.hpp
  )


set(SOURCE_FILES_LIBAMC
    particle_amc.cpp
    bulk_amc_simulation.cpp
    valley_model.cpp
    intervalley_phonon.cpp
    device_amc_simulation.cpp
    amc_transport_kernel.cpp
    amc_material_model.cpp
    amc_scattering_model.cpp
    amc_self_consistent_device_simulation_3d.cpp
)

add_library(lib_amc STATIC ${SOURCE_FILES_LIBAMC} ${HEADER_FILES_LIBAMC})
add_library(uepm::amc ALIAS lib_amc)

target_include_directories(lib_amc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link with reuired libraries
target_link_libraries(lib_amc PUBLIC libepp lib_bzmesh libfile libmesh libphysics uepm::poisson)
target_link_libraries(lib_amc PUBLIC Eigen3::Eigen rapidcsv fmt::fmt)
target_include_directories(lib_amc PUBLIC ${GMSH_INC})
target_link_libraries(lib_amc PUBLIC ${GMSH_LIB})


  target_link_libraries(lib_amc PUBLIC OpenMP::OpenMP_CXX)

# All users of this library will need at least C++20
target_compile_features(lib_amc PUBLIC cxx_std_20)