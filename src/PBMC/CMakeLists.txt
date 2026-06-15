set(HEADER_FILES_LIBPBMC
    valley_model.hpp
    pbmc_particle.hpp
    bulk_pbmc_simulation.hpp
    scattering_events.hpp
    intervalley_phonon.hpp
    scattering_channels.hpp
    device_pbmc_simulation.hpp
    pbmc_transport_kernel.hpp
    pbmc_material_model.hpp
    pbmc_scattering_model.hpp
    pbmc_device_history.hpp
    pbmc_self_consistent_device_simulation_base.hpp
    pbmc_self_consistent_device_simulation_3d.hpp
    pbmc_self_consistent_device_simulation_2d.hpp
    pbmc_quench_circuit.hpp
    pbmc_avalanche_detector.hpp
    pbmc_quench_detector.hpp
    pbmc_device_setup.hpp
    pbmc_device_config.hpp
    pbmc_device_runner.hpp
    pbmc_run_manifest.hpp
  )


set(SOURCE_FILES_LIBPBMC
    pbmc_particle.cpp
    bulk_pbmc_simulation.cpp
    valley_model.cpp
    intervalley_phonon.cpp
    device_pbmc_simulation.cpp
    pbmc_transport_kernel.cpp
    pbmc_material_model.cpp
    pbmc_scattering_model.cpp
    pbmc_self_consistent_device_simulation_base.cpp
    pbmc_self_consistent_device_simulation_3d.cpp
    pbmc_self_consistent_device_simulation_2d.cpp
    pbmc_quench_circuit.cpp
    pbmc_device_setup.cpp
    pbmc_device_config.cpp
    pbmc_device_runner.cpp
    pbmc_run_manifest.cpp
)

add_library(lib_PBMC STATIC ${SOURCE_FILES_LIBPBMC} ${HEADER_FILES_LIBPBMC})
add_library(uepm::PBMC ALIAS lib_PBMC)

target_include_directories(lib_PBMC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link with reuired libraries
target_link_libraries(lib_PBMC PUBLIC libepp lib_bzmesh libfile libmesh libphysics uepm::poisson)
target_link_libraries(lib_PBMC PUBLIC Eigen3::Eigen rapidcsv fmt::fmt yaml-cpp::yaml-cpp)

target_link_libraries(lib_PBMC PUBLIC OpenMP::OpenMP_CXX)

# All users of this library will need at least C++20
target_compile_features(lib_PBMC PUBLIC cxx_std_20)
target_compile_definitions(
  lib_PBMC
  PUBLIC
    UEPM_PROJECT_VERSION="${PROJECT_VERSION}"
    UEPM_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
    UEPM_COMPILER="${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}"
)
