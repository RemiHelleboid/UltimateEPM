set(HEADER_FILES_MMMC
    mmmc_policy.hpp
    mmmc_particle_transfer.hpp
    mmmc_particle_pools.hpp
    self_consistent_device_mmmc_simulation_2d.hpp
    mmmc_device_config.hpp
    mmmc_device_runner.hpp
)

set(SOURCE_FILES_MMMC
    mmmc_policy.cpp
    mmmc_particle_transfer.cpp
    mmmc_particle_pools.cpp
    self_consistent_device_mmmc_simulation_2d.cpp
    mmmc_device_config.cpp
    mmmc_device_runner.cpp
)

add_library(lib_MMMC STATIC ${SOURCE_FILES_MMMC} ${HEADER_FILES_MMMC})
add_library(uepm::MMMC ALIAS lib_MMMC)

target_include_directories(lib_MMMC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(lib_MMMC PUBLIC uepm::PBMC uepm::ADMC uepm::mesh uepm::common)
target_compile_features(lib_MMMC PUBLIC cxx_std_20)
