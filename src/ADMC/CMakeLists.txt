set(HEADER_FILES_ADMC
    admc_carrier.hpp
    admc_mobility.hpp
    admc_transport.hpp
    bulk_admc_simulation.hpp
)

set(SOURCE_FILES_ADMC
    admc_mobility.cpp
    admc_transport.cpp
    bulk_admc_simulation.cpp
)

add_library(lib_ADMC STATIC ${SOURCE_FILES_ADMC} ${HEADER_FILES_ADMC})
add_library(uepm::ADMC ALIAS lib_ADMC)

target_include_directories(lib_ADMC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(lib_ADMC PUBLIC uepm::common uepm::physics yaml-cpp::yaml-cpp)
target_compile_features(lib_ADMC PUBLIC cxx_std_20)
