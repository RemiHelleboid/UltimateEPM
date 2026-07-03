set(HEADER_FILES_LIBBZMESH
    gmsh_guard.hpp
    band_catalog.hpp
    bz_domain.hpp
    bz_dos.hpp
    fermi_level.hpp
    mesh_tetra.hpp
    mesh_vertex.hpp
    bz_mesh.hpp
    vector_bz.hpp
    iso_triangle.hpp
    bz_states.hpp
    electron_phonon.hpp
    impact_ionization.hpp
    impact_ionization_state.hpp
    dielectric_mesh.hpp
    screened_coulomb_kernel.hpp
    bbox_mesh.hpp
    octree_bz.hpp
    reciprocal_space.hpp
    tetra_energy_index.hpp
    elph_common.hpp
    elph_deformation_potential.hpp
    elph_dispersion.hpp
    elph_overlap_integral.hpp)

set(SOURCE_FILES_LIBBZMESH
    mesh_tetra.cpp
    band_catalog.cpp
    bz_domain.cpp
    bz_dos.cpp
    bz_mesh_export.cpp
    bz_mesh.cpp
    bz_states.cpp
    electron_phonon.cpp
    impact_ionization.cpp
    impact_ionization_state.cpp
    fermi_level.cpp
    dielectric_mesh.cpp
    screened_coulomb_kernel.cpp
    octree_bz.cpp
    reciprocal_space.cpp
    tetra_energy_index.cpp)

add_library(lib_bzmesh STATIC ${SOURCE_FILES_LIBBZMESH}
                              ${HEADER_FILES_LIBBZMESH})
add_library(uepm::bzmesh ALIAS lib_bzmesh)

target_include_directories(lib_bzmesh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link with reuired libraries
target_link_libraries(lib_bzmesh PUBLIC libepp uepm::common)
target_link_libraries(lib_bzmesh PUBLIC Eigen3::Eigen rapidcsv fmt::fmt)
target_include_directories(lib_bzmesh PUBLIC ${GMSH_INC})
target_link_libraries(lib_bzmesh PUBLIC ${GMSH_LIB})
target_link_libraries(lib_bzmesh PUBLIC OpenMP::OpenMP_CXX)

# All users of this library will need at least C++20
target_compile_features(lib_bzmesh PUBLIC cxx_std_20)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(lib_bzmesh PRIVATE -Wall -Wextra -Wpedantic)
endif()
