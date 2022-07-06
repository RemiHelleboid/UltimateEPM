set(HEADER_FILES_LIBEPP
    SymmetryPoints.h
    Options.h
    Material.h
    Hamiltonian.h
    Vector3D.h
    Vector3D.inl
    BandStructure.h
    Pseudopotential.h)

set(SOURCE_FILES_LIBEPP
    Material.cpp
    SymmetryPoints.cpp
    Options.cpp
    Hamiltonian.cpp
    BandStructure.cpp
    Pseudopotential.cpp
    )

add_library(libepp STATIC ${SOURCE_FILES_LIBEPP} ${HEADER_FILES_LIBEPP})

target_include_directories(libepp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(libepp PUBLIC pthread Eigen3::Eigen)
# target_link_libraries(libepp PUBLIC ${wxWidgets_LIBRARIES})
# target_link_libraries(libepp PUBLIC ${VTK_LIBRARIES})
# target_link_libraries(libepp PUBLIC ${GTK_LIBRARIES})
# target_link_libraries(libepp PUBLIC ${GTK_USE_FILE})

# All users of this library will need at least C++20
target_compile_features(libepp PRIVATE cxx_std_20)

# add_executable(EPseudopotentialApp EPseudopotentialApp.cpp)
# target_link_libraries(EPseudopotentialApp PUBLIC libepp pthread ${wxWidgets_LIBRARIES}
#                                                  ${VTK_LIBRARIES})
# target_compile_features(EPseudopotentialApp PRIVATE cxx_std_20)


add_executable(EmpiricalPseudoPotentialMain EmpiricalPseudoPotentialMain.cpp)
target_link_libraries(EmpiricalPseudoPotentialMain PUBLIC libepp pthread )
