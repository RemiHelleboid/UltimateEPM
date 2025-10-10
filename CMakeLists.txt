cmake_minimum_required(VERSION 3.15)

project(
  EmpiricalPseudopotential
  VERSION 1.0.0
  DESCRIPTION "Empirical Pseudo Potential Code"
  LANGUAGES C CXX)

include(CMakePrintHelpers)
include(FetchContent)
include(CTest)
include(cmake/CPM.cmake)

# ---------------------------- Options ----------------------------
option(ENABLE_MPI_BUILD     "Build with MPI"                               ON)
option(ENABLE_CLANG_TIDY    "Enable clang-tidy"                            OFF)
option(ENABLE_BUILD_DOC     "Build documentation (Doxygen)"                OFF)
option(ENABLE_BUILD_TEST    "Build tests (doctest/CTest)"                  ON)
option(ENABLE_MARCH_NATIVE  "Enable -march=native and aggressive flags"    ON)

# ---------------------------- Toolchain / Lang ----------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_compile_definitions(PROJECT_SRC_DIR=\"${PROJECT_SOURCE_DIR}\")

# ---------------------------- Build type default ----------------------------
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(default_build_type "RelWithDebInfo")
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}'")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# ---------------------------- Clang-Tidy ----------------------------
if(ENABLE_CLANG_TIDY)
  message(STATUS "clang-tidy enabled")
  set(CMAKE_CXX_CLANG_TIDY clang-tidy -checks=-*,readability-*)
endif()

# ---------------------------- CPU/Vectorization flags ----------------------------
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-march=native -O3 -fno-trapping-math -fno-math-errno
                         -fomit-frame-pointer -fstrict-aliasing -fopenmp-simd
                        -ftree-vectorize)
  endif()
endif()

# ---------------------------- Platform warnings ----------------------------
if(MSVC)
  add_compile_options(/W4 $<$<CONFIG:RELEASE>:/O2>)
else()
  # add_compile_options(-Wall -Wextra -Werror -pedantic)  # uncomment if desired
endif()

# ---------------------------- Conda ignore (optional) ----------------------------
if(DEFINED ENV{CONDA_PREFIX})
  list(APPEND CMAKE_IGNORE_PATH
       $ENV{CONDA_PREFIX}/bin $ENV{CONDA_PREFIX}/lib $ENV{CONDA_PREFIX}/include
       $ENV{CONDA_PREFIX}/Library/bin $ENV{CONDA_PREFIX}/Library/lib $ENV{CONDA_PREFIX}/Library/include)
endif()

# ---------------------------- Dependencies (FetchContent) ----------------------------
# yaml-cpp (pin release)
message(STATUS "Fetching yaml-cpp")


FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG 65c1c270dbe7eec37b2df2531d7497c4eea79aee
)

FetchContent_MakeAvailable(yaml-cpp)

# Eigen (pin release)
CPMAddPackage(
  NAME Eigen
  VERSION 5.0.0
  URL https://gitlab.com/libeigen/eigen/-/archive/5.0.0/eigen-5.0.0.tar.gz
  # Eigen's CMakelists are not intended for library use
  DOWNLOAD_ONLY YES 
)

if(Eigen_ADDED)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  target_include_directories(Eigen3::Eigen INTERFACE ${Eigen_SOURCE_DIR})
endif()



# fmt (pin release)
message(STATUS "Fetching fmt")
FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 11.0.2
)
FetchContent_MakeAvailable(fmt)

# rapidcsv (header-only)
message(STATUS "Fetching rapidcsv")
FetchContent_Declare(
  rapidcsv
  GIT_REPOSITORY https://github.com/d99kris/rapidcsv.git
  GIT_TAG master
  CONFIGURE_COMMAND "" BUILD_COMMAND ""
)
FetchContent_MakeAvailable(rapidcsv)
include_directories(${rapidcsv_SOURCE_DIR})

# doctest (for tests)
if(ENABLE_BUILD_TEST AND BUILD_TESTING)
  message(STATUS "Fetching doctest")
  FetchContent_Declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG v2.4.12
  )
  FetchContent_MakeAvailable(doctest)
endif()

# ---------------------------- OpenMP ----------------------------
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
  set(OPENMP_IMPORTED_TARGET OpenMP::OpenMP_CXX CACHE STRING "" FORCE)
    message(STATUS "OpenMP found")
  else()
    message(STATUS "OpenMP not found")
endif()

# ---------------------------- MPI ----------------------------
if(ENABLE_MPI_BUILD)
  find_package(MPI)
  if(MPI_CXX_FOUND)
    add_compile_definitions(USE_MPI=1)
    add_definitions(-DOMPI_SKIP_MPICXX)
    include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
    set(MPI_IMPORTED_TARGET MPI::MPI_CXX CACHE STRING "" FORCE)
    message(STATUS "MPI found")
  else()
    add_compile_definitions(USE_MPI=0)
    message(STATUS "MPI not found")
  endif()
else()
  add_compile_definitions(USE_MPI=0)
endif()

# ---------------------------- GMSH ----------------------------
find_library(GMSH_LIB gmsh HINTS $ENV{GMSH_LIB})
if(NOT GMSH_LIB)
  message(FATAL_ERROR "Could not find libgmsh")
endif()
find_path(GMSH_INC gmsh.h HINTS $ENV{GMSH_INC})
if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find gmsh.h")
endif()
cmake_print_variables(GMSH_LIB GMSH_INC)


# ---------------------------- BLAS / LAPACK ----------------------------
find_package(BLAS)
find_package(LAPACK)
set(EXTRA_LIBS "")
if(BLAS_FOUND)
  list(APPEND EXTRA_LIBS ${BLAS_LIBRARIES})
  message(STATUS "BLAS: ${BLAS_LIBRARIES}")
endif()
if(LAPACK_FOUND)
  list(APPEND EXTRA_LIBS ${LAPACK_LIBRARIES})
  message(STATUS "LAPACK: ${LAPACK_LIBRARIES}")
endif()

# ---------------------------- Doxygen ----------------------------
find_package(Doxygen)
if(DOXYGEN_FOUND AND ENABLE_BUILD_DOC)
  set(DOXYGEN_IN  ${PROJECT_SOURCE_DIR}/doc/Doxyfile)
  set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  set(MAINPAGE_IN ${PROJECT_SOURCE_DIR}/doc/mainpage.md)
  set(MAINPAGE_OUT ${CMAKE_CURRENT_BINARY_DIR}/mainpage.md)
  configure_file(${MAINPAGE_IN} ${MAINPAGE_OUT} @ONLY)
  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
  add_custom_target(doc_doxygen ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
else()
  if(ENABLE_BUILD_DOC)
    message(STATUS "Doxygen not found; docs will not be built")
  endif()
endif()

# ---------------------------- Third-party include shortcuts ----------------------------
include_directories(${GMSH_INC})
# tclap (header-only)
include_directories(${PROJECT_SOURCE_DIR}/external/tclap-1.4.0/include)

# ---------------------------- Subprojects ----------------------------
add_subdirectory(src/EPP)
add_subdirectory(src/BZ_MESH)
add_subdirectory(src/FBMC)
add_subdirectory(apps)

# ---------------------------- Tests ----------------------------
enable_testing()
if(ENABLE_BUILD_TEST AND BUILD_TESTING)
  add_subdirectory(tests)
endif()
