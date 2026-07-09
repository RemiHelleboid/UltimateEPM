cmake_minimum_required(VERSION 3.15)

project(
  EmpiricalPseudopotential
  VERSION 1.0.0
  DESCRIPTION "Empirical Pseudopotential Code"
  LANGUAGES C CXX)

include(CMakePrintHelpers)
include(FetchContent)
include(CTest)
include(cmake/CPM.cmake)
include(cmake/FindGmsh.cmake)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# ---------------------------- Options ----------------------------
option(ENABLE_MPI_BUILD     "Build with MPI"                            ON)
option(ENABLE_CLANG_TIDY    "Enable clang-tidy"                         OFF)
option(ENABLE_BUILD_DOC     "Build documentation (Doxygen)"             OFF)
option(ENABLE_BUILD_TEST    "Build tests (doctest/CTest)"               OFF)
option(ENABLE_MARCH_NATIVE  "Enable -march=native / native arch flags"  ON)
option(ENABLE_FAST_MATH     "Enable unsafe fast-math optimizations"     OFF)
option(IGNORE_CONDA_PREFIX  "Ignore active Conda environment paths"     ON)

# ---------------------------- Toolchain / Language ----------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_compile_definitions(PROJECT_SRC_DIR=\"${PROJECT_SOURCE_DIR}\")

# ---------------------------- Project data ----------------------------
set(UEPM_MATERIAL_DATA_DIR
    "${PROJECT_SOURCE_DIR}/data/materials"
    CACHE PATH "Root directory of the UltimateEPM material repository")

if(NOT IS_DIRECTORY "${UEPM_MATERIAL_DATA_DIR}")
  message(FATAL_ERROR
          "UEPM material repository does not exist: ${UEPM_MATERIAL_DATA_DIR}")
endif()

message(STATUS "Material repository: ${UEPM_MATERIAL_DATA_DIR}")

# ---------------------------- Build type default ----------------------------
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(default_build_type "RelWithDebInfo")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}'")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
               "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ---------------------------- Clang-Tidy ----------------------------
if(ENABLE_CLANG_TIDY)
  find_program(CLANG_TIDY_EXE NAMES clang-tidy)
  if(CLANG_TIDY_EXE)
    message(STATUS "clang-tidy enabled: ${CLANG_TIDY_EXE}")
    set(CMAKE_CXX_CLANG_TIDY
        "${CLANG_TIDY_EXE};-checks=-*,readability-*,bugprone-*,performance-*")
  else()
    message(WARNING "ENABLE_CLANG_TIDY=ON but clang-tidy was not found")
  endif()
endif()

# ---------------------------- CPU / optimization flags ----------------------------
#
# Keep this global for now so existing subdirectory CMake files do not need changes.
# Long term, move this to a project_options INTERFACE target.
#
if(MSVC)
  add_compile_options(/W4)

  if(ENABLE_MARCH_NATIVE)
    add_compile_options(
      $<$<CONFIG:Release>:/arch:AVX2>
      $<$<CONFIG:RelWithDebInfo>:/arch:AVX2>)
  endif()

else()
  # Warnings are intentionally not -Werror here, so external/compiler-version noise
  # does not break development builds.
  add_compile_options(-Wall -Wextra)

  add_compile_options(
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:RelWithDebInfo>:-O3>)

  if(ENABLE_MARCH_NATIVE)
    message(STATUS "Enabling native architecture flags for Release/RelWithDebInfo")
    add_compile_options(
      $<$<CONFIG:Release>:-march=native>
      $<$<CONFIG:RelWithDebInfo>:-march=native>)
  endif()

  if(ENABLE_FAST_MATH)
    message(WARNING "ENABLE_FAST_MATH=ON: this can change floating-point behavior")
    add_compile_options(
      $<$<CONFIG:Release>:-ffast-math>
      $<$<CONFIG:RelWithDebInfo>:-ffast-math>)
  endif()
endif()

# ---------------------------- Conda ignore ----------------------------
if(IGNORE_CONDA_PREFIX AND DEFINED ENV{CONDA_PREFIX})
  message(STATUS "Ignoring active Conda prefix: $ENV{CONDA_PREFIX}")
  list(APPEND CMAKE_IGNORE_PATH
       "$ENV{CONDA_PREFIX}/bin"
       "$ENV{CONDA_PREFIX}/lib"
       "$ENV{CONDA_PREFIX}/include"
       "$ENV{CONDA_PREFIX}/Library/bin"
       "$ENV{CONDA_PREFIX}/Library/lib"
       "$ENV{CONDA_PREFIX}/Library/include")
endif()

# ---------------------------- Dependencies ----------------------------

# ---------------------------- Dependencies ----------------------------

CPMAddPackage("gh:jbeder/yaml-cpp#yaml-cpp-0.9.0@0.9.0")

CPMAddPackage(
  NAME Eigen
  VERSION 5.0.0
  URL https://gitlab.com/libeigen/eigen/-/archive/5.0.0/eigen-5.0.0.tar.gz
  DOWNLOAD_ONLY YES)

if(Eigen_ADDED AND NOT TARGET Eigen3::Eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  target_include_directories(Eigen3::Eigen SYSTEM INTERFACE "${Eigen_SOURCE_DIR}")
endif()

CPMAddPackage("gh:fmtlib/fmt#12.1.0@12.1.0")
CPMAddPackage("gh:d99kris/rapidcsv#v8.99@v8.99")

if(ENABLE_BUILD_TEST AND BUILD_TESTING)
  CPMAddPackage("gh:doctest/doctest#v2.5.2@v2.5.2")
endif()

find_package(OpenMP REQUIRED COMPONENTS CXX)
set(OPENMP_IMPORTED_TARGET OpenMP::OpenMP_CXX CACHE STRING "" FORCE)
set(USE_OPENMP_ACCELERATION ON CACHE BOOL "Enable OpenMP acceleration" FORCE)


# ---------------------------- MPI ----------------------------
if(ENABLE_MPI_BUILD)
  find_package(MPI)

  if(MPI_CXX_FOUND)
    add_compile_definitions(USE_MPI=1 OMPI_SKIP_MPICXX)
    include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})

    # Preserve the variable used by existing subdirectory CMake files.
    set(MPI_IMPORTED_TARGET MPI::MPI_CXX CACHE STRING "" FORCE)

    message(STATUS "MPI found")
  else()
    add_compile_definitions(USE_MPI=0)
    set(ENABLE_MPI_BUILD OFF CACHE BOOL "" FORCE)
    message(STATUS "MPI not found; building without MPI")
  endif()
else()
  add_compile_definitions(USE_MPI=0)
  message(STATUS "MPI disabled")
endif()

# ---------------------------- GMSH ----------------------------
find_package(Gmsh REQUIRED)
set(GMSH_LIB "${GMSH_LIBRARY}" CACHE FILEPATH "Path to libgmsh" FORCE)
set(GMSH_INC "${GMSH_INCLUDE_DIR}" CACHE PATH "Path to gmsh include directory" FORCE)
include_directories(SYSTEM "${GMSH_INC}")


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
if(ENABLE_BUILD_DOC)
  find_package(Doxygen)

  if(DOXYGEN_FOUND)
    set(DOXYGEN_IN  "${PROJECT_SOURCE_DIR}/doc/Doxyfile")
    set(DOXYGEN_OUT "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")
    set(MAINPAGE_IN "${PROJECT_SOURCE_DIR}/doc/mainpage.md")
    set(MAINPAGE_OUT "${CMAKE_CURRENT_BINARY_DIR}/mainpage.md")

    configure_file("${MAINPAGE_IN}" "${MAINPAGE_OUT}" @ONLY)
    configure_file("${DOXYGEN_IN}" "${DOXYGEN_OUT}" @ONLY)

    add_custom_target(doc_doxygen ALL
      COMMAND "${DOXYGEN_EXECUTABLE}" "${DOXYGEN_OUT}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      COMMENT "Generating API documentation with Doxygen"
      VERBATIM)
  else()
    message(STATUS "Doxygen not found; docs will not be built")
  endif()
endif()

# ---------------------------- Third-party include shortcuts ----------------------------
#
# These are intentionally kept global for now so existing subdirectory CMake files
# do not need changes. Long term, replace with target_include_directories().
#
include_directories("${PROJECT_SOURCE_DIR}/external/tclap-1.4.0/include")
include_directories("${PROJECT_SOURCE_DIR}/external/fast-cpp-csv-parser")
include_directories("${rapidcsv_SOURCE_DIR}")

# ---------------------------- Subprojects ----------------------------
add_subdirectory(src/Common)
add_subdirectory(src/Physics)
add_subdirectory(src/EPP)
add_subdirectory(src/Parser)
add_subdirectory(src/BZ_MESH)
add_subdirectory(src/FBMC)
add_subdirectory(src/PBMC)
add_subdirectory(src/ADMC)
add_subdirectory(src/DEVICE_MESH)
add_subdirectory(src/Poisson)
add_subdirectory(src/MMMC)
add_subdirectory(apps)

# ---------------------------- Tests ----------------------------
if(ENABLE_BUILD_TEST AND BUILD_TESTING)
  add_subdirectory(tests)
endif()
