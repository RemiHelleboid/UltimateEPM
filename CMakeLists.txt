cmake_minimum_required(VERSION 3.15...3.25)

include(CMakePrintHelpers)


project(
  EmpiricalPseudopotential
  VERSION "1.0.0"
  DESCRIPTION "Empirical Pseudo Potential Code"
  LANGUAGES C CXX)

add_compile_definitions(CMAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# -----------------------------------------------------------------------------
# Set the project options.

option(ENABLE_OPEN_MP "Build test programs" ON)
option(ENABLE_CLANG_TIDY "Enable clang tidy at compilation" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Handling build type nicely Set a default build type if none was specified
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(default_build_type "RelWithDebInfo")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(
    STATUS
      "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "${default_build_type}"
      CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

if(ENABLE_CLANG_TIDY)
  message("Clang tidy at compilation is enabled")
  set(CMAKE_CXX_CLANG_TIDY clang-tidy -checks=-*,readability-*)
endif()

# -----------------------------------------------------------------------------
# FIND / GET RELEVANT LIBRARIES
include(FetchContent)
# Find Eigen library
# Find Eigen library
FetchContent_Declare(
  Eigen
  GIT_REPOSITORY "https://gitlab.com/libeigen/eigen.git"
  GIT_TAG "3.4.0"
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)

find_package(Eigen3 NO_MODULE)
if(NOT Eigen3_FOUND)
  set(EIGEN_BUILD_DOC OFF)
  set(BUILD_TESTING OFF)
  set(EIGEN_BUILD_PKGCONFIG OFF)
  message("Fetching Eigen3 lib ...")
  add_compile_definitions(EIGEN_USE_LAPACKE=0)
  FetchContent_MakeAvailable(Eigen)
endif()

# Find OpenMP
if(ENABLE_OPEN_MP)
  message("OpenMP optimization activated.")
  find_package(OpenMP)
  if(OpenMP_FOUND)
    set(USE_OPENMP_ACCELERATION 1)
    add_compile_definitions(USE_OPEN_MP=1)
  endif()
else()
  set(USE_OPENMP_ACCELERATION 0)
  add_compile_definitions(USE_OPEN_MP=0)
endif(ENABLE_OPEN_MP)

# FIND GMSH
# Adds GMSH
find_library(GMSH_LIB gmsh HINTS $ENV{GMSH_LIB})
if(NOT GMSH_LIB)
  message(FATAL_ERROR "Could not find libgmsh")
endif()
message("GMSH.a FOUND" ${GMSH_LIB})

find_path(GMSH_INC gmsh.h HINTS $ENV{GMSH_INC})
if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find gmsh.h")
endif()
message("GMSH.H FOUND" ${GMSH_INC})

cmake_print_variables(GMSH_LIB)
cmake_print_variables(GMSH_INC)

if(GMSH_LIB MATCHES ".a") # FIXME - generalize this
  find_library(BLAS_LIB blas)
  if(BLAS_LIB)
    list(APPEND EXTRA_LIBS ${BLAS_LIB})
    message("BLAS FOUND : " ${BLAS_LIB})
  endif()
  find_library(LAPACK_LIB lapack)
  if(LAPACK_LIB)
    list(APPEND EXTRA_LIBS ${LAPACK_LIB})
    message("LAPACK FOUND : " ${BLAS_LIB})
  endif()
endif()


# tclap library for argument parsing
include_directories(${PROJECT_SOURCE_DIR}/external/tclap/include)

# -----------------------------------------------------------------------------
# The compiled library code is here
add_subdirectory(src/EPP)
add_subdirectory(src/BZ_MESH)
add_subdirectory(apps/)
