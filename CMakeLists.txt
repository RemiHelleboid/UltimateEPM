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
option(ENABLE_BUILD_DOC "Enable the documentation build (doxygen)" OFF)

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
# Find Eigen library Find Eigen library
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
  add_compile_definitions(EIGEN_USE_LAPACK=0)
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

# FIND GMSH Adds GMSH
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

FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY "https://github.com/jbeder/yaml-cpp"
  GIT_TAG "328d2d85e833be7cb5a0ab246cc3f5d7e16fc67a")

# Get yaml_cpp lib
find_package(yaml-cpp NO_MODULE)
if(NOT yaml-cpp_FOUND)
  message("Fetching yaml-cpp lib ...")
  FetchContent_MakeAvailable(yaml-cpp)
endif()

# find python libraries
find_package(
  Python3
  COMPONENTS Interpreter Development NumPy
  REQUIRED)
find_package(PythonLibs 3.0 REQUIRED)
include_directories(${PYTHON3_INCLUDE_DIRS} ${NumPy_INCLUDE_DIRS})

# populate matplotlib repository
include(FetchContent)
FetchContent_Declare(
  matplotlib
  GIT_REPOSITORY https://github.com/lava/matplotlib-cpp.git
  GIT_TAG f23347fca25219d1c42cbb91608b5556814bf572)
FetchContent_GetProperties(matplotlib)
if(NOT matplotlib_POPULATED)
  FetchContent_Populate(matplotlib)
endif()
include_directories(SYSTEM ${matplotlib_SOURCE_DIR})

# check if Doxygen is installed
find_package(Doxygen)
if(DOXYGEN_FOUND AND ENABLE_BUILD_DOC)
  # set input and output files
  set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/doc/Doxyfile)
  set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  set(MAINPAGE_IN ${PROJECT_SOURCE_DIR}/doc/mainpage.md)
  set(MAINPAGE_OUT ${CMAKE_CURRENT_BINARY_DIR}/mainpage.md)

  # request to configure the file
  configure_file(${MAINPAGE_IN} ${MAINPAGE_OUT} @ONLY)
  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
  message("Doxygen build started")

  # note the option ALL which allows to build the docs together with the
  # application
  add_custom_target(
    doc_doxygen ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)
else(DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif()

# tclap library for argument parsing
include_directories(${PROJECT_SOURCE_DIR}/external/tclap/include)

# -----------------------------------------------------------------------------
# The compiled library code is here
add_subdirectory(src/EPP)
add_subdirectory(src/BZ_MESH)
add_subdirectory(apps/)
