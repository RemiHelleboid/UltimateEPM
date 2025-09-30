cmake_minimum_required(VERSION 3.15)

include(CMakePrintHelpers)

project(
  EmpiricalPseudopotential
  VERSION "1.0.0"
  DESCRIPTION "Empirical Pseudo Potential Code"
  LANGUAGES C CXX)

add_compile_definitions(CMAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# -----------------------------------------------------------------------------
# Set the project options.
option(ENABLE_OPEN_MP "Build programs with suppot of OpenMP." ON)
option(ENABLE_MPI_BUILD "Build some programs with MPI parralelization." ON)
option(ENABLE_CLANG_TIDY "Enable clang tidy at compilation" OFF)
option(ENABLE_BUILD_DOC "Enable the documentation build (doxygen)" OFF)
option(ENABLE_BUILD_TEST "Enable the documentation build (doxygen)" ON)
option(ENABLE_MARCH_NATIVE "Enable optimization for the local machine." ON)

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

if (ENABLE_MARCH_NATIVE AND CMAKE_BUILD_TYPE STREQUAL "Release")
  message("Optimization for the local machine is enabled")
  add_compile_options(-march=native)
  add_compile_options(-O3)
  add_compile_options(-ffast-math -fno-trapping-math -fno-math-errno)
  add_compile_options(-funroll-loops)
  add_compile_options(-fomit-frame-pointer)
  add_compile_options(-fstrict-aliasing)
  add_compile_options(-fopenmp-simd)
  add_compile_options(-ftree-vectorize)
endif()


# Flags for the compiler
if(MSVC)
  add_compile_options("/W4" "$<$<CONFIG:RELEASE>:/O2>")
else()
  # add_compile_options("-Wall" "-Wextra" "-Werror" "-pedantic")

  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    add_compile_options("-stdlib=libc++")
  else()
    # nothing special for gcc at the moment
  endif()
endif()

# exclude Anaconda directories from search
if(DEFINED ENV{CONDA_PREFIX})
  set(h5_ignore_path
    $ENV{CONDA_PREFIX}/bin $ENV{CONDA_PREFIX}/lib $ENV{CONDA_PREFIX}/include
    $ENV{CONDA_PREFIX}/Library/bin $ENV{CONDA_PREFIX}/Library/lib $ENV{CONDA_PREFIX}/Library/include
  )
  list(APPEND CMAKE_IGNORE_PATH ${h5_ignore_path})
endif()


# -----------------------------------------------------------------------------
# FIND / GET RELEVANT LIBRARIES
include(FetchContent)


# Get yaml_cpp lib
message("Fetching yaml-cpp lib ...")
FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG master
)
FetchContent_MakeAvailable(yaml-cpp)


# Find Eigen library Find Eigen library
FetchContent_Declare(
  Eigen
  GIT_REPOSITORY "https://gitlab.com/libeigen/eigen.git"
  GIT_TAG "3.4.0"
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)

# find_package(Eigen3)

# if(NOT Eigen3_FOUND)
  set(EIGEN_BUILD_DOC OFF)
  set(BUILD_TESTING OFF)
  set(EIGEN_BUILD_PKGCONFIG OFF)
  message("Fetching Eigen3 lib ...")
  # add_compile_definitions(EIGEN_USE_LAPACK=0)
  FetchContent_MakeAvailable(Eigen)
  find_package(Eigen3 REQUIRED)
# endif()

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

find_path(GMSH_INC gmsh.h HINTS $ENV{GMSH_INC})

if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find gmsh.h")
endif()

message("---------------- GMSH.a FOUND" ${GMSH_LIB})
message("---------------- GMSH.H FOUND" ${GMSH_INC})

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
  DocTest
  GIT_REPOSITORY "https://github.com/onqtam/doctest"
  GIT_TAG "1da23a3e8119ec5cce4f9388e91b065e20bf06f5")

FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY "https://github.com/jbeder/yaml-cpp"
  GIT_TAG "master")

FetchContent_Declare(
  rapidcsv
  GIT_REPOSITORY "https://github.com/d99kris/rapidcsv.git"
  GIT_TAG "master"
  CONFIGURE_COMMAND "" BUILD_COMMAND "")

# Get fast-parse-csv lib
find_package(rapidcsv NO_MODULE)

if(NOT rapidcsv_FOUND)
  message("Fetching rapidcsv lib ...")
  FetchContent_MakeAvailable(rapidcsv)
  include_directories(${rapidcsv_SOURCE_DIR})
endif()

# Get yaml_cpp lib
# find_package(yaml-cpp NO_MODULE)
# FetchContent_MakeAvailable(yaml-cpp)


# if(NOT yaml-cpp_FOUND)
#   message("Fetching yaml-cpp lib ...")
#   set(YAML_CPP_BUILD_TESTS OFF)
#   set(YAML_CPP_BUILD_TOOLS OFF)
#   set(YAML_CPP_BUILD_CONTRIB OFF)
#   FetchContent_MakeAvailable(yaml-cpp)
# endif()

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



FetchContent_Declare(
  fmt
  GIT_REPOSITORY "https://github.com/fmtlib/fmt.git"
  GIT_TAG master)
FetchContent_MakeAvailable(fmt)



# MPI SUPPORT
if(ENABLE_MPI_BUILD)
  message("MPI support activated.")
  find_package(MPI)

  if(MPI_FOUND)
    set(USE_MPI_ACCELERATION 1)
    add_compile_definitions(USE_MPI=1)
    add_definitions(-DOMPI_SKIP_MPICXX)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
    message("MPI ACCELERATION ENABLED.")
  endif()
else()
  set(USE_MPI_ACCELERATION 0)
  add_compile_definitions(USE_MPI=0)
  message("MPI support deactivated. Not found.")
endif(ENABLE_MPI_BUILD)

# The tests are here
if(ENABLE_BUILD_TEST)
  message("Fetching DocTest lib ...")
  FetchContent_MakeAvailable(DocTest)
  add_subdirectory(tests/)
endif(ENABLE_BUILD_TEST)

# tclap library for argument parsing
include_directories(${PROJECT_SOURCE_DIR}/external/tclap-1.4.0/include)
# -----------------------------------------------------------------------------
# The compiled library code is here
add_subdirectory(src/EPP)
add_subdirectory(src/BZ_MESH)
add_subdirectory(src/FBMC)
# The applications are here
add_subdirectory(apps/)
