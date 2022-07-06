cmake_minimum_required(VERSION 3.15...3.25)

project(
  EmpiricalPseudopotential
  VERSION "1.0.0"
  DESCRIPTION "Empirical Pseudo Potential Code"
  LANGUAGES C CXX)

# -----------------------------------------------------------------------------
# Set the project options.

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

# -----------------------------------------------------------------------------
# FIND / GET RELEVANT LIBRARIES
include(FetchContent)
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

# -----------------------------------------------------------------------------
# The compiled library code is here
add_subdirectory(EmpiricalPseudopotential/)
