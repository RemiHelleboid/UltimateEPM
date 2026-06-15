# cmake/FindGmsh.cmake
#
# Defines:
#   Gmsh::Gmsh
#   GMSH_FOUND
#   GMSH_INCLUDE_DIR
#   GMSH_LIBRARY
#
# Also supports:
#   GMSH_DIR
#   GMSH_INC
#   GMSH_LIB

include(FindPackageHandleStandardArgs)

find_path(GMSH_INCLUDE_DIR
  NAMES gmsh.h
  HINTS
    ENV GMSH_INC
    ENV GMSH_DIR
  PATH_SUFFIXES
    include
)

find_library(GMSH_LIBRARY
  NAMES gmsh
  HINTS
    ENV GMSH_LIB
    ENV GMSH_DIR
  PATH_SUFFIXES
    lib
    lib64
)

find_package_handle_standard_args(Gmsh
  REQUIRED_VARS
    GMSH_LIBRARY
    GMSH_INCLUDE_DIR
)

if(GMSH_FOUND AND NOT TARGET Gmsh::Gmsh)
  add_library(Gmsh::Gmsh UNKNOWN IMPORTED)
  set_target_properties(Gmsh::Gmsh PROPERTIES
    IMPORTED_LOCATION "${GMSH_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${GMSH_INCLUDE_DIR}"
  )
endif()

mark_as_advanced(GMSH_INCLUDE_DIR GMSH_LIBRARY)