set(HEADER_FILES_LIBFILE
  "finite_element.hpp" "finite_element2d.hpp" "finite_element3d.hpp"
  "poisson_solver_2d.hpp" "poisson_solver_3d.hpp"
  "nonlinear_poisson_solver.hpp")

set(SOURCE_FILES_LIBFILE
  "finite_element.cpp" "finite_element2d.cpp" "finite_element3d.cpp"
  "poisson_solver_2d.cpp" "poisson_solver_3d.cpp"
  "nonlinear_poisson_solver.cpp")

# Add the STF library for reading STF files
add_library(libpoisson STATIC ${SOURCE_FILES_LIBFILE} ${HEADER_FILES_LIBFILE})
add_library(uepm::poisson ALIAS libpoisson)

target_include_directories(libpoisson PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Add the STF library for reading STF files
target_link_libraries(libpoisson PRIVATE libmesh)
target_link_libraries(libpoisson PUBLIC pthread Eigen3::Eigen)
target_link_libraries(libpoisson PUBLIC libphysics)

if(USE_OPENMP_ACCELERATION)
  target_link_libraries(libmesh PUBLIC OpenMP::OpenMP_CXX)
endif(USE_OPENMP_ACCELERATION)
