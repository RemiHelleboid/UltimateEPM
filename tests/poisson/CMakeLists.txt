add_executable(test_fem_2d test_fem_2d.cpp)
target_link_libraries(test_fem_2d PRIVATE libutils plog doctest::doctest libmesh libpoisson libfile)
add_test(test_fem_2d test_fem_2d)

add_executable(test_poisson_2d test_poisson_2d.cpp)
target_link_libraries(test_poisson_2d PUBLIC libutils plog doctest::doctest libmesh libpoisson libfile)
add_test(test_poisson_2d test_poisson_2d)

add_executable(test_poisson_3d test_poisson_3d.cpp)
target_link_libraries(test_poisson_3d PUBLIC libutils plog doctest::doctest libmesh libpoisson libfile)
add_test(test_poisson_3d test_poisson_3d)

add_executable(test_fem_3d test_fem_3d.cpp)
target_link_libraries(test_fem_3d PRIVATE libutils plog doctest::doctest libmesh libpoisson libfile)
add_test(test_fem_3d test_fem_3d)
