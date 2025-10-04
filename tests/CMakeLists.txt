# include(${doctest_SOURCE_DIR}/scripts/cmake/doctest.cmake OPTIONAL)

file(GLOB TEST_SOURCES *.cpp)
message("TEST_SOURCES: ${TEST_SOURCES}")

message(STATUS "TESTS:")
foreach(TEST_SRC ${TEST_SOURCES})
  get_filename_component(TEST_NAME ${TEST_SRC} NAME_WE)
  message(STATUS "   └─ ${TEST_NAME}")

  add_executable(${TEST_NAME} ${TEST_SRC})
  target_link_libraries(${TEST_NAME} PRIVATE
    doctest::doctest
    libepp
    lib_bzmesh
    ${EXTRA_LIBS}      # BLAS/LAPACK etc., if you collected them
  )

  if(COMMAND doctest_discover_tests)
    doctest_discover_tests(${TEST_NAME} DISCOVERY_MODE POST_BUILD)
  else()
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
  endif()
endforeach()
