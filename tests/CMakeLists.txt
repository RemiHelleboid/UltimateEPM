include(${doctest_SOURCE_DIR}/scripts/cmake/doctest.cmake OPTIONAL)

set(TEST_DIRS
  "${CMAKE_CURRENT_SOURCE_DIR}/bz_mesh"
  "${CMAKE_CURRENT_SOURCE_DIR}/epp"
  "${CMAKE_CURRENT_SOURCE_DIR}/fbmc"
  "${CMAKE_CURRENT_SOURCE_DIR}/poisson"
  "${CMAKE_CURRENT_SOURCE_DIR}/PBMC"
  "${CMAKE_CURRENT_SOURCE_DIR}/ADMC"
  "${CMAKE_CURRENT_SOURCE_DIR}/MMMC"
  "${CMAKE_CURRENT_SOURCE_DIR}/physics"
)

foreach(dir IN LISTS TEST_DIRS)
  file(GLOB TEST_SOURCES CONFIGURE_DEPENDS "${dir}/*.cpp")

  foreach(src IN LISTS TEST_SOURCES)
    get_filename_component(base "${src}" NAME_WE)
    get_filename_component(dir_name "${dir}" NAME)  # unit/integration/...
    # Unique target: e.g. unit_math, integration_db
    set(exe "${dir_name}_${base}")

    add_executable(${exe} "${src}")
    target_link_libraries(${exe} PRIVATE doctest::doctest libepp lib_bzmesh libpoisson uepm::mesh libfile ${EXTRA_LIBS})
    if(dir_name STREQUAL "PBMC")
      target_link_libraries(${exe} PRIVATE uepm::PBMC)
    endif()
    if(dir_name STREQUAL "ADMC")
      target_link_libraries(${exe} PRIVATE uepm::ADMC)
    endif()
    if(dir_name STREQUAL "MMMC")
      target_link_libraries(${exe} PRIVATE uepm::MMMC)
    endif()
    if(dir_name STREQUAL "fbmc")
      target_link_libraries(${exe} PRIVATE uepm::fbmc)
    endif()
    if(dir_name STREQUAL "physics")
      target_link_libraries(${exe} PRIVATE uepm::physics)
    endif()

    if(COMMAND doctest_discover_tests)
      doctest_discover_tests(${exe}
        DISCOVERY_MODE POST_BUILD
        TEST_PREFIX "${dir_name}::"
        PROPERTIES LABELS "${dir_name}"
      )
    else()
      add_test(NAME "${dir_name}_${base}" COMMAND ${exe})
      set_tests_properties("${dir_name}_${base}" PROPERTIES LABELS "${dir_name}")
    endif()
  endforeach()
endforeach()
