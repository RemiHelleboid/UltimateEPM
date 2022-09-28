file(GLOB BENCHMARKS *.cpp)

message("TESTS: ")

foreach(BENCH ${BENCHMARKS})
    get_filename_component(EXAMPLENAME ${BENCH} NAME_WE)
    message("   --- Test: ${EXAMPLENAME}")
    add_executable(${EXAMPLENAME} ${BENCH})
    target_link_libraries(${EXAMPLENAME} doctest::doctest libepp lib_bzmesh)
    add_test(${EXAMPLENAME} ${EXAMPLENAME})
endforeach()