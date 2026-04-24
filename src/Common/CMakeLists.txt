# Header-only
add_library(uepm_common INTERFACE)
add_library(uepm::common ALIAS uepm_common)

target_include_directories(uepm_common INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>)

target_link_libraries(uepm_common INTERFACE rapidcsv fmt::fmt)


target_compile_features(uepm_common INTERFACE cxx_std_20)

