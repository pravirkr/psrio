if(NOT TARGET CLI11::CLI11)
  CPMAddPackage(
    NAME CLI11
    VERSION 2.6.2
    URL https://github.com/CLIUtils/CLI11/archive/refs/tags/v2.6.2.tar.gz
    URL_HASH SHA256=c6ea6b2e5608b3ea8617999bd5f47420c71b2ebdb8dc4767c1034d1da5785711
    DOWNLOAD_ONLY YES
    EXCLUDE_FROM_ALL YES
    SYSTEM YES
  )
  if(CLI11_ADDED)
    add_library(CLI11::CLI11 INTERFACE IMPORTED GLOBAL)
    target_include_directories(CLI11::CLI11 SYSTEM INTERFACE "${CLI11_SOURCE_DIR}/include")
  endif()
endif()
if(NOT TARGET CLI11::CLI11)
  message(FATAL_ERROR "CLI11::CLI11 target not available after CPM setup.")
endif()

# dmt needs OpenMP and FFTW. Pass -DCPM_dmt_SOURCE=<path> to build a local
# checkout instead of the pinned commit.
if(NOT TARGET dmt::dmt)
  CPMAddPackage(
    NAME dmt
    GIT_REPOSITORY https://github.com/pravirkr/dmt.git
    GIT_TAG 64fc478835586ae8b91b4b96f5edf9eba021b9fa
    OPTIONS "DMT_BUILD_PYTHON OFF" "DMT_BUILD_TESTING OFF" "DMT_BUILD_BENCHMARKS OFF"
            "DMT_BUILD_DOCS OFF"
    EXCLUDE_FROM_ALL YES
    SYSTEM YES
  )
endif()
if(NOT TARGET dmt::dmt)
  message(FATAL_ERROR "dmt::dmt target not available after CPM setup.")
endif()

# The commands live in a library so the tests can drive them in-process.
add_library(
  psrio_cli_commands STATIC
  cli.cpp
  common.cpp
  cmd_inspect.cpp
  cmd_slice.cpp
  cmd_splice.cpp
  cmd_decimate.cpp
  cmd_convert.cpp
  cmd_digifil.cpp
  cmd_dedisperse.cpp
)
target_include_directories(psrio_cli_commands PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(
  psrio_cli_commands
  PUBLIC psrio::psrio CLI11::CLI11
  PRIVATE dmt::dmt psrio_developer_flags
)
target_compile_definitions(psrio_cli_commands PRIVATE PSRIO_CLI_VERSION="${PROJECT_VERSION}")
if(TARGET psrio::fbh5)
  target_link_libraries(psrio_cli_commands PUBLIC psrio::fbh5)
  target_compile_definitions(psrio_cli_commands PRIVATE PSRIO_CLI_HAS_HDF5=1)
endif()

add_executable(psrio_cli psrio.cpp)
set_target_properties(psrio_cli PROPERTIES OUTPUT_NAME psrio)
target_link_libraries(psrio_cli PRIVATE psrio_cli_commands psrio_developer_flags)

if(PSRIO_INSTALL)
  install(TARGETS psrio_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
