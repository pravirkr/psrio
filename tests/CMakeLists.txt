include(FetchContent)

# Catch2: prefer a system package (apt/brew) when it is at least 3.8; otherwise
# pin 3.8.0. Older releases are rejected: before 3.5.3 the assertion macros trip
# clang-tidy's bugprone-chained-comparison, and before 3.8 catch_discover_tests
# reports a SKIP as a failure. Ubuntu 24.04 ships 3.4.0.
find_package(
  Catch2
  3.8
  CONFIG
  QUIET
)
if(Catch2_FOUND)
  message(STATUS "Using system Catch2 ${Catch2_VERSION}")
else()
  if(NOT DEFINED PSRIO_CATCH2_VERSION)
    set(PSRIO_CATCH2_VERSION "3.8.0")
  endif()
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v${PSRIO_CATCH2_VERSION}
    GIT_SHALLOW TRUE
  )
  set(CATCH_INSTALL_DOCS
      OFF
      CACHE BOOL "" FORCE
  )
  set(CATCH_INSTALL_EXTRAS
      OFF
      CACHE BOOL "" FORCE
  )
  FetchContent_MakeAvailable(Catch2)
endif()

if(NOT TARGET Catch2::Catch2WithMain)
  message(
    FATAL_ERROR "Catch2 is required for PSRIO_BUILD_TESTS=ON but was not found or downloaded."
  )
endif()

list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)

set(TARGET_TESTS ${PROJECT_NAME}_tests)

file(
  GLOB
  TEST_SOURCES
  CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)
# Opt-in formats and the CLI compile only in their own test binaries.
list(FILTER TEST_SOURCES EXCLUDE REGEX
     "test_fbh5\\.cpp$|test_psrdada\\.cpp$|test_cli_commands\\.cpp$"
)

add_executable(${TARGET_TESTS} ${TEST_SOURCES})

target_link_libraries(
  ${TARGET_TESTS} PRIVATE ${PROJECT_NAME}::${PROJECT_NAME} Catch2::Catch2WithMain
                          psrio_developer_flags
)

target_compile_definitions(
  ${TARGET_TESTS}
  PRIVATE "PSRIO_TEST_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/data\""
          "PSRIO_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\""
)

include(Catch)
catch_discover_tests(${TARGET_TESTS})

# Standalone header check: compile every public header on its own as a one-line
# translation unit so missing self-includes fail the build immediately.
file(
  GLOB_RECURSE
  PSRIO_PUBLIC_HEADERS
  CONFIGURE_DEPENDS
  RELATIVE ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/include/psrio/*.hpp
)
list(FILTER PSRIO_PUBLIC_HEADERS EXCLUDE REGEX
     "^psrio/formats/fbh5\\.hpp$|^psrio/formats/psrdada\\.hpp$"
)

set(PSRIO_HEADER_CHECK_DIR "${CMAKE_CURRENT_BINARY_DIR}/header_checks")
set(PSRIO_HEADER_CHECK_SOURCES "")
foreach(header IN LISTS PSRIO_PUBLIC_HEADERS)
  string(MAKE_C_IDENTIFIER "${header}" _header_id)
  set(_check_src "${PSRIO_HEADER_CHECK_DIR}/${_header_id}.cpp")
  file(
    GENERATE
    OUTPUT "${_check_src}"
    CONTENT "#include \"${header}\" // IWYU pragma: associated\n"
  )
  list(APPEND PSRIO_HEADER_CHECK_SOURCES "${_check_src}")
endforeach()

add_library(psrio_header_checks OBJECT ${PSRIO_HEADER_CHECK_SOURCES})
target_link_libraries(
  psrio_header_checks PRIVATE ${PROJECT_NAME}::${PROJECT_NAME} psrio_developer_flags
)
add_dependencies(${TARGET_TESTS} psrio_header_checks)

if(PSRIO_WITH_HDF5)
  set(_fbh5_check_src "${PSRIO_HEADER_CHECK_DIR}/psrio_formats_fbh5_hpp.cpp")
  file(
    GENERATE
    OUTPUT "${_fbh5_check_src}"
    CONTENT "#include \"psrio/formats/fbh5.hpp\" // IWYU pragma: associated\n"
  )
  add_executable(
    psrio_fbh5_tests ${_fbh5_check_src} ${CMAKE_CURRENT_SOURCE_DIR}/test_fbh5.cpp
  )
  target_link_libraries(
    psrio_fbh5_tests PRIVATE psrio::fbh5 Catch2::Catch2WithMain psrio_developer_flags
  )
  target_compile_definitions(
    psrio_fbh5_tests PRIVATE "PSRIO_TEST_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/data\""
  )
  catch_discover_tests(
    psrio_fbh5_tests
    PROPERTIES
    ENVIRONMENT "HDF5_DISABLE_VERSION_CHECK=2"
  )
endif()

if(TARGET psrio_cli_commands)
  add_executable(psrio_cli_tests ${CMAKE_CURRENT_SOURCE_DIR}/test_cli_commands.cpp)
  target_link_libraries(
    psrio_cli_tests PRIVATE psrio_cli_commands dmt::dmt Catch2::Catch2WithMain
                            psrio_developer_flags
  )
  target_compile_definitions(
    psrio_cli_tests PRIVATE "PSRIO_TEST_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/data\""
  )
  catch_discover_tests(
    psrio_cli_tests
    PROPERTIES
    ENVIRONMENT "HDF5_DISABLE_VERSION_CHECK=2"
  )
endif()

if(PSRIO_WITH_PSRDADA)
  set(_psrdada_check_src "${PSRIO_HEADER_CHECK_DIR}/psrio_formats_psrdada_hpp.cpp")
  file(
    GENERATE
    OUTPUT "${_psrdada_check_src}"
    CONTENT "#include \"psrio/formats/psrdada.hpp\" // IWYU pragma: associated\n"
  )
  add_executable(
    psrio_psrdada_tests ${_psrdada_check_src} ${CMAKE_CURRENT_SOURCE_DIR}/test_psrdada.cpp
  )
  target_link_libraries(
    psrio_psrdada_tests PRIVATE psrio::psrdada Catch2::Catch2WithMain psrio_developer_flags
  )
  catch_discover_tests(psrio_psrdada_tests)
endif()
